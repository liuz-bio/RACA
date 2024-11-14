#include <iostream>
#include <fstream>
#include <string>
#include <coroutine>
#include <optional>
#include <string_view>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <future>        
#include <functional>    
#include <filesystem>
#include <memory>      
#include <algorithm>   
#include <vector>
#include <mutex>
#include <sstream>
#include <cmath>
#include <unordered_set>
#include <set>
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/tbx.h"
#include <stdexcept>
#include <cstdlib>
#include <atomic>
#include <optional>
#include "argparse.hpp"

std::mutex file_mutex;

class ThreadPool {
public:
    ThreadPool(size_t numThreads) : stop(false), tasksLeft(0), tasksInProgress(0) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queueMutex);
                        this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty()) {
                            return;
                        }
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                        ++tasksInProgress;
                    }

                    task();

                    {
                        std::unique_lock<std::mutex> lock(this->queueMutex);
                        --tasksInProgress;
                        --tasksLeft;
                        if (tasksInProgress == 0 && tasksLeft == 0) {
                            tasksFinished.notify_all();
                        }
                    }
                }
            });
        }
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<typename std::invoke_result<F, Args...>::type> {
        using return_type = typename std::invoke_result<F, Args...>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queueMutex);

            if (stop) {
                throw std::runtime_error("enqueue on stopped ThreadPool");
            }

            tasks.emplace([task]() { (*task)(); });
            ++tasksLeft;
        }
        condition.notify_one();
        return res;
    }

    void wait() {
        std::unique_lock<std::mutex> lock(queueMutex);
        tasksFinished.wait(lock, [this] { return tasksLeft == 0 && tasksInProgress == 0; });
    }

    ~ThreadPool() {
        wait(); // Ensure all tasks are completed before stopping the pool
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queueMutex;
    std::condition_variable condition;
    std::condition_variable tasksFinished;
    std::atomic<bool> stop;
    std::atomic<size_t> tasksLeft;
    std::atomic<size_t> tasksInProgress;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct FastaEntry {
    std::string id;
    std::string sequence;
};

struct FastaReader {
    struct promise_type {
        std::optional<FastaEntry> current_value;
        
        FastaReader get_return_object() {
            return FastaReader{ std::coroutine_handle<promise_type>::from_promise(*this) };
        }

        std::suspend_always initial_suspend() { return {}; }
        std::suspend_always final_suspend() noexcept { return {}; }
        std::suspend_always yield_value(FastaEntry value) {
            current_value = std::move(value);
            return {};
        }
        void return_void() {}
        void unhandled_exception() { std::terminate(); }
    };

    using handle_type = std::coroutine_handle<promise_type>;
    handle_type coro;

    FastaReader(handle_type h) : coro(h) {}
    FastaReader(FastaReader&) = delete;
    FastaReader(FastaReader&& other) : coro(other.coro) {
        other.coro = nullptr;
    }

    ~FastaReader() {
        if (coro) coro.destroy();
    }

    struct iterator {
        handle_type coro;

        iterator(handle_type h) : coro(h) {}
        iterator& operator++() {
            coro.resume();
            if (coro.done()) {
                coro = nullptr;
            }
            return *this;
        }

        FastaEntry operator*() {
            return *coro.promise().current_value;
        }

        bool operator==(std::default_sentinel_t) const { return !coro || coro.done(); }
        bool operator!=(std::default_sentinel_t) const { return coro && !coro.done(); }
    };

    iterator begin() {
        coro.resume();
        return iterator{ coro };
    }

    std::default_sentinel_t end() {
        return {};
    }
};

FastaReader read_fasta(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    FastaEntry entry;
    bool first = true;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!first) {
                co_yield entry;
            }
            entry.id = line.substr(1);
            entry.sequence.clear();
            first = false;
        } else {
            entry.sequence += line;
        }
    }

    if (!first) {
        co_yield entry;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
class SuffixArray {
public:
    using size_type = unsigned;
    using pointer = size_type*;
    using const_pointer = const size_type*;

private:
    // 比较 s 的两个长度相同的子串是否相同
    template<typename S>
    inline static bool substring_equal(const S& s, size_type p1, size_type p2, size_type len) {
        for (size_type i = 0;i < len;++i)
            if (s[p1 + i] != s[p2 + i])
                return false;
        return true;
    }

    // 诱导排序
    // 诱导排序需要使用 L 型桶和 S 型桶
    // 由于一轮操作中只会用到其中一种
    // 因此采用了互相推导的机制
    // 传入时只要求 lbuk 的内容正确
    template<typename S>
    inline static void induced_sort(
        const S& s,
        pointer sa,
        bool* type,
        pointer pos,
        pointer lbuk,
        pointer sbuk,
        size_type n,
        size_type m,
        size_type n0) {
        std::fill_n(sa, n, 0);
        // 用 lbuk 的信息填充 sbuk 的信息
        std::copy_n(lbuk + 1, m - 1, sbuk);
        sbuk[m - 1] = n;
        for (size_type i = n0;i-- > 0;)
            sa[--sbuk[s[pos[i]]]] = pos[i];
        // sbuk 的内容被修改，用 lbuk 的内容修正
        std::copy_n(lbuk + 1, m - 1, sbuk);
        sbuk[m - 1] = n;
        sa[lbuk[s[n - 1]]++] = n - 1;
        for (size_type i = 0;i < n;++i)
            if (sa[i] > 0 && !type[sa[i] - 1])
                sa[lbuk[s[sa[i] - 1]]++] = sa[i] - 1;
        // lbuk 的内容被修改，用 sbuk 的内容修正
        std::copy_n(sbuk, m - 1, lbuk + 1);
        lbuk[0] = 0;
        for (size_type i = n;i-- > 0;)
            if (sa[i] > 0 && type[sa[i] - 1])
                sa[--sbuk[s[sa[i] - 1]]] = sa[i] - 1;
        // 诱导排序结束，此时 lbuk 的内容是正确的
    }

public:
    // SA-IS 算法，构造后缀数组
    // 为了优化空间，用了一些重用空间的技巧
    template<typename S>
    static void sais(
        const S& s,
        pointer sa,
        bool* type,
        pointer len,
        pointer pos,
        pointer lbuk,
        pointer sbuk,
        size_type n,
        size_type m) {
        // 计算每个后缀的类型
        // 发现 LMS 字符则保存 LMS 字符的位置
        // 并计算以这个字符为起始的 LMS 子串的长度
        // LMS 字符位置保存在 pos 中，LMS 子串长度保存在 len 中
        // 保存长度的好处是在求出所有 LMS 子串的排序后的重命名阶段，比较时长度不相等可以直接判断不等
        pos += n / 2;
                type[n - 1] = false;
                std::fill_n(len, n, 0);
                const pointer npos = pos;
                for (size_type p = n - 1, i = n - 1;i-- > 0;) {
                        type[i] = s[i] != s[i + 1] ? s[i] < s[i + 1] : type[i + 1];
                        if (!type[i] && type[i + 1]) {
                                len[i + 1] = p - i;
                                *--pos = i + 1;
                                p = i + 1;
                        }
                }
                const size_type n0 = npos - pos;
        // 统计每个字符出现次数，计算 lbuk
        std::fill_n(lbuk, m, 0);
                for (size_type i = 0;i < n;++i)
                        ++lbuk[s[i]];
                for (size_type sum = 0, i = 0;i < m;++i) {
                        const size_type cur = lbuk[i];
                        lbuk[i] = sum;
                        sum += cur;
                }
        // 进行第一轮诱导排序，求出所有 LMS 子串的排序
        induced_sort(s, sa, type, pos, lbuk, sbuk, n, m, n0);
        // 根据排名对所有 LMS 子串重命名
        // 重命名后的排名记录在 len 里
        // 覆盖之前记录的长度信息
        // 因为长度信息已经用不到了
        // 可以重用记录长度的数组节约空间
        size_type m0 = -1;
        size_type ppos = -1, plen = 0;
        for (size_type i = 0;i < n;++i) {
            if (len[sa[i]] == 0) continue;
            if (len[sa[i]] != plen || !substring_equal(s, sa[i], ppos, plen)) ++m0;
            plen = len[sa[i]];
            len[sa[i]] = m0;
            ppos = sa[i];
        }
        // 重命名后的 s0、sa0 数组长度小于 n / 2（LMS 子串的性质）
        // 因此长度和小于 n
        // 可以直接塞到 sa 数组里
        // 同样是重用数组节约空间
        const pointer s0 = sa;
        const pointer sa0 = sa + n0;
        for (size_type i = 0;i < n0;++i)
            s0[i] = len[pos[i]];
        // 递归计算 LMS 后缀的排序
        // 当 LMS 子串各不相同时可以直接得出结果，不必递归
        // 注意此处 lbuk、sbuk 参数传的分别是 sbuk、sbuk + n0
        // 这样可以避免 lbuk 的内容被破坏
        // 由于 LMS 子串的性质
        // 并不需要 sbuk 特意留出空间
        // 另外注意 type、pos 这两个参数的传法
        // 也是为了避免当前的这两个数据被破坏
        // 避免重复计算
        // 但没有可以互相推导的目标
        // 只能在开空间时特意留出足够空间
        if (++m0 < n0)
            sais(s0, sa0, type + n, len, npos, sbuk, sbuk + n0, n0, m0);
        else for (size_type i = 0;i < n0;++i)
            sa0[s0[i]] = i;
        // LMS 后缀需要按顺序加入 sa 数组中
        // 因此需要根据 sa0 的结果调整顺序
        for (size_type i = 0;i < n0;++i)
            npos[i] = pos[sa0[i]];
        // 第二次诱导排序，得到最终的后缀数组
        induced_sort(s, sa, type, npos, lbuk, sbuk, n, m, n0);
    }

private:
    std::unique_ptr<size_type[]> data;

public:
    const_pointer sa, rk, ht;

public:
    // 包装 SA-IS 算法，在 SA-IS 算法计算 sa 的基础上，计算出 rk、ht
    template<typename S>
    SuffixArray(const S& s, size_type n, size_type m)
        : data(std::make_unique<size_type[]>(3 * n)) {
        const auto type = std::make_unique<bool[]>(2 * n);
        const auto lbuk = std::make_unique<size_type[]>(m);
        const auto sbuk = std::make_unique<size_type[]>(std::max(n, m));
        const pointer sa = data.get(), rk = sa + n, ht = rk + n;
        sais(s, sa, type.get(), rk, ht, lbuk.get(), sbuk.get(), n, m);
        for (size_type i = 0;i < n;++i)
            rk[sa[i]] = i;
        for (size_type k = 0, i = 0;i < n;++i) {
            if (rk[i] == 0) continue;
            if (k > 0) --k;
            const size_type j = sa[rk[i] - 1];
            const size_type l = n - std::max(i, j);
            for (;k < l;++k) if (s[i + k] != s[j + k]) break;
            ht[rk[i]] = k;
        }
        this->sa = sa;
        this->rk = rk;
        this->ht = ht;
    }

    inline size_type suffix(size_type p) const {
        return sa[p];
    }

    inline size_type rank(size_type p) const {
        return rk[p];
    }

    inline size_type height(size_type p) const {
        return ht[p];
    }
};

struct dupSeq {
    std::vector<int> inv_lens;
    std::vector<int> dup_lens;
    std::vector<std::vector<SuffixArray::size_type>> inv_poss;
    std::vector<std::vector<SuffixArray::size_type>> dup_poss;
    std::vector<std::string> inv_seqs;
    std::vector<std::string> dup_seqs;

};

class SortVector {
private:
    int threshold;

public:
    // 构造函数，初始化阈值
    SortVector(int t = INT_MIN) : threshold(t) {}

    // 设置阈值
    void setThreshold(int t) {
        threshold = t;
    }

    std::vector<int> sortArray(std::vector<int>& nums) {
        // 使用unordered_set进行去重
        std::unordered_set<int> unique_nums(nums.begin(), nums.end());
        std::vector<int> result(unique_nums.begin(), unique_nums.end());

        // 查找最大值和最小值
        int max = *std::max_element(result.begin(), result.end());
        int min = *std::min_element(result.begin(), result.end());

        // 创建计数数组
        std::vector<int> count(max - min + 1, 0);
        
        // 填充计数数组
        for(const int &el: result) {
            ++count[el - min];
        }

        // 累加计数数组
        int precount = 0;
        for(int i = 0; i < count.size(); ++i) {
            int temp = count[i];
            count[i] = precount;
            precount += temp;
        }

        // 根据计数数组构建结果数组
        std::vector<int> sorted_result(result.size());
        for(const int &el: result) {
            sorted_result[count[el - min]++] = el;
        }

        // 从大到小排序
        std::sort(sorted_result.begin(), sorted_result.end(), std::greater<int>());

        // 过滤掉小于阈值的元素
        std::vector<int> final_result;
        for (const int &el : sorted_result) {
            if (el > threshold) {
                final_result.push_back(el);
            }
        }

        return final_result;
    }
};

std::pair<SuffixArray::size_type, SuffixArray::size_type> min_max(SuffixArray::size_type a, SuffixArray::size_type b){
    if (a >b){
        return {b, a};
    }else{
        return {a, b};
    };
}

class Solution {
public:
    //std::string longestDupSubstring(std::string s) {
    dupSeq longestDupSubstring(int nu_i, const std::string& s, int win, int min_repeat_size, int min_dist) {
        const int n = s.size();
        SuffixArray sa(s, n, 128);
        //std::vector<int> poss;
	std::vector<std::vector<SuffixArray::size_type>> inv_poss;
	std::vector<std::vector<SuffixArray::size_type>> dup_poss;
	//std::vector<int> inv_lens;
	//std::vector<int> dup_lens;

	dupSeq dp;
        for (int i = 1;i < n;++i) {
	    int tmp = std::abs(static_cast<int>(sa.sa[i]) - static_cast<int>(sa.sa[i-1]));
            if (sa.ht[i] >= min_repeat_size) {
                auto [min_val, max_val]  = min_max(sa.sa[i], sa.sa[i-1]);
                if (max_val >= win){
                    if (min_val < win) {
			//std::cout << "Solution: " << n << " " << chrom << " " << sa.sa[i] << " " << sa.sa[i-1] << " | " << sa.sa[i+1] << " " << sa.ht[i] << " " << nu_i << " " << min_val << " " << max_val << std::endl;
			SuffixArray::size_type a = nu_i+min_val;
			SuffixArray::size_type b = nu_i+n-max_val-sa.ht[i];
			if (a <= b){
		            inv_poss.push_back({a, b, 1, 2, sa.ht[i]});
			}else{
			    inv_poss.push_back({b, a, 2, 1, sa.ht[i]});
			}
			//auto [a, b] = min_max(min_val, 2*win-max_val-1-sa.ht[i]);
			//SuffixArray::size_type a = min_val;
			//SuffixArray::size_type b = 2*win-max_val-1-sa.ht[i];
		        //inv_poss.push_back({a, b});
			//inv_lens.push_back(sa.ht[i]);
		    };
		}else if (tmp >= min_dist){
		    dup_poss.push_back({nu_i+min_val, nu_i+max_val, sa.ht[i]});
                    //dup_lens.push_back(sa.ht[i]);	
		};
	    }else {
	        continue;
	    };
	};
        //dp.inv_lens = inv_lens;
	//dp.dup_lens = dup_lens;	
       	// 使用 std::sort 和 lambda 表达式排序
        std::sort(inv_poss.begin(), inv_poss.end(),
              [](const std::vector<SuffixArray::size_type>& a, const std::vector<SuffixArray::size_type>& b) {
                  return a[0] < b[0];
              });
        dp.inv_poss = std::move(inv_poss);
        std::sort(dup_poss.begin(), dup_poss.end(),
              [](const std::vector<SuffixArray::size_type>& a, const std::vector<SuffixArray::size_type>& b) {
                  return a[0] < b[0];
              });
        dp.dup_poss = std::move(dup_poss);
	return dp;
    };
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> sortArrayWithIndex(const std::vector<int>& nums) {
        int n = nums.size();
	std::cout << "g1_1" << std::endl;
	std::cout << n << std::endl;
        std::vector<std::pair<int, int>> indexed_nums(n); // 用于存储元素及其原始下标
        
        // 绑定元素及其原始下标
        for (int i = 0; i < n; ++i) {
            indexed_nums[i] = {nums[i], i};
        }
        std::cout << "g1_2" << std::endl;
        // 对 indexed_nums 按元素值进行排序
        std::sort(indexed_nums.begin(), indexed_nums.end());

        // 提取排序后每个元素的原始下标
        std::vector<int> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = indexed_nums[i].second;
        }

        return result;
    };


void replace_invalid_chars(std::string& sequence) {
    const std::string valid_chars = "ATCGN";
    for (char& ch : sequence) {
        if (valid_chars.find(ch) == std::string::npos) {
            ch = 'N';
        }
    }
}

std::string reverse_and_concatenate(std::string& input) {
    replace_invalid_chars(input);
    std::string reversed = input;
    std::ranges::reverse(reversed);
    std::string result = input + "X" + reversed;
    return result;
}

// 函数：获取DNA碱基的互补碱基
char getComplement(char base) {
    switch(base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N'; // 处理非标准碱基
    }
}

// 函数：计算给定DNA序列的反向互补序列
std::string getReverseComplement(const std::string& dna) {
    std::string reverseComplement;
    // 反向遍历原始DNA序列，同时获取互补碱基
    for (auto it = dna.rbegin(); it != dna.rend(); ++it) {
        reverseComplement += getComplement(*it);
    }
    return dna+'X'+reverseComplement;
}


void write_bgzip_file(const std::string filename, const char* data) {
    const char* tmp_file = filename.c_str();
    BGZF* fp = bgzf_open(tmp_file, "a");
    if (fp == nullptr) {
        std::cerr << "Failed to open BGZF file for writing" << std::endl;
        return;
    }

    int length = strlen(data);
    if (bgzf_write(fp, data, length) < 0) {
        std::cerr << "Failed to write data to BGZF file" << std::endl;
    }

    if (bgzf_close(fp) < 0) {
        std::cerr << "Failed to close BGZF file" << std::endl;
    }
}

void create_tabix_index_oldx(const std::string filename) {
     // 定义要执行的命令
    std::string tmp = "tabix "+filename; // 对于Windows系统，可以使用 "dir"
    const char* command = tmp.c_str();
    // 使用 system() 函数执行命令
    int result = system(command);
    
    // 检查命令执行的结果
    if (result == 0) {
        std::cout << "Command executed successfully." << std::endl;
    } else {
        std::cerr << "Command execution failed with error code: " << result << std::endl;
    }

}

void create_tabix_index(const std::string filename) {
    const char* tmp_file = filename.c_str();
    std::cout << "Entering create_tabix_index function" << std::endl;

    // Configuration for BED files (UCSC)
    tbx_conf_t conf = tbx_conf_bed;

    std::cout << "Configuration set, attempting to build index for file: " << filename << std::endl;

    int result = tbx_index_build(tmp_file, 0, &conf);
    if (result < 0) {
        std::cerr << "Failed to create Tabix index" << std::endl;
    } else {
        std::cout << "Successfully created Tabix index" << std::endl;
    }

    std::cout << "Exiting create_tabix_index function" << std::endl;
}


std::string conect_to_string(const std::string& contig, SuffixArray::size_type pos1, SuffixArray::size_type pos2, SuffixArray::size_type len){
    std::stringstream ss;
    //ss << contig << "\t" << i << "\t" << i+win << "\t" << pos1 << "," << pos2 << "\t" << len;
    ss << contig << "\t" << pos1 << "\t" << pos2+len << "\t" << len;
    std::string result = ss.str();
    return result;
}

std::string conect_to_inv_string(const std::string& contig, SuffixArray::size_type pos1, SuffixArray::size_type pos2, SuffixArray::size_type r1, SuffixArray::size_type r2, SuffixArray::size_type len){
    std::stringstream ss;
    //ss << contig << "\t" << i << "\t" << i+win << "\t" << pos1 << "," << pos2 << "\t" << len;
    ss << contig << "\t" << pos1 << "\t" << pos2+len << "\t" << r1 << "\t" << r2 <<"\t" << len << "\t" << (pos1+pos2+len)/2;
    std::string result = ss.str();
    return result;
}

struct winData {
    std::vector<std::vector<SuffixArray::size_type>> lens;
    std::vector<std::vector<std::vector<SuffixArray::size_type>>> positions;
    std::vector<std::vector<int>> usetosortpos;
    std::string contig;
    int seq_size;
    int win;
    std::string filename1;
    std::string filename2;
};

void dealDUP(const winData& win_dat ){
    std::ostringstream dat1;
    std::ostringstream dat2;
    for (size_t i=0; i < win_dat.seq_size; ++i){
        if (win_dat.positions[i].empty()){
            continue;
        };
        //std::vector<int> sortedindex = sortArrayWithIndex(win_dat.usetosortpos[i]);
        int head_j = 0;
        std::string tmpline = conect_to_string(win_dat.contig, win_dat.positions[i][head_j][0], win_dat.positions[i][head_j][1], win_dat.positions[i][head_j][2]);
        dat1 << tmpline << "\n";
        bool flag = false;
	if (win_dat.positions.size() >1){
            for (size_t ii=1; ii<win_dat.positions[i].size(); ++ii){
                int j = ii;
                std::vector<SuffixArray::size_type> head = win_dat.positions[i][head_j];
                std::vector<SuffixArray::size_type> next = win_dat.positions[i][j];
                std::string tmp = conect_to_string(win_dat.contig, win_dat.positions[i][j][0], win_dat.positions[i][j][1], win_dat.positions[i][j][2]);
                dat1 << tmp << "\n";
                if ((head[0] <= next[1]) && (head[1] >= next[0])){
                    int head_abs = std::abs(static_cast<int>(head[0]) - static_cast<int>(head[1]));
                    int next_abs = std::abs(static_cast<int>(next[0]) - static_cast<int>(next[1]));
                    if ((win_dat.positions[i][j][2] > win_dat.positions[i][head_j][2]) || ((win_dat.positions[i][j][2] == win_dat.positions[i][head_j][2]) && (next_abs > head_abs))) {
                        tmpline = tmp;
                        head_j = j;
                    };
                    flag = false;
                }else{
                    dat2 << tmpline << "\n";
                    head_j = j;
                    tmpline = conect_to_string(win_dat.contig, win_dat.positions[i][j][0], win_dat.positions[i][j][1], win_dat.positions[i][j][2]);
                    flag = true;
                };
            };
            if (flag){
                int index_j = win_dat.positions[i].size()-1;
                tmpline = conect_to_string(win_dat.contig, win_dat.positions[i][index_j][0], win_dat.positions[i][index_j][1], win_dat.positions[i][index_j][2]);
            };
            dat2 << tmpline << "\n";
	}else{
	    dat2 << tmpline << "\n";
	};
    };

    std::string dat1_str = dat1.str();
    std:: string dat2_str = dat2.str();
    {
        std::lock_guard<std::mutex> guard(file_mutex); // 自动加锁和解锁
        write_bgzip_file(win_dat.filename1, dat1_str.c_str());
        write_bgzip_file(win_dat.filename2, dat2_str.c_str());
    };
};

void dealINV(const winData& win_dat ){
    std::ostringstream dat1;
    std::ostringstream dat2;
    for (size_t i=0; i < win_dat.seq_size; ++i){
        if (win_dat.positions[i].empty()){
            continue;
        };
	//std::cout << "g1" << "\t" << win_dat.positions.size() << "\t" << i << "\t"  << win_dat.positions[i].empty() << "\t" << win_dat.positions[i].size() <<  std::endl;
        //std::vector<int> sortedindex = sortArrayWithIndex(win_dat.usetosortpos[i]);
        int head_j = 0;
	//std::cout << "g11" << std::endl;
        std::string tmpline = conect_to_inv_string(win_dat.contig, win_dat.positions[i][head_j][0], win_dat.positions[i][head_j][1], win_dat.positions[i][head_j][2], win_dat.positions[i][head_j][3], win_dat.positions[i][head_j][4]);
        dat1 << tmpline << "\n";
	std::cout << "g2" << std::endl;
        bool flag = false;
	if (win_dat.positions.size() >1){ 
            for (size_t ii=1; ii<win_dat.positions[i].size(); ++ii){
                int j = ii;
		//std::cout << "gg " << head_j << " " << j << " " << win_dat.positions[i].size() << std::endl;
                std::vector<SuffixArray::size_type> head = win_dat.positions[i][head_j];
                std::vector<SuffixArray::size_type> next = win_dat.positions[i][j];
		//std::cout << win_dat.positions[i][j][4] << " " << win_dat.positions[i][head_j][4] <<std::endl;
                std::string tmp = conect_to_inv_string(win_dat.contig, win_dat.positions[i][j][0], win_dat.positions[i][j][1], win_dat.positions[i][j][2], win_dat.positions[i][j][3], win_dat.positions[i][j][4]);
                dat1 << tmp << "\n";
		//std::cout << "g3" << std::endl;
                if (head[0] == head[1] && next[0] == next[1]) {
                    if ((head[1] - head[0]) == (next[1] - next[0])){
                        int head_abs = std::abs(static_cast<int>(head[0]) - static_cast<int>(head[1]));
                        int next_abs = std::abs(static_cast<int>(next[0]) - static_cast<int>(next[1]));
			//std::cout << "g4" << std::endl;
                        if (win_dat.positions[i][j][4] > win_dat.positions[i][head_j][4]) {
                            tmpline = tmp;
                            head_j = j;
                        };
                        flag = false;
                    };
                    flag = false;
                }else if ((next[0] - head[0]) <= 50) {
		    int head_abs = std::abs(static_cast<int>(head[0]) - static_cast<int>(head[1]));
                    int next_abs = std::abs(static_cast<int>(next[0]) - static_cast<int>(next[1]));
                    //g5 0 1413 1412 5 5
		    //std::cout << "g5 " << i << " " << j << " " << head_j << " " << win_dat.positions[i][j].size() << " " << win_dat.positions[i][head_j].size()<<std::endl;
		    //std::cout << win_dat.positions[i][j][4] << " " << win_dat.positions[i][head_j][4] <<std::endl;
                    if (win_dat.positions[i][j][4] > win_dat.positions[i][head_j][4]) {
                        tmpline = tmp;
                        head_j = j;
                    };
                    flag = false;
                }else{
		    //std::cout << "g6" << std::endl;
                    dat2 << tmpline << "\n";
                    head_j = j;
                    tmpline = conect_to_inv_string(win_dat.contig, win_dat.positions[i][j][0], win_dat.positions[i][j][1], win_dat.positions[i][j][2], win_dat.positions[i][j][3], win_dat.positions[i][j][4]);
                    flag = true;
                };
            };
            if (flag){
		//std::cout << "g7" << std::endl;
                int index_j = win_dat.positions[i].size()-1;
                tmpline = conect_to_inv_string(win_dat.contig, win_dat.positions[i][index_j][0], win_dat.positions[i][index_j][1], win_dat.positions[i][index_j][2], win_dat.positions[i][index_j][3], win_dat.positions[i][index_j][4]);
            };
            dat2 << tmpline << "\n";
	}else{
	    dat2 << tmpline << "\n";
	};
    };
    std::cout << "g3" << std::endl;
    std::string dat1_str = dat1.str();
    std:: string dat2_str = dat2.str();
    {
        std::lock_guard<std::mutex> guard(file_mutex); // 自动加锁和解锁
        write_bgzip_file(win_dat.filename1, dat1_str.c_str());
        write_bgzip_file(win_dat.filename2, dat2_str.c_str());
    };
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void processWindow_genome(const std::string& sequence, const std::string& contig, int win, int min_repeat_size,  int step, int min_dist, const std::string inv_file1, const std::string inv_file2, const std::string dup_file1, const std::string dup_file2) {
    Solution solution;
    std::vector<std::vector<int>> inv_usetosortpos;
    std::vector<std::vector<int>> dup_usetosortpos;
    std::vector<std::vector<SuffixArray::size_type>> inv_lens;
    std::vector<std::vector<SuffixArray::size_type>> dup_lens;
    std::vector<std::vector<std::vector<SuffixArray::size_type>>> inv_positions;
    std::vector<std::vector<std::vector<SuffixArray::size_type>>> dup_positions;
    
    int seq_size = sequence.length();
    winData win_dat_inv;
    winData win_dat_dup;
    win_dat_inv.win = win;
    win_dat_dup.win = win;
    win_dat_inv.contig = contig;
    win_dat_dup.contig = contig;
    win_dat_inv.filename1 = inv_file1;
    win_dat_inv.filename2 = inv_file2;
    win_dat_dup.filename1 = dup_file1;
    win_dat_dup.filename2 = dup_file2;

    int nums = 0;
    for (size_t i = 0; i <= seq_size; i += step) {
	nums += 1;
        std::cout << win << std::endl;
        int end = std::min(static_cast<int>(i) + win, seq_size);
        std::string seq = std::string(sequence.substr(i, end));
        std::cout << "1" << std::endl;
        for (char& c : seq) {
            c = std::toupper(c);
        };
        std::cout << "2" << std::endl;
        //std::string seq_inv = reverse_and_concatenate(seq);
	std::string seq_inv = getReverseComplement(seq);
        std::cout << "a" << std::endl;

        dupSeq dp = solution.longestDupSubstring(i, seq_inv, win+1, min_repeat_size, min_dist);

        std::vector<std::string> line;
        std::vector<int> inv_tmppos;
	std::vector<int> dup_tmppos;
        std::vector<std::vector<SuffixArray::size_type>> inv_pos;
	std::vector<std::vector<SuffixArray::size_type>> dup_pos;
        std::vector<SuffixArray::size_type> inv_len;
	std::vector<SuffixArray::size_type> dup_len;
        std::cout << "b" << std::endl;

        for (size_t ii = 0; ii < dp.inv_poss.size(); ++ii) {
            //std::stringstream ss;
            //ss << contig << "\t" << i << "\t" << i+win << "\t" << dp.poss[ii][0] << "," << dp.poss[ii][1] << "\t" << dp.lens[ii];
            //std::string result = ss.str();
            //std::cout << result << "\n" << std::endl;
            //line.push_back(result);
            //inv_tmppos.push_back(static_cast<int>(dp.inv_poss[ii][0]+(dp.inv_poss[ii][0] + dp.inv_poss[ii][1])/2));
            inv_pos.push_back(dp.inv_poss[ii]);
            //inv_len.push_back(dp.inv_lens[ii]);
        };

	for (size_t ii = 0; ii < dp.dup_poss.size(); ++ii) {
            //dup_tmppos.push_back(static_cast<int>(dp.dup_poss[ii][0]+(dp.dup_poss[ii][0] + dp.dup_poss[ii][1])/2));
            dup_pos.push_back(dp.dup_poss[ii]);
            //dup_len.push_back(dp.dup_lens[ii]);
        };
        std::cout << "c" << std::endl;
        //inv_usetosortpos.push_back(inv_tmppos);
	//dup_usetosortpos.push_back(dup_tmppos);
        inv_positions.push_back(inv_pos);
	dup_positions.push_back(dup_pos);
        //inv_lens.push_back(inv_len);
	//dup_lens.push_back(dup_len);
	std::cout << "e" << std::endl;
    };

    win_dat_inv.seq_size = nums;
    win_dat_dup.seq_size = nums;

    std::cout << "f" << std::endl; 
    //win_dat_inv.usetosortpos = std::move(inv_usetosortpos);
    win_dat_inv.positions = std::move(inv_positions);
    //win_dat_inv.lens = std::move(inv_lens);
    //win_dat_dup.usetosortpos = std::move(dup_usetosortpos);
    win_dat_dup.positions = std::move(dup_positions);
    //win_dat_dup.lens = std::move(dup_lens);
    std::cout << "g" << std::endl;
    dealINV(win_dat_inv);
    std::cout << "h" << std::endl;
    dealDUP(win_dat_dup);
    std::cout << "i" << std::endl;
};

void is_overwrite(std::vector<std::string> files, bool overwrite){
    if (overwrite){
        for (size_t i = 0;i<files.size();i++){
	    std::filesystem::path file_path{files[i]};
	    if (std::filesystem::exists(file_path)) {
                try {
                    if (std::filesystem::remove(file_path)) {
                        std::cout << "File deleted successfully.\n";
                    }
                } catch (const std::filesystem::filesystem_error& e) {
                    std::cerr << "Error: " << e.what() << '\n';
                }
            }

	}
    }
}

void printNestedVectors(const std::vector<std::vector<std::vector<SuffixArray::size_type>>>& nestedVectors) {
    for (const auto& outerVector : nestedVectors) {
        for (const auto& middleVector : outerVector) {
            for (const auto& innerValue : middleVector) {
                std::cout << innerValue << "\t";
            }
            std::cout << "\n";  // 分隔每个中间向量
        }
        //std::cout <<  std::endl;  // 每个外层向量输出为一行
    }
}

void processWindow_seq(const std::string& sequence, int win, int min_repeat_size, int step, int min_dist, bool inv_or_dup) {
    Solution solution;
    std::vector<std::vector<int>> inv_usetosortpos;
    std::vector<std::vector<int>> dup_usetosortpos;
    std::vector<std::vector<SuffixArray::size_type>> inv_lens;
    std::vector<std::vector<SuffixArray::size_type>> dup_lens;
    std::vector<std::vector<std::vector<SuffixArray::size_type>>> inv_positions;
    std::vector<std::vector<std::vector<SuffixArray::size_type>>> dup_positions;

    int seq_size = sequence.length();

    int nums = 0;
    for (size_t i = 0; i <= seq_size; i += step) {
        int end = std::min(static_cast<int>(i) + win, seq_size);
   	std::string seq = std::string(sequence.substr(i,end));
        for (char& c : seq) {
            c = std::toupper(c);
        };
        std::string seq_inv = getReverseComplement(seq);
	dupSeq dp = solution.longestDupSubstring(i, seq_inv, win+1, min_repeat_size, min_dist);

        std::vector<std::vector<SuffixArray::size_type>> inv_pos;
        std::vector<std::vector<SuffixArray::size_type>> dup_pos;
        for (size_t ii = 0; ii < dp.inv_poss.size(); ++ii) {
            inv_pos.push_back(dp.inv_poss[ii]);
        };

        for (size_t ii = 0; ii < dp.dup_poss.size(); ++ii) {
            dup_pos.push_back(dp.dup_poss[ii]);
        };
        inv_positions.push_back(inv_pos);
        dup_positions.push_back(dup_pos);
    };

    if (inv_or_dup){
        printNestedVectors(inv_positions);
    }else{
        printNestedVectors(dup_positions);
    };

};

std::string readSingleLineStringFromFile(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;

    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return "";
    }

    // 读取文件中的第一行
    if (std::getline(file, line)) {
        file.close();
        return line;
    } else {
        std::cerr << "文件为空或读取失败: " << filename << std::endl;
        file.close();
        return "";
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("SeqRepeater");

    program.add_argument("-p", "--prefix")
          //.required()
          .help("The output file prefix.");

    program.add_argument("-g", "--genome")
          //.required()
          .help("The fasta file of genome.");

    program.add_argument("-q", "--seq")
          //.required()
          .help("The seqence.");

    program.add_argument("-r","--inv_or_dup")
          .default_value(false)
          .implicit_value(true)
	  .help("The seqence.");

    program.add_argument("-d","--min_dist")
	  .scan<'d', int>()
	  .nargs(1)
	  .required()
          .default_value(5)
          .help("The min distance between two repeat seq.");

    program.add_argument("-w", "--win")
          .scan<'d', int>()
          .default_value(100000)
          .help("The window size.");

    program.add_argument("-s", "--step")
          .scan<'d', int>()
          .default_value(50000)
          .help("The window move step size.");

    program.add_argument("-t", "--thread")
          .scan<'d', int>()
          .default_value(20)
          .help("Number of threads.");

    program.add_argument("-m", "--min_size")
          .scan<'d', int>()
          .default_value(1)
          .help("Min repeat seq size.");

    program.add_argument("-f", "--force")
      .default_value(true)
      .implicit_value(false)
      .help("Force overwrite the files without asking. default: true");

     // 手动添加帮助信息
    /*program.add_description("This program does something.");
    program.add_argument("-h", "--help")
           .action([&program](const std::string& value) {
               std::cout << "Usage: program_name [options]" << std::endl;
               std::cout << "Options:" << std::endl;
               std::cout << "  -d, --min_dist INT  The min distance between two repeat seq. (required)" << std::endl;
               std::cout << "  -h, --help          Show this help message and exit" << std::endl;
               std::exit(0);
           })
           .help("Show this help message and exit");
    */

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    //std::cout  << "ddddddd" << std::endl;
     // Check if both --seq and --genome are provided
    if (program.present("--seq") && program.present("--genome")) {
        std::cerr << "Error: --seq and --genome cannot be used together." << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    // Check if neither --seq nor --genome are provided
    if (!program.present("--seq") && !program.present("--genome")) {
        std::cerr << "Error: Either --seq or --genome must be provided." << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    // 使用长名称来获取参数值
    //std::string filename = program.get<std::string>("--genome");
    //std::string prefix = program.get<std::string>("--prefix");
    //std::string seq = program.get<std::string>("--seq");
    //int win = program.get<int>("--win");
    //int thread_size = program.get<int>("--thread");
    //int min_repeat_size = program.get<int>("--minSize");
    //bool overwrite = program.get<bool>("--force");
    //std::cout  << "ddddddd" << std::endl;
    std::optional<std::string> genome_tmp = program.present("--genome") ? std::optional<std::string>(program.get<std::string>("--genome")) : std::nullopt;
    std::optional<std::string> prefix_tmp = program.present("--prefix") ? std::optional<std::string>(program.get<std::string>("--prefix")) : std::nullopt;
    std::optional<std::string> seq_tmp = program.present("--seq") ? std::optional<std::string>(program.get<std::string>("--seq")) : std::nullopt;
    int win = program.get<int>("--win");
    int step = program.get<int>("--step");
    int thread_size = program.get<int>("--thread");
    int min_repeat_size = program.get<int>("--min_size");
    int min_dist = program.get<int>("--min_dist");
    bool overwrite = program.get<bool>("--force");
    bool inv_or_dup = program.get<bool>("--inv_or_dup");
    
    //std::cout <<  win  << " " << thread_size << " " <<  std::endl;
    
    if (seq_tmp.has_value()) {
        //std::string seq = seq_tmp.value_or("No seq");
	std::string seq = readSingleLineStringFromFile(seq_tmp.value_or("No seq"));
        processWindow_seq(seq, win, min_repeat_size, step, min_dist, inv_or_dup);
    }else{
	std::string genome = genome_tmp.value_or("No genome file");
        std::string prefix = prefix_tmp.value_or("No prefix"); 
        ThreadPool pool(thread_size);

        const std::string inv_file1 = prefix + ".inv1.bed.gz";
        const std::string inv_file2 = prefix + ".inv2.bed.gz";
        const std::string dup_file1 = prefix + ".dup1.bed.gz";
        const std::string dup_file2 = prefix + ".dup2.bed.gz";
        is_overwrite(std::vector<std::string> {inv_file1, inv_file2, dup_file1, dup_file2}, overwrite);

        for (const auto& entry : read_fasta(genome)) {
    	    std::cout << entry.id << std::endl;
            pool.enqueue(processWindow_genome, entry.sequence, entry.id, win, min_repeat_size, min_dist, step, inv_file1, inv_file2, dup_file1, dup_file2);
        };

        // 等待所有线程池中的任务完成
        pool.wait();

        create_tabix_index(inv_file1);
        create_tabix_index(inv_file2);
        create_tabix_index(dup_file1);
        create_tabix_index(dup_file2);
    };

    return 0;
};

