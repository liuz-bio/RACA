from dataclasses import dataclass

@dataclass(frozen=True)
class ConFig:
    genome: int
    winSize: int
    bamFiles: list
    sampleIDs: list
    sampleToBam: dict
