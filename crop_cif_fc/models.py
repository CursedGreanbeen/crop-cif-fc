from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal


ChainType = Literal["H", "K", "L", "AG", "UNKNOWN"]
PlanAction = Literal["keep_all", "crop", "manual_review"]


@dataclass(frozen=True)
class ChainAnnotation:
    chain_type: ChainType
    range_start: int | None = None
    range_end: int | None = None

@dataclass(frozen=True)
class ChainCropPlan:
    action: PlanAction
    status: str
    chain_type: ChainType
    fasta_start: int | None = None
    fasta_end: int | None = None
    keep_cif_indices: tuple[int, ...] = ()
    missing_positions: tuple[int, ...] = ()
    reason: str | None = None


@dataclass(frozen=True)
class ResidueRef:
    chain_id: str
    cif_index: int
    name: str
    auth_seq_id: int
    insertion_code: str
    one_letter: str
    def __str__(self) -> str:
        return f"{self.auth_seq_id}{self.insertion_code.strip()}"
