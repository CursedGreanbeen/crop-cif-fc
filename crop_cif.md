# crop_cif_fc Fixes

## Imports & Module Structure

**crop_planner.py**
- Changed `from models import` → `from .models import` (relative import)
- Removed unused `ResidueRef` import

**pipeline.py**
- Added `build_residue_refs` to imports from `.cif_utils`
- Added `ProcessingResult` import (if needed)

## models.py

**ChainAnnotation** — removed field:
- `range_source: str | None = None` — deleted from dataclass

**ChainCropPlan** — removed field:
- `range_source` — deleted from all constructor calls

## pipeline.py

**process_structure()**
- Removed `range_source=row.get("range_source")` from ChainAnnotation construction
- Changed residue handling: now calls `build_residue_refs(chain_id, residues)` before `map_fasta_to_cif()`

## cif_utils.py

**load_structure()**
- Changed signature: `cif_path: str` → `cif_path: str | Path`
- Added `str()` conversion: `gemmi.read_structure(str(cif_path))`

**build_residue_refs()**
- Fixed field name: `residue_name=` → `name=` (matches ResidueRef dataclass)

**apply_structure_plans()**
- Fixed residue removal: `chain.remove_residue(residue.seqid)` → `del chain[cif_index]`
- Delete in reverse order to avoid index shifting

## anarci_annotator.py

**Added missing functions:**
- `is_antibody_chain_type(chain_type)` — checks if type is H, K, or L
- `chain_role(chain_type)` — returns "heavy", "light", or "other"

## run_crop_cif_fc.py

**Path fix:**
- `FASTA_DIR` — corrected to `testing/fasta-cropped-by-cif-test`

## Unused Functions (can be removed)

- `is_antibody_chain_type()` — defined but never called
- `annotate_fasta()` — not used (annotate_record is used instead)
- `get_chain_sequence()` — not called
- `get_residue_identifier()` — not called
- `build_residue_refs()` — now used after fix
- `remove_chain()` — not called
- `update_sequence()` — not called
- `write_fasta_file()` — not called
