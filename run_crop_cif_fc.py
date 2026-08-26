from pathlib import Path

from crop_cif_fc.pipeline import process_structure


CIF_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/testing/CIFs-filtered-test")
FASTA_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/testing/fasta-cropped-by-cif-test")
OUTPUT_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/testing/CIFs-fv-only")
REPORT_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/results/reports")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for cif_path in sorted(CIF_DIR.glob("*.cif")):
    pdb_id = cif_path.stem
    fasta_path = FASTA_DIR / f"{pdb_id}.fasta"

    if not fasta_path.exists():
        print(f"[{pdb_id}] SKIP: no matching FASTA file")
        continue

    output_cif_path = OUTPUT_DIR / f"{pdb_id}_fv.cif"
    report_path = REPORT_DIR / f"{pdb_id}_report.tsv"

    try:
        process_structure(
            cif_path=cif_path,
            fasta_path=fasta_path,
            output_path=output_cif_path,
            report_path=report_path,
        )
    except Exception as error:
        print(f"[{pdb_id}] FAILED: {error}")
