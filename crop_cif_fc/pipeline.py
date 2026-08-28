from pathlib import Path

from .cif_utils import (
    load_structure,
    get_chain_residues,
    get_protein_chain_ids,
    apply_structure_plans,
    save_structure,
    build_residue_refs
)

from .fasta_manager import FastaRecord, read_fasta_file
from .annotation import annotate_record, expand_annotations_by_chain
from .models import ChainAnnotation
from .sequence_mapping import map_fasta_to_cif
from .crop_planner import plan_structure_crop
from .report import build_structure_report, write_report
from .chain_splitter import split_structure_chains


def process_structure(
    cif_path: str | Path,
    fasta_path: str | Path,
    output_path: str | Path,
    report_path: str | Path,
) -> None:
    cif_path = Path(cif_path)
    fasta_path = Path(fasta_path)
    output_path = Path(output_path)
    pdb_id = cif_path.stem

    print(f"[{pdb_id}] loading structure and FASTA...")

    structure = load_structure(cif_path)
    fasta_records: dict[str, FastaRecord] = read_fasta_file(fasta_path)

    split_report, new_chain_annotations, new_chain_mappings = split_structure_chains(
        structure=structure,
        fasta_records=fasta_records,
    )
    consumed_chain_ids = set(split_report.keys())

    annotations = []
    unique_records = list({id(r): r for r in fasta_records.values()}.values())

    for record in unique_records:
        if set(record.chain_ids) & consumed_chain_ids:
            continue
        annotations.append(annotate_record(record))

    expanded_annotations = expand_annotations_by_chain(annotations)
    chain_annotations = {
        row["chain_id"]: ChainAnnotation(
            chain_type=row["chain_type"],
            range_start=row.get("start"),
            range_end=row.get("end"),
        )
        for row in expanded_annotations
    }

    # Merge in annotations for the newly split chains.
    chain_annotations.update(new_chain_annotations)

    chain_mappings = {}

    for chain_id in get_protein_chain_ids(structure):
        chain_mappings[chain_id] = []

    for chain_id, annotation in chain_annotations.items():
        if chain_id in new_chain_mappings:
            continue  # already built by split_structure_chains, trivial 1:1

        if annotation.chain_type in {"H", "K", "L"}:
            residues = get_chain_residues(structure, chain_id)
            residue_refs = build_residue_refs(chain_id, residues)
            fasta_record = fasta_records[chain_id]

            mapping, alignment = map_fasta_to_cif(
                reference_sequence=fasta_record.sequence,
                cif_residues=residue_refs,
            )
            chain_mappings[chain_id] = mapping

    chain_mappings.update(new_chain_mappings)

    plans = plan_structure_crop(
        chain_annotations=chain_annotations,
        chain_mappings=chain_mappings,
    )

    for chain_id, plan in plans.items():
        print(f"[{pdb_id}] chain {chain_id}: {plan.action} ({plan.status})")

    apply_structure_plans(structure=structure, plans=plans)
    save_structure(structure=structure, path=output_path)

    write_report(
        rows=build_structure_report(
            pdb_id=pdb_id,
            annotations=chain_annotations,
            plans=plans,
            mappings=chain_mappings,
        ),
        path=report_path,
    )

    print(f"[{pdb_id}] done → {output_path.name}, {report_path.name}")
   