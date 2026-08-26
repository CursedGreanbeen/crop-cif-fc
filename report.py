import csv
from pathlib import Path
from typing import Any

from .models import (
    ChainAnnotation,
    ChainCropPlan
)


def build_chain_report(
    pdb_id: str,
    chain_id: str,
    annotation: ChainAnnotation | None,
    plan: ChainCropPlan,
    mapping: list[dict],
) -> dict:
    mapped = [
        item for item in mapping
        if item["in_cif"] and item["cif_residue"] is not None
    ]

    kept_indices = set(plan.keep_cif_indices)

    kept_residues = [
        item["cif_residue"]
        for item in mapped
        if plan.action == "keep_all"
        or item["cif_index"] in kept_indices
    ]

    fasta_range = (
        f"{plan.fasta_start}-{plan.fasta_end}"
        if plan.fasta_start is not None and plan.fasta_end is not None
        else ""
    )

    cif_original_range = (
        f"{mapped[0]['cif_residue']}-{mapped[-1]['cif_residue']}"
        if mapped else ""
    )

    cif_kept_range = (
        f"{kept_residues[0]}-{kept_residues[-1]}"
        if kept_residues else ""
    )

    return {
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        "chain_type": (
            annotation.chain_type
            if annotation is not None
            else plan.chain_type
        ),
        "plan": plan.action,
        "reason": plan.reason,
        "fasta_range": fasta_range,
        "cif_original_range": cif_original_range,
        "cif_kept_range": cif_kept_range,
        "status": plan.status,
    }


def build_structure_report(
    pdb_id: str,
    annotations: dict[str, ChainAnnotation],
    plans: dict[str, ChainCropPlan],
    mappings: dict[str, list[dict]],
) -> list[dict]:
    chain_ids = sorted(set(annotations) | set(plans) | set(mappings))

    return [
        build_chain_report(
            pdb_id=pdb_id,
            chain_id=chain_id,
            annotation=annotations.get(chain_id),
            plan=plans[chain_id],
            mapping=mappings.get(chain_id, []),
        )
        for chain_id in chain_ids
        if chain_id in plans
    ]


REPORT_FIELDS = [
    "pdb_id",
    "chain_id",
    "chain_type",
    "plan",
    "reason",
    "fasta_range",
    "cif_original_range",
    "cif_kept_range",
    "status",
]


def write_report(
    rows: list[dict[str, Any]],
    path: str | Path,
    include_antigen: bool = False,
) -> None:
    """Write chain-processing results to a tab-separated report."""

    output_path = Path(path)

    filtered_rows = [
        row
        for row in rows
        if include_antigen or row.get("chain_type") != "AG"
    ]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=REPORT_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in filtered_rows:
            writer.writerow({
                field: row.get(field, "")
                for field in REPORT_FIELDS
            })