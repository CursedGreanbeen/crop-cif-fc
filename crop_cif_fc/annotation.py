from __future__ import annotations

from pathlib import Path
from typing import Any

from .anarci_annotator import chain_role, run_anarci
from .fasta_manager import FastaRecord, read_fasta_file


def annotate_record(record: FastaRecord) -> dict[str, Any]:
    """Запускает ANARCI для одной FASTA-записи."""
    domains, status = run_anarci(record.sequence)

    chain_type, start, end = domains[0] if domains else (None, None, None)
    role = chain_role(chain_type)

    return {
        "header": record.header,
        "sequence": record.sequence,
        "chain_ids": record.chain_ids,
        "chain_type": chain_type,
        "role": role,
        "start": start,
        "end": end,
        "status": status,
        "domains": domains,
    }


def expand_annotations_by_chain(
    annotations: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """
    Делает отдельную запись для каждой цепи.
    """
    expanded: list[dict[str, Any]] = []

    for annotation in annotations:
        for chain_id in annotation["chain_ids"]:
            row = dict(annotation)
            row["chain_id"] = chain_id
            expanded.append(row)

    return expanded


