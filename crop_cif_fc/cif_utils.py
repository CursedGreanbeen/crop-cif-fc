from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import re
import subprocess, sys, site
import gemmi

from .models import (
    ResidueRef,
    ChainAnnotation,
    ChainCropPlan
)


def load_structure(cif_path: str | Path) -> gemmi.Structure:
    """Читает CIF и возвращает структуру только с первой моделью."""
    structure = gemmi.read_structure(str(cif_path))
    structure.setup_entities()
    while len(structure) > 1:
        del structure[1]
    return structure


def save_structure(
    structure: gemmi.Structure,
    path: str | Path,
) -> None:
    """Save a Gemmi structure as an mmCIF file."""

    output_path = Path(path)
    document = structure.make_mmcif_document()
    document.write_file(str(output_path))


def get_protein_chain_ids(
    structure: gemmi.Structure,
) -> list[str]:
    """Return IDs of chains containing at least one protein residue."""

    if not structure:
        raise ValueError("Structure contains no models")

    model = structure[0]
    chain_ids = []

    for chain in model:
        if any(
            residue.het_flag == "A"
            and gemmi.find_tabulated_residue(residue.name).is_amino_acid()
            for residue in chain
        ):
            chain_ids.append(chain.name)

    return chain_ids


def get_chain_residues(
    structure: gemmi.Structure,
    chain_id: str,
) -> list[gemmi.Residue]:
    """Return ordered, observed protein residues from one chain."""

    if not structure:
        raise ValueError("Structure contains no models")

    model = structure[0]

    for chain in model:
        if chain.name != chain_id:
            continue

        residues = []

        for residue in chain:
            if residue.het_flag != "A":
                continue

            residue_info = gemmi.find_tabulated_residue(residue.name)

            if residue_info.is_amino_acid():
                residues.append(residue)

        return residues

    raise KeyError(f"Chain not found: {chain_id}")


def residue_to_one_letter(residue: gemmi.Residue) -> str:
    """Конвертирует остаток (включая нестандартные) в one-letter код. Неизвестные -> 'X'."""
    info = gemmi.find_tabulated_residue(residue.name)
    if info is None:
        return 'X'
    code = info.one_letter_code
    if not code or code.strip() == '':
        return 'X'
    return code.upper()


def apply_structure_plans(
    structure: gemmi.Structure,
    plans: dict[str, ChainCropPlan],
) -> gemmi.Structure:
    """Apply chain crop plans to the first model in place."""

    if not structure:
        raise ValueError("Structure contains no models")

    model = structure[0]

    for chain in model:
        plan = plans.get(chain.name)

        if plan is None:
            continue

        if plan.action in {"keep_all", "manual_review"}:
            continue

        if plan.action != "crop":
            raise ValueError(f"Unknown plan action: {plan.action}")

        residues = get_chain_residues(structure, chain.name)
        keep_indices = set(plan.keep_cif_indices)
    
        drop_seqids = {
            (residue.seqid.num, residue.seqid.icode)
            for cif_index, residue in enumerate(residues)
            if cif_index not in keep_indices
        }
    
        for i in range(len(chain) - 1, -1, -1):
            real_residue = chain[i]
            if (real_residue.seqid.num, real_residue.seqid.icode) in drop_seqids:
                del chain[i]
    
    return structure


def build_residue_refs(
    chain_id: str,
    residues: list[gemmi.Residue],
) -> list[ResidueRef]:
    """Build immutable references for observed residues in CIF order."""

    return [
        ResidueRef(
            chain_id=chain_id,
            cif_index=index,
            name=residue.name,
            auth_seq_id=residue.seqid.num,
            insertion_code=residue.seqid.icode,
            one_letter=residue_to_one_letter(residue)
        )
        for index, residue in enumerate(residues)
    ]