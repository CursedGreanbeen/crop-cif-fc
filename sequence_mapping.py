from Bio.Align import PairwiseAligner
from Bio.Align import Alignment

from .cif_utils import residue_to_one_letter
from .models import ResidueRef


def map_fasta_to_cif(
    reference_sequence: str,
    cif_residues: list[ResidueRef],
) -> tuple[list[dict], object]:
    """
    Map positions of a reference FASTA sequence to observed CIF residues.

    Parameters
    ----------
    reference_sequence
        Full continuous sequence from FASTA.

    cif_residues
        Observed CIF residues in chain order, as ResidueRef.

    Returns
    -------
    mapping
        One record for every FASTA position.

    alignment
        Selected global alignment, retained for diagnostics.
    """

    reference_sequence = reference_sequence.strip().upper()

    if not reference_sequence:
        raise ValueError("Reference sequence is empty")

    if not cif_residues:
        raise ValueError("CIF residue list is empty")

    cif_sequence = "".join(ref.one_letter for ref in cif_residues)

    if "-" in reference_sequence or "-" in cif_sequence:
        raise ValueError(
            "Sequences must not contain gap characters; "
            "missing CIF residues should be absent from cif_sequence"
        )

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(reference_sequence, cif_sequence)

    if len(alignments) == 0:
        raise ValueError("Unable to align reference and CIF sequences")

    alignment = alignments[0]

    mapping = build_position_mapping(
        reference_sequence=reference_sequence,
        cif_sequence=cif_sequence,
        alignment=alignment,
        cif_residues=cif_residues,
    )

    return mapping, alignment


def build_position_mapping(
    reference_sequence: str,
    cif_sequence: str,
    alignment,
    cif_residues: list,
) -> list[dict]:
    """
    Map every reference FASTA position to a CIF residue or None.

    cif_residues must be in exactly the same order as cif_sequence.
    """

    reference_sequence = reference_sequence.strip().upper()
    cif_sequence = cif_sequence.strip().upper()

    if len(cif_sequence) != len(cif_residues):
        raise ValueError(
            "cif_residues length must match cif_sequence length"
        )

    mapping = [
        {
            "fasta_position": fasta_index,
            "reference_aa": reference_sequence[fasta_index],
            "in_cif": False,
            "cif_index": None,
            "cif_residue": None,
            "match": None,
        }
        for fasta_index in range(len(reference_sequence))
    ]

    # aligned contains only residue-to-residue alignment blocks.
    # Gaps are therefore left as unmapped positions.
    for (ref_start, ref_end), (cif_start, cif_end) in zip(
        alignment.aligned[0],
        alignment.aligned[1],
    ):
        block_length = ref_end - ref_start

        if block_length != cif_end - cif_start:
            raise ValueError(
                "Invalid alignment block: reference and CIF block lengths differ"
            )

        for offset in range(block_length):
            fasta_index = ref_start + offset
            cif_index = cif_start + offset

            reference_aa = reference_sequence[fasta_index]
            cif_aa = cif_sequence[cif_index]

            mapping[fasta_index].update(
                {
                    "in_cif": True,
                    "cif_index": cif_index,
                    "cif_residue": cif_residues[cif_index],
                    "match": reference_aa == cif_aa,
                }
            )

    return mapping

