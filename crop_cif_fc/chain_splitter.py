"""chain_splitter.py

Detects when a single physical CIF chain actually contains more than one
Ig domain, so it can later be physically split into separate gemmi chains.

Both paths return the same normalized shape: a list of domains, each with
CIF-residue-list indices (positions in `residue_refs`, i.e. exactly what
get_chain_residues()/build_residue_refs() already return) - ready for the
splitting step to consume without caring which path found them.
"""
import gemmi

from .anarci_annotator import run_anarci
from .annotation import annotate_record
from .sequence_mapping import map_fasta_to_cif
from .cif_utils import get_chain_residues, residue_to_one_letter, build_residue_refs
from .models import ChainAnnotation


def build_cif_chain_sequence(residue_refs: list) -> str:
    chars = []
    for residue in residue_refs:
        info = gemmi.find_tabulated_residue(residue.name)
        chars.append(info.one_letter_code.upper() if info and info.one_letter_code.strip() else "X")
    return "".join(chars)


def detect_chain_domains(chain_id, residue_refs, fasta_records):
    """
    Returns (domains, source, status).

    domains: list of {"chain_type": str, "cif_start": int, "cif_end": int}
             start/end are inclusive indices into `residue_refs` (0-based) -
             directly usable, no further coordinate translation needed.
    source:  "matched" or "orphan", kept for the report / debugging.
    status:  ANARCI status string ("OK", "NO_DOMAIN", ...).
    """
    if chain_id in fasta_records:
        return _detect_from_fasta_reference(chain_id, residue_refs, fasta_records)
    return _detect_from_cif_sequence(residue_refs)


def _fasta_positions_to_cif_range(mapping: list[dict], fasta_start: int, fasta_end: int):
    """
    Translate an ANARCI domain range (0-based, inclusive, in FASTA-position
    space) into a 0-based inclusive range of CIF residue indices, using the
    existing FASTA<->CIF alignment mapping.
    """
    cif_indices = [
        mapping[pos]["cif_index"]
        for pos in range(fasta_start, fasta_end + 1)
        if mapping[pos]["in_cif"]
    ]

    if not cif_indices:
        # Whole domain fell in an unresolved stretch -- nothing to split on.
        return None

    return min(cif_indices), max(cif_indices)


def _detect_from_fasta_reference(chain_id, residue_refs, fasta_records):
    """
    'Matched' path: chain_id still has its own FASTA record. Run ANARCI on
    the FASTA reference sequence, then project each reported
    domain from FASTA-position space into CIF-index space so downstream
    splitting only ever has to think in terms of residue_refs positions.
    """
    record = fasta_records[chain_id]
    raw_domains, status = run_anarci(record.sequence)

    if not raw_domains:
        # Nothing to split; existing AG/manual_review
        # handling downstream takes over as before.
        return [], "matched", status

    refs = build_residue_refs(chain_id, residue_refs)
    mapping, _alignment = map_fasta_to_cif(record.sequence, refs)

    domains = []
    for chain_type, fasta_start, fasta_end in raw_domains:
        cif_range = _fasta_positions_to_cif_range(mapping, fasta_start, fasta_end)
        if cif_range is None:
            continue
        cif_start, cif_end = cif_range
        domains.append({
            "chain_type": chain_type,
            "cif_start": cif_start,
            "cif_end": cif_end,
        })

    return domains, "matched", status


def _detect_from_cif_sequence(residue_refs):
    """
    'Orphan' path: no FASTA record left for this physical chain_id (it was
    already pre-split externally, e.g. A -> C, D). Build the sequence
    straight from observed CIF residues and run ANARCI on that instead --
    positions returned by ANARCI are then ALREADY 0-based indices into
    residue_refs, no alignment/translation needed.
    """
    cif_sequence = "".join(residue_to_one_letter(r) for r in residue_refs)
    raw_domains, status = run_anarci(cif_sequence)

    if not raw_domains:
        return [], "orphan", status

    domains = [
        {"chain_type": chain_type, "cif_start": start, "cif_end": end}
        for chain_type, start, end in raw_domains
    ]

    return domains, "orphan", status


def resolve_unclaimed_residues(domains: list[dict]):
    """
    The function detects whether an interior gap
    existed, so the chain can be flagged for manual review.
    """
    if len(domains) < 2:
        return domains, False

    domains = sorted(domains, key=lambda d: d["cif_start"])

    flagged = any(
        domains[i + 1]["cif_start"] > domains[i]["cif_end"] + 1
        for i in range(len(domains) - 1)
    )

    return domains, flagged


def split_chain(structure: "gemmi.Structure", chain_id: str, domains: list[dict], new_chain_ids: list[str]) -> None:
    """
    Physically splits one physical chain into len(domains) new chains.

    domains:        output of Block 2/3 -- {"chain_type", "cif_start", "cif_end"},
                     cif_start/cif_end are 0-based inclusive indices into
                     residue_refs / get_chain_residues(model, chain_id) order.
    new_chain_ids:   pre-allocated new IDs, same length/order as domains
                     (domain i -> new_chain_ids[i]).
    """
    model = structure[0]
    source = model[chain_id]
    residues = get_chain_residues(structure, chain_id)

    for domain, new_id in zip(domains, new_chain_ids):
        new_chain = model.add_chain(new_id)
        for idx in range(domain["cif_start"], domain["cif_end"] + 1):
            residue = residues[idx]
            residue.subchain = new_id
            new_chain.add_residue(residue)
            
    claimed_seqids = {
        (residues[idx].seqid.num, residues[idx].seqid.icode)
        for domain in domains
        for idx in range(domain["cif_start"], domain["cif_end"] + 1)
    }
    for i in range(len(source) - 1, -1, -1):
        r = source[i]
        if (r.seqid.num, r.seqid.icode) in claimed_seqids:
            del source[i]

    model.remove_chain(chain_id)


def allocate_new_ids(used_chain_ids: set[str], n_needed: int) -> list[str]:
    """
    Ported from the external FASTA-splitter, adapted to hand back N ids at
    once instead of one pair. Same convention: next free single-letter
    chain_id after all currently used ones (A, B -> C, D, ...).
    Mutates used_chain_ids in place as it allocates, so a single shared
    set threaded through the whole structure prevents any collision
    across multiple chains needing a split in the same file.
    """
    new_ids = []
    letter_codes = [
        ord(cid) for cid in used_chain_ids if len(cid) == 1 and cid.isalpha()
    ]
    next_code = max(letter_codes, default=ord("@")) + 1

    while len(new_ids) < n_needed:
        while chr(next_code) in used_chain_ids:
            next_code += 1
        new_id = chr(next_code)
        used_chain_ids.add(new_id)
        new_ids.append(new_id)
        next_code += 1

    return new_ids


def build_split_chain_annotation(chain_type: str, new_chain_residues: list) -> tuple["ChainAnnotation", list[dict]]:
    """
    Simplified re-annotation for a freshly split chain

    Always builds the sequence straight from observed CIF residues,
    same as the orphan path -- gapless by construction, matched-path
    reference-sequence nuance deliberately deferred.
    """
    annotation = ChainAnnotation(
        chain_type=chain_type,
        range_start=0,
        range_end=len(new_chain_residues) - 1,
    )
    mapping = [
        {
            "fasta_position": i,
            "reference_aa": residue_to_one_letter(residue),
            "in_cif": True,
            "cif_index": i,
            "cif_residue": residue,
            "match": True,
        }
        for i, residue in enumerate(new_chain_residues)
    ]
    return annotation, mapping


def split_structure_chains(structure, fasta_records):
    model = structure[0]
    original_chain_ids = [chain.name for chain in model]
    used_chain_ids = set(original_chain_ids)

    split_report = {}
    new_chain_annotations = {}
    new_chain_mappings = {}

    for chain_id in original_chain_ids:
        residue_refs = get_chain_residues(structure, chain_id)
        if not residue_refs:
            continue

        domains, source, status = detect_chain_domains(chain_id, residue_refs, fasta_records)
        if len(domains) < 2:
            continue

        domains, flagged = resolve_unclaimed_residues(domains)
        new_chain_ids = allocate_new_ids(used_chain_ids, len(domains))
        split_chain(structure, chain_id, domains, new_chain_ids)

        for new_id, domain in zip(new_chain_ids, domains):
            new_residues = get_chain_residues(structure, new_id)
            annotation, mapping = build_split_chain_annotation(domain["chain_type"], new_residues)
            new_chain_annotations[new_id] = annotation
            new_chain_mappings[new_id] = mapping

        split_report[chain_id] = {
            "new_chain_ids": new_chain_ids,
            "domain_types": [d["chain_type"] for d in domains],
            "flagged": flagged,
            "source": source,
            "status": status,
        }

    return split_report, new_chain_annotations, new_chain_mappings
