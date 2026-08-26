from .models import (
    ChainAnnotation,
    ChainCropPlan
)


def plan_chain_crop(
    annotation: ChainAnnotation,
    mapping: list[dict],
) -> ChainCropPlan:
    chain_type = annotation.chain_type

    if chain_type == "AG":
        return ChainCropPlan(
            action="keep_all",
            status="antigen_unchanged",
            chain_type=chain_type,
        )

    if chain_type not in {"H", "K", "L"}:
        return ChainCropPlan(
            action="manual_review",
            status="unknown_chain_type",
            chain_type=chain_type,
            reason="Chain is not recognized as H, K, L, or AG",
        )

    start = annotation.range_start
    end = annotation.range_end

    if start is None or end is None:
        return ChainCropPlan(
            action="manual_review",
            status="missing_domain_range",
            chain_type=chain_type,
            reason="No usable VH/VL range",
        )

    if not 0 <= start <= end < len(mapping):
        return ChainCropPlan(
            action="manual_review",
            status="invalid_domain_range",
            chain_type=chain_type,
            fasta_start=start,
            fasta_end=end,
            reason="Domain range is outside the FASTA sequence",
        )

    selected = mapping[start:end + 1]

    missing_positions = tuple(
        item["fasta_position"]
        for item in selected
        if not item["in_cif"]
    )

    boundary_missing = (
        not mapping[start]["in_cif"]
        or not mapping[end]["in_cif"]
    )

    keep_cif_indices = tuple(
        item["cif_index"]
        for item in mapping[start:end + 1]
        if item["in_cif"]
    )

    if boundary_missing:
        action = "manual_review"
        status = "boundary_residue_missing"
        reason = "Domain boundary residue is absent from CIF"
    else:
        action = "crop"
        status = "ready_for_crop"
        reason = None

    return ChainCropPlan(
        action=action,
        status=status,
        chain_type=chain_type,
        fasta_start=start,
        fasta_end=end,
        missing_positions=missing_positions,
        keep_cif_indices=keep_cif_indices,
        reason=reason,
    )


def plan_structure_crop(
    chain_annotations: dict[str, ChainAnnotation],
    chain_mappings: dict[str, list[dict]],
) -> dict[str, ChainCropPlan]:
    """
    Build a crop plan for every chain.

    Chains present only in CIF are marked for manual review.
    Chains present only in annotations are also marked for manual review.
    """

    plans: dict[str, ChainCropPlan] = {}

    all_chain_ids = set(chain_annotations) | set(chain_mappings)

    for chain_id in sorted(all_chain_ids):
        annotation = chain_annotations.get(chain_id)
        mapping = chain_mappings.get(chain_id)

        if annotation is None:
            plans[chain_id] = ChainCropPlan(
                action="manual_review",
                status="missing_annotation",
                chain_type="UNKNOWN",
                reason="Chain has no FASTA/annotation record",
            )
            continue

        if mapping is None:
            plans[chain_id] = ChainCropPlan(
                action="manual_review",
                status="missing_cif_mapping",
                chain_type=annotation.chain_type,
                    reason="Chain has no CIF sequence mapping",
            )
            continue

        plans[chain_id] = plan_chain_crop(annotation, mapping)

    return plans