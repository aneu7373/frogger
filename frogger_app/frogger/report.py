def build_report_row(construct, result, codonopt_version, seed_codonopt, seed_reassembly):
    return {
        "construct_id": construct["construct_id"],
        "frag_sources": construct["frag_sources_str"],
        "len_nt": len(construct["core_orf_seq"]),
        "len_aa": result["aa_len"],
        "passes_all": result["passes_all"],
        "fail_reasons": ";".join(result["fail_reasons"]),
        "gc_fraction": f"{result['gc_fraction']:.5f}",
        "max_homopolymer": result["max_homopolymer"],
        "n_hits": len(result["hits"]),
        "codonopt_version": codonopt_version,
        "seed_codonopt": int(seed_codonopt),
        "seed_reassembly": int(seed_reassembly),
    }


def build_hit_rows(construct, result):
    rows = []
    for h in result["hits"]:
        rows.append(
            {
                "construct_id": construct["construct_id"],
                "frag_sources": construct["frag_sources_str"],
                "motif": h["motif"],
                "start": h["start"],
                "end": h["end"],
                "context": h["context"],
            }
        )
    return rows
