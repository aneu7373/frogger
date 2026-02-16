from collections import defaultdict
import random
import hashlib


def generate_constructs(
    fragments,
    max_constructs: int,
    sample_n: int,
    seed: int,
    deduplicate: bool = True,
):
    """
    Reassemble core fragments into post-ligation ORFs.
    Maintains fragment order: position 1 then 2 then ... k.
    Allows mixing sources across genes, within each fragment position.

    Deterministic sampling when combinations exceed max_constructs.
    """
    by_pos = defaultdict(list)
    for f in fragments:
        by_pos[int(f["frag_index"])].append(f)

    positions = sorted(by_pos.keys())
    if not positions:
        return []

    total = 1
    for p in positions:
        total *= len(by_pos[p])

    rng = random.Random(int(seed))
    target = total if total <= int(max_constructs) else min(int(sample_n), int(max_constructs))

    constructs = []
    seen = set()

    def prov_tuple(choice):
        return tuple((c["gene_id"], int(c["frag_index"])) for c in choice)

    tries = 0
    retry_cap = target * 50

    while len(constructs) < target and tries < retry_cap:
        tries += 1
        choice = [rng.choice(by_pos[p]) for p in positions]
        prov = prov_tuple(choice)
        if deduplicate and prov in seen:
            continue
        seen.add(prov)

        core = "".join(c["core_seq"] for c in choice)
        prov_str = ",".join(f"{g}:{i}" for g, i in prov)
        cid = hashlib.sha256(prov_str.encode("utf-8")).hexdigest()[:12]

        constructs.append(
            {
                "construct_id": f"c_{cid}",
                "frag_sources": prov,
                "frag_sources_str": prov_str,
                "core_orf_seq": core,
            }
        )

    return constructs
