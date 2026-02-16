from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Iterable, List, Optional, Sequence, Tuple


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_ORDER)

CODON_TO_AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAA": "*", "TAG": "*", "TGA": "*",
}

MUT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


@dataclass(frozen=True)
class Variant:
    gene_id: str
    pos1: int
    wt_aa: str
    target_aa: str
    mutation_label: str
    codon_used: str
    protein: str
    cds: str


def _normalize_dna(seq: str) -> str:
    return seq.upper().replace("U", "T").replace(" ", "").replace("\n", "").replace("\r", "")


def translate_cds(cds: str) -> str:
    cds = _normalize_dna(cds)
    if len(cds) % 3 != 0:
        raise ValueError("CDS length is not a multiple of 3.")
    prot = []
    for i in range(0, len(cds), 3):
        codon = cds[i:i + 3]
        aa = CODON_TO_AA.get(codon)
        if aa is None:
            raise ValueError(f"Invalid/unknown codon: {codon}")
        prot.append(aa)
    return "".join(prot)


def _contains_any_motif(seq: str, motifs: Sequence[str]) -> bool:
    s = seq.upper()
    for m in motifs:
        mm = str(m).upper()
        if mm and mm in s:
            return True
    return False


def _as_list_of_codon_freq_pairs(x: Any) -> List[Tuple[str, float]]:
    """
    Normalize various codon table representations to List[(codon, freq)].
    Supports frogger.common.codon_table.CodonEntry objects (codon, fraction).
    """
    if x is None:
        return []

    if isinstance(x, dict):
        out = [(str(cod), float(freq)) for cod, freq in x.items()]
        out.sort(key=lambda t: t[1], reverse=True)
        return out

    if isinstance(x, list):
        if not x:
            return []
        first = x[0]

        if isinstance(first, tuple) and len(first) == 2:
            return [(str(c), float(f)) for (c, f) in x]

        # frogger CodonEntry-like
        if hasattr(first, "codon") and (hasattr(first, "fraction") or hasattr(first, "freq")):
            out = []
            for e in x:
                cod = getattr(e, "codon")
                frac = getattr(e, "fraction", getattr(e, "freq", 0.0))
                out.append((str(cod), float(frac)))
            out.sort(key=lambda t: t[1], reverse=True)
            return out

        return [(str(c), 0.0) for c in x]

    return []


def _get_codon_pairs_for_aa(codon_table: Any, aa: str) -> List[Tuple[str, float]]:
    """
    Supports:
      1) dict-like: codon_table[aa] -> list[(codon,freq)] or other representations
      2) object with method get_codons(aa) or codons_for_aa(aa)
      3) object with mapping attribute aa_to_codons or table
    """
    aa = aa.upper()

    # dict-like
    try:
        val = codon_table[aa]  # type: ignore[index]
        pairs = _as_list_of_codon_freq_pairs(val)
        if pairs:
            return pairs
    except Exception:
        pass

    # method patterns
    for meth_name in ("get_codons", "codons_for_aa", "get_for_aa"):
        if hasattr(codon_table, meth_name):
            meth = getattr(codon_table, meth_name)
            try:
                val = meth(aa)
                pairs = _as_list_of_codon_freq_pairs(val)
                if pairs:
                    return pairs
            except Exception:
                pass

    # attribute patterns
    for attr in ("aa_to_codons", "table", "mapping"):
        if hasattr(codon_table, attr):
            m = getattr(codon_table, attr)
            try:
                val = m[aa]
                pairs = _as_list_of_codon_freq_pairs(val)
                if pairs:
                    return pairs
            except Exception:
                pass

    raise TypeError(
        "Unsupported codon_table type for mutation codon selection. "
        f"type={type(codon_table)}. Expected dict-like mapping or CodonTable exposing "
        "one of: __getitem__[aa], get_codons(aa), codons_for_aa(aa), aa_to_codons[aa], table[aa]. "
        f"Available attrs: {sorted([a for a in dir(codon_table) if not a.startswith('_')])}"
    )


def _codon_choices_for_aa(target_aa: str, codon_table: Any) -> List[str]:
    """
    Returns codons in preferred order (highest frequency first).
    """
    aa = target_aa.upper()
    if aa not in AA_SET:
        raise ValueError(f"Invalid amino acid '{aa}' for codon selection.")

    pairs = _get_codon_pairs_for_aa(codon_table, aa)
    if not pairs:
        raise ValueError(f"No codon choices found for amino acid '{aa}' in codon_table.")

    out = [_normalize_dna(cod) for cod, _freq in pairs]

    # de-dup preserve order
    seen = set()
    uniq = []
    for c in out:
        if c not in seen:
            uniq.append(c)
            seen.add(c)
    return uniq


def choose_mutation_codon(
    baseline_cds: str,
    codon_index0: int,
    target_aa: str,
    forbidden_motifs: Sequence[str],
    codon_table: Any,
) -> Tuple[str, str]:
    baseline_cds = _normalize_dna(baseline_cds)
    start = 3 * codon_index0
    end = start + 3
    if start < 0 or end > len(baseline_cds):
        raise IndexError("codon_index0 is out of bounds for baseline_cds.")

    choices = _codon_choices_for_aa(target_aa, codon_table)
    for codon in choices:
        candidate = baseline_cds[:start] + codon + baseline_cds[end:]
        if not _contains_any_motif(candidate, forbidden_motifs):
            return codon, candidate

    # fallback: best codon even if it introduces a forbidden motif
    top = choices[0]
    candidate = baseline_cds[:start] + top + baseline_cds[end:]
    return top, candidate


def parse_explicit_mutation_string(mut: str) -> Tuple[str, int, str]:
    m = MUT_RE.match(mut.strip())
    if not m:
        raise ValueError(f"Invalid mutation string '{mut}'. Expected format like 'T24V'.")

    wt = m.group(1).upper()
    pos1 = int(m.group(2))
    tgt = m.group(3).upper()

    if wt not in AA_SET:
        raise ValueError(f"Invalid WT amino acid '{wt}' in '{mut}'. Must be one of {AA_ORDER}.")
    if tgt not in AA_SET:
        raise ValueError(f"Invalid target amino acid '{tgt}' in '{mut}'. Must be one of {AA_ORDER}.")
    if pos1 < 1:
        raise ValueError(f"Invalid position {pos1} in '{mut}'. Positions are 1-based and must be >= 1.")

    return wt, pos1, tgt


def _normalize_positions_to_mutate(
    aa_len: int,
    positions_to_mutate: Optional[Sequence[int]],
    *,
    allow_start_mutation: bool,
) -> List[int]:
    """
    Normalize mutation positions (1-based) and enforce default policy:
      - Full GSSM default excludes position 1 (start codon).
      - If caller passes explicit positions, we still exclude position 1 unless allow_start_mutation=True.

    Returns sorted unique positions.
    """
    if positions_to_mutate is None:
        positions = list(range(2, aa_len + 1)) if aa_len >= 2 else []
    else:
        positions = [int(p) for p in positions_to_mutate]

    # Validate bounds first (helps user pinpoint bad configs)
    for p in positions:
        if p < 1 or p > aa_len:
            raise ValueError(f"Requested mutation position {p} outside protein length {aa_len}.")

    if not allow_start_mutation:
        positions = [p for p in positions if p != 1]

    return sorted(set(positions))


def generate_saturation_variants(
    *,
    gene_id: str,
    baseline_cds: str,
    forbidden_motifs: Sequence[str],
    codon_table: Any,
    include_wt: bool = True,
    positions_to_mutate: Optional[Sequence[int]] = None,
    allow_start_mutation: bool = False,
) -> Iterable[Variant]:
    baseline_cds = _normalize_dna(baseline_cds)
    baseline_prot = translate_cds(baseline_cds)
    aa_len = len(baseline_prot)

    # Cleanest architecture: enforce "do not mutate position 1" policy HERE,
    # regardless of what the pipeline passes, unless explicitly overridden.
    positions = _normalize_positions_to_mutate(
        aa_len=aa_len,
        positions_to_mutate=positions_to_mutate,
        allow_start_mutation=allow_start_mutation,
    )

    for pos1 in positions:
        wt = baseline_prot[pos1 - 1]
        codon_index0 = pos1 - 1

        if include_wt:
            mut_label = f"{wt}{pos1}{wt}"
            yield Variant(
                gene_id=gene_id,
                pos1=pos1,
                wt_aa=wt,
                target_aa=wt,
                mutation_label=mut_label,
                codon_used="WT",
                protein=baseline_prot,
                cds=baseline_cds,
            )

        for target in AA_ORDER:
            if target == wt:
                continue
            codon_used, mutated_cds = choose_mutation_codon(
                baseline_cds=baseline_cds,
                codon_index0=codon_index0,
                target_aa=target,
                forbidden_motifs=forbidden_motifs,
                codon_table=codon_table,
            )
            mutated_prot = translate_cds(mutated_cds)
            mut_label = f"{wt}{pos1}{target}"
            yield Variant(
                gene_id=gene_id,
                pos1=pos1,
                wt_aa=wt,
                target_aa=target,
                mutation_label=mut_label,
                codon_used=codon_used,
                protein=mutated_prot,
                cds=mutated_cds,
            )


def generate_explicit_variants(
    *,
    gene_id: str,
    baseline_cds: str,
    forbidden_motifs: Sequence[str],
    codon_table: Any,
    explicit_mutations: Sequence[str],
    require_wt_match: bool = True,
) -> Iterable[Variant]:
    baseline_cds = _normalize_dna(baseline_cds)
    baseline_prot = translate_cds(baseline_cds)
    aa_len = len(baseline_prot)

    seen = set()
    muts: List[str] = []
    for m in explicit_mutations:
        mm = str(m).strip()
        if mm and mm not in seen:
            muts.append(mm)
            seen.add(mm)

    for mstr in muts:
        wt, pos1, tgt = parse_explicit_mutation_string(mstr)
        if pos1 > aa_len:
            raise ValueError(f"{gene_id}: explicit mutation {mstr} exceeds protein length {aa_len}.")

        # Explicit mutations are explicit: allow pos1==1 if user requested it.
        baseline_wt = baseline_prot[pos1 - 1]
        if require_wt_match and baseline_wt != wt:
            raise ValueError(
                f"{gene_id}: explicit mutation {mstr} WT='{wt}' does not match baseline AA at position {pos1} "
                f"(baseline has '{baseline_wt}')."
            )

        codon_index0 = pos1 - 1

        if tgt == baseline_wt:
            mut_label = f"{baseline_wt}{pos1}{baseline_wt}"
            yield Variant(
                gene_id=gene_id,
                pos1=pos1,
                wt_aa=baseline_wt,
                target_aa=baseline_wt,
                mutation_label=mut_label,
                codon_used="WT",
                protein=baseline_prot,
                cds=baseline_cds,
            )
            continue

        codon_used, mutated_cds = choose_mutation_codon(
            baseline_cds=baseline_cds,
            codon_index0=codon_index0,
            target_aa=tgt,
            forbidden_motifs=forbidden_motifs,
            codon_table=codon_table,
        )
        mutated_prot = translate_cds(mutated_cds)
        mut_label = f"{baseline_wt}{pos1}{tgt}"
        yield Variant(
            gene_id=gene_id,
            pos1=pos1,
            wt_aa=baseline_wt,
            target_aa=tgt,
            mutation_label=mut_label,
            codon_used=codon_used,
            protein=mutated_prot,
            cds=mutated_cds,
        )
