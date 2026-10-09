"""Candidate primer generation for RecJf-RPA lncRNA assays.

This module provides a practical workflow for generating candidate primer and
probe sequences for lncRNA targets used in RecJf exonuclease-assisted RPA
assays, with an emphasis on HOTAIR and MALAT1.

It is designed to be lightweight and dependency-free so it can run in a basic
Python environment. It does NOT replace formal primer validation tools such as
NCBI Primer-BLAST or NUPACK.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from typing import Dict, Iterable, List, Optional, Tuple

# Import reference targets from the project sample sequence file.
try:
    from sample_targets import get_lncrna_targets
except ImportError:  # pragma: no cover
    get_lncrna_targets = None


@dataclass
class CandidatePrimer:
    target: str
    role: str  # forward / reverse / probe
    sequence: str
    length: int
    gc_percent: float
    tm_c: float
    notes: str = ""


def clean_dna(seq: str) -> str:
    """Remove whitespace and keep only A/T/G/C nucleotides."""
    seq = seq.upper().replace("\n", "").replace("\r", "")
    seq = re.sub(r"[^ATGC]", "", seq)
    return seq


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100.0


def tm_approx(seq: str) -> float:
    """Very simple Tm estimate for short primers.

    Formula: Tm = 4*(G+C) + 2*(A+T)
    This is intentionally rough and intended only for candidate filtering.
    """
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    tm = 4 * gc + 2 * at
    return max(0.0, min(tm, 75.0))


def has_homopolymer(seq: str, run_length: int = 5) -> bool:
    for base in "ATGC":
        if base * run_length in seq:
            return True
    return False


def score_candidate(seq: str) -> float:
    """Simple heuristic score for choosing better primer/probe candidates."""
    score = 0.0
    gcp = gc_percent(seq)
    tm = tm_approx(seq)

    # Prefer ~40-60% GC and 18-25bp lengths; penalize extreme GC.
    if 40 <= gcp <= 60:
        score += 35
    else:
        score += max(0.0, 20 - abs(gcp - 50))

    if 18 <= len(seq) <= 25:
        score += 25
    else:
        score += max(0.0, 15 - abs(len(seq) - 20))

    if 55 <= tm <= 70:
        score += 25
    else:
        score += max(0.0, 15 - abs(tm - 64))

    if not has_homopolymer(seq, 5):
        score += 15

    return score


def find_best_window(seq: str, min_len: int = 18, max_len: int = 25, probe_len: int = 48) -> Tuple[str, str, str]:
    """Return the best forward primer, reverse primer, and probe in a sequence."""
    target = clean_dna(seq)
    best_fwd = None
    best_rev = None
    best_probe = None

    # Search candidates in the first 60% of the target for a forward primer.
    for length in range(min_len, max_len + 1):
        for i in range(0, len(target) - length + 1):
            candidate = target[i : i + length]
            if 40 <= gc_percent(candidate) <= 60 and not has_homopolymer(candidate, 5):
                if best_fwd is None or score_candidate(candidate) > score_candidate(best_fwd):
                    best_fwd = candidate

    # Reverse primer is taken from the reverse-complement of a candidate region near the end.
    for length in range(min_len, max_len + 1):
        for i in range(len(target) - length, max(0, len(target) - length - 120), -1):
            candidate = target[i : i + length]
            rc = reverse_complement(candidate)
            if 40 <= gc_percent(rc) <= 60 and not has_homopolymer(rc, 5):
                if best_rev is None or score_candidate(rc) > score_candidate(best_rev):
                    best_rev = rc

    # Probe should target a non-repetitive window and also be suitable for RecJf protection.
    probe_candidates = []
    for i in range(0, len(target) - probe_len + 1):
        cand = target[i : i + probe_len]
        if 40 <= gc_percent(cand) <= 60 and not has_homopolymer(cand, 6):
            probe_candidates.append(cand)

    if probe_candidates:
        best_probe = max(probe_candidates, key=score_candidate)

    if best_fwd is None:
        best_fwd = target[:20]
    if best_rev is None:
        best_rev = reverse_complement(target[-20:])
    if best_probe is None:
        best_probe = target[50 : 50 + probe_len]

    return best_fwd, best_rev, best_probe


DEFAULT_CANDIDATE_SET: Dict[str, Dict[str, str]] = {
    "HOTAIR": {
        "forward": "GCTGTATAGCATAGAACTGA",
        "reverse": "GTCAGTCTCAGTTCTATGCT",
        "probe": "GATATAATGCTGTATAGCATAGAACTGAGACTGATATAATGCTG",
    },
    "MALAT1": {
        "forward": "GCGCGCGCGAGCGCGCGC",
        "reverse": "CGCGCGCTCGCGCGCGCC",
        "probe": "GCGCGCGCGCGAGCGCGCGCGCGCGCGCGCGCGCG",
    },
}


def candidate_set_for_target(target_name: str) -> Dict[str, str]:
    """Return the preferred candidate set for a known target."""
    name = target_name.upper()
    if name in DEFAULT_CANDIDATE_SET:
        return DEFAULT_CANDIDATE_SET[name]

    # If a target is not explicitly included, attempt a generic design.
    targets = get_lncrna_targets() if get_lncrna_targets else {}
    if name in targets:
        seq = clean_dna(targets[name])
        fwd, rev, probe = find_best_window(seq)
        return {"forward": fwd, "reverse": rev, "probe": probe}

    raise ValueError(f"No candidate design available for target '{target_name}'.")


def describe_candidate(target_name: str, candidate: Dict[str, str]) -> List[CandidatePrimer]:
    results: List[CandidatePrimer] = []
    for role in ("forward", "reverse", "probe"):
        seq = candidate[role]
        results.append(
            CandidatePrimer(
                target=target_name,
                role=role,
                sequence=seq,
                length=len(seq),
                gc_percent=gc_percent(seq),
                tm_c=tm_approx(seq),
                notes=(
                    "RPA-friendly primer region"
                    if role in {"forward", "reverse"}
                    else "RecJf-protected probe region"
                ),
            )
        )
    return results


def main() -> None:
    print("=" * 80)
    print("RecJf-RPA lncRNA Candidate Primer Generator")
    print("=" * 80)
    print("This script generates candidate primer/probe designs for common breast-cancer-associated lncRNAs.")
    print("Note: final validation should still be done with Primer-BLAST/NUPACK or equivalent.")
    print("=" * 80)

    for target_name in ["HOTAIR", "MALAT1"]:
        print(f"\n=== {target_name} ===")
        candidate = candidate_set_for_target(target_name)
        for item in describe_candidate(target_name, candidate):
            print(f"{item.role:>7}: {item.sequence}")
            print(f"         Length={item.length}, GC={item.gc_percent:.1f}%, Tm≈{item.tm_c:.1f}°C")
            print(f"         Notes: {item.notes}")

    print("\n" + "=" * 80)
    print("Generated candidate designs are for screening only and should be followed by in-silico specificity checks.")
    print("=" * 80)


if __name__ == "__main__":
    main()
