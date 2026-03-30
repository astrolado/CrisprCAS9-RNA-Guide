from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import gzip


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def _clean_chrom_name(name: str) -> str:
    name = name.strip()
    if not name:
        return name
    if name.startswith(">"):
        name = name[1:]
    return name.split()[0]


def _canonical_chrom_name(name: str) -> str:
    """
    Return a standardized chromosome name for alias handling.

    Examples:
    - chr17 -> chr17
    - 17 -> chr17
    - X -> chrX
    - MT -> chrM
    - chrMT -> chrM
    - M -> chrM
    """
    name = _clean_chrom_name(name)
    if not name:
        return name

    lower = name.lower()
    if lower.startswith("chr"):
        base = name[3:]
    else:
        base = name

    base_upper = base.upper()

    if base_upper in {"M", "MT"}:
        return "chrM"

    return f"chr{base}"


def _build_chrom_aliases(raw_name: str) -> List[str]:
    """
    Generate all common aliases for a chromosome/contig name.
    """
    raw_name = _clean_chrom_name(raw_name)
    canonical = _canonical_chrom_name(raw_name)

    aliases = {raw_name, canonical}

    if canonical.startswith("chr"):
        no_chr = canonical[3:]
        aliases.add(no_chr)
        aliases.add("chr" + no_chr)

        no_chr_upper = no_chr.upper()
        aliases.add(no_chr_upper)
        aliases.add("chr" + no_chr_upper)

    if canonical == "chrM":
        aliases.update({"M", "MT", "chrM", "chrMT"})

    return [a for a in aliases if a]


class GenomeStore:
    """
    Simple in-memory genome store.

    Features:
    - load FASTA
    - fetch genomic intervals
    - cache whole chromosomes
    - exact and 1-mismatch PAM-aware guide counting for SpCas9-like systems
    - robust chromosome alias handling (chr17 <-> 17, chrM <-> MT, etc.)

    Assumptions:
    - guide length is 20
    - PAM is 3' of protospacer on target strand
    - SpCas9 PAM default is NGG
    """

    def __init__(self, fasta_path: str | Path) -> None:
        self.fasta_path = Path(fasta_path)
        self.chrom_sequences: Dict[str, str] = {}
        self.chrom_aliases: Dict[str, str] = {}
        self._off_target_cache: Dict[Tuple[str, Tuple[str, ...]], Dict[str, int | bool]] = {}

    def _register_chromosome(self, raw_name: str, sequence: str) -> None:
        """
        Store sequence under its original FASTA contig name and register aliases.
        """
        raw_name = _clean_chrom_name(raw_name)
        if not raw_name:
            return

        self.chrom_sequences[raw_name] = sequence.upper()

        for alias in _build_chrom_aliases(raw_name):
            self.chrom_aliases[alias] = raw_name

    def resolve_chrom_name(self, chrom: str) -> str:
        chrom = _clean_chrom_name(chrom)

        # exact raw name first
        if chrom in self.chrom_sequences:
            return chrom

        # alias match
        resolved = self.chrom_aliases.get(chrom)
        if resolved is not None:
            return resolved

        # canonical fallback
        canonical = _canonical_chrom_name(chrom)
        resolved = self.chrom_aliases.get(canonical)
        if resolved is not None:
            return resolved

        available_preview = ", ".join(list(self.chrom_sequences.keys())[:10])
        raise KeyError(f"Chromosome not found: {chrom}. Available examples: {available_preview}")

    def load(self) -> None:
        current_name: Optional[str] = None
        chunks: List[str] = []

        opener = gzip.open if self.fasta_path.suffix == ".gz" else open
        with opener(self.fasta_path, "rt") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if current_name is not None:
                        self._register_chromosome(current_name, "".join(chunks))
                    current_name = _clean_chrom_name(line[1:])
                    chunks = []
                else:
                    chunks.append(line)

        if current_name is not None:
            self._register_chromosome(current_name, "".join(chunks))

        if not self.chrom_sequences:
            raise ValueError(f"No chromosome sequences loaded from {self.fasta_path}")

    def get_sequence(self, chrom: str, start: int, end: int) -> str:
        """
        Fetch sequence using 1-based inclusive genomic coordinates.
        Accepts chr17 or 17 style chromosome names.
        """
        resolved_chrom = self.resolve_chrom_name(chrom)
        seq = self.chrom_sequences[resolved_chrom]

        start = max(1, start)
        end = min(len(seq), end)

        if end < start:
            return ""

        # 1-based inclusive -> Python slice
        return seq[start - 1:end]

    def chromosome_names(self) -> List[str]:
        return list(self.chrom_sequences.keys())

    def chromosome_length(self, chrom: str) -> int:
        resolved_chrom = self.resolve_chrom_name(chrom)
        return len(self.chrom_sequences[resolved_chrom])

    @staticmethod
    def _pam_matches(pam_seq: str, pam_patterns: Sequence[str]) -> bool:
        pam_seq = pam_seq.upper()
        for pattern in pam_patterns:
            pattern = pattern.upper()
            if len(pattern) != len(pam_seq):
                continue
            ok = True
            for base, pat in zip(pam_seq, pattern):
                if pat == "N":
                    continue
                if base != pat:
                    ok = False
                    break
            if ok:
                return True
        return False

    @staticmethod
    def _mismatch_count(a: str, b: str, max_mismatches: Optional[int] = None) -> int:
        mismatches = 0
        for x, y in zip(a, b):
            if x != y:
                mismatches += 1
                if max_mismatches is not None and mismatches > max_mismatches:
                    return mismatches
        return mismatches

    def count_guide_offtargets(
        self,
        guide: str,
        pam_patterns: Sequence[str] = ("NGG",),
        guide_length: int = 20,
        max_mismatches: int = 1,
    ) -> Dict[str, int | bool]:
        """
        Count genome-wide guide matches with PAM filtering.

        Returns:
        {
            "genome_match_count": exact PAM-aware matches,
            "off_target_count": exact matches beyond the first,
            "exact_pam_match_count": exact PAM-aware matches,
            "one_mismatch_pam_match_count": 1-mismatch PAM-aware matches,
            "is_unique": exact_pam_match_count == 1,
        }
        """
        guide = guide.upper()
        cache_key = (guide, tuple(p.upper() for p in pam_patterns))
        if cache_key in self._off_target_cache:
            return dict(self._off_target_cache[cache_key])

        if len(guide) != guide_length:
            raise ValueError(f"Expected guide length {guide_length}, got {len(guide)} for {guide}")

        if not pam_patterns:
            raise ValueError("pam_patterns cannot be empty")

        exact_count = 0
        one_mismatch_count = 0
        pam_len = len(pam_patterns[0])

        for chrom_seq in self.chrom_sequences.values():
            seq_len = len(chrom_seq)
            if seq_len < guide_length + pam_len:
                continue

            max_i = seq_len - guide_length - pam_len + 1

            # plus-strand targets
            for i in range(max_i):
                protospacer = chrom_seq[i:i + guide_length]
                pam = chrom_seq[i + guide_length:i + guide_length + pam_len]

                if "N" in protospacer or "N" in pam:
                    continue
                if not self._pam_matches(pam, pam_patterns):
                    continue

                mm = self._mismatch_count(guide, protospacer, max_mismatches=max_mismatches)
                if mm == 0:
                    exact_count += 1
                elif mm == 1 and max_mismatches >= 1:
                    one_mismatch_count += 1

            # minus-strand targets
            rc_seq = reverse_complement(chrom_seq)
            for i in range(max_i):
                protospacer = rc_seq[i:i + guide_length]
                pam = rc_seq[i + guide_length:i + guide_length + pam_len]

                if "N" in protospacer or "N" in pam:
                    continue
                if not self._pam_matches(pam, pam_patterns):
                    continue

                mm = self._mismatch_count(guide, protospacer, max_mismatches=max_mismatches)
                if mm == 0:
                    exact_count += 1
                elif mm == 1 and max_mismatches >= 1:
                    one_mismatch_count += 1

        result: Dict[str, int | bool] = {
            "genome_match_count": exact_count,
            "off_target_count": max(0, exact_count - 1),
            "exact_pam_match_count": exact_count,
            "one_mismatch_pam_match_count": one_mismatch_count,
            "is_unique": exact_count == 1,
        }
        self._off_target_cache[cache_key] = dict(result)
        return result

    def preload_guides_offtargets(
        self,
        guides: Iterable[str],
        pam_patterns: Sequence[str] = ("NGG",),
        guide_length: int = 20,
        max_mismatches: int = 1,
    ) -> None:
        unique_guides = sorted({g.upper() for g in guides if len(g) == guide_length})
        for guide in unique_guides:
            self.count_guide_offtargets(
                guide=guide,
                pam_patterns=pam_patterns,
                guide_length=guide_length,
                max_mismatches=max_mismatches,
            )