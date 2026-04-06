from __future__ import annotations

from typing import Iterable, Sequence
from tqdm import tqdm


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


class GenomeStore:
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.sequences: dict[str, str] = {}
        self.loaded = False

        self._exact_match_cache: dict[str, int] = {}
        self._one_mismatch_cache: dict[str, int] = {}
        self._pam_offtarget_cache: dict[tuple[str, tuple[str, ...], int], dict] = {}

        self.load()

    def load(self) -> None:
        if self.loaded:
            return
        self.load_fasta(self.fasta_path)
        self.loaded = True

    def load_fasta(self, path: str) -> None:
        current_chrom = None
        chunks = []

        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if current_chrom is not None:
                        self.sequences[current_chrom] = "".join(chunks).upper()

                    header = line[1:].split()[0]
                    current_chrom = header
                    chunks = []
                else:
                    chunks.append(line)

        if current_chrom is not None:
            self.sequences[current_chrom] = "".join(chunks).upper()

    def chromosome_names(self) -> list[str]:
        return list(self.sequences.keys())

    def has_chromosome(self, chrom: str) -> bool:
        try:
            self._resolve_chrom_name(chrom)
            return True
        except KeyError:
            return False

    def _resolve_chrom_name(self, chrom: str) -> str:
        if chrom in self.sequences:
            return chrom

        candidates = []
        if chrom.startswith("chr"):
            candidates.append(chrom[3:])
        else:
            candidates.append("chr" + chrom)

        for candidate in candidates:
            if candidate in self.sequences:
                return candidate

        raise KeyError(
            f"Chromosome not found: {chrom}. "
            f"Available examples: {list(self.sequences.keys())[:10]}"
        )

    def get_sequence(self, chrom: str, start: int, end: int) -> str:
        resolved_chrom = self._resolve_chrom_name(chrom)
        seq = self.sequences[resolved_chrom]

        if start < 1 or end > len(seq) or start > end:
            raise ValueError(
                f"Invalid coordinates for {resolved_chrom}: {start}-{end} "
                f"(length {len(seq)})"
            )

        return seq[start - 1:end]

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
    def _count_mismatches(a: str, b: str, max_mismatches: int | None = None) -> int:
        mismatches = 0
        for x, y in zip(a, b):
            if x != y:
                mismatches += 1
                if max_mismatches is not None and mismatches > max_mismatches:
                    return mismatches
        return mismatches

    def _iter_pam_hits_on_plus_strand(
        self,
        chrom_name: str,
        chrom_seq: str,
        guide: str,
        pam_patterns: Sequence[str],
        max_mismatches: int,
    ):
        guide_len = len(guide)
        pam_len = len(pam_patterns[0])
        max_i = len(chrom_seq) - guide_len - pam_len + 1

        for i in range(max_i):
            protospacer = chrom_seq[i:i + guide_len]
            pam = chrom_seq[i + guide_len:i + guide_len + pam_len]

            if "N" in protospacer or "N" in pam:
                continue
            if not self._pam_matches(pam, pam_patterns):
                continue

            mismatches = self._count_mismatches(protospacer, guide, max_mismatches)
            if mismatches > max_mismatches:
                continue

            guide_start = i + 1
            guide_end = guide_start + guide_len - 1
            pam_start = guide_end + 1
            pam_end = pam_start + pam_len - 1

            yield {
                "chrom": chrom_name,
                "strand": "+",
                "guide_start": guide_start,
                "guide_end": guide_end,
                "pam_start": pam_start,
                "pam_end": pam_end,
                "protospacer": protospacer,
                "pam": pam,
                "mismatches": mismatches,
            }

    def _iter_pam_hits_on_minus_strand(
        self,
        chrom_name: str,
        chrom_seq: str,
        guide: str,
        pam_patterns: Sequence[str],
        max_mismatches: int,
    ):
        rc_seq = reverse_complement(chrom_seq)
        guide_len = len(guide)
        pam_len = len(pam_patterns[0])
        chrom_len = len(chrom_seq)
        max_i = len(rc_seq) - guide_len - pam_len + 1

        for i in range(max_i):
            protospacer = rc_seq[i:i + guide_len]
            pam = rc_seq[i + guide_len:i + guide_len + pam_len]

            if "N" in protospacer or "N" in pam:
                continue
            if not self._pam_matches(pam, pam_patterns):
                continue

            mismatches = self._count_mismatches(protospacer, guide, max_mismatches)
            if mismatches > max_mismatches:
                continue

            orig_window_end = chrom_len - i
            guide_end = orig_window_end - pam_len
            guide_start = guide_end - guide_len + 1
            pam_start = guide_end + 1
            pam_end = pam_start + pam_len - 1

            yield {
                "chrom": chrom_name,
                "strand": "-",
                "guide_start": guide_start,
                "guide_end": guide_end,
                "pam_start": pam_start,
                "pam_end": pam_end,
                "protospacer": protospacer,
                "pam": pam,
                "mismatches": mismatches,
            }

    def _scan_one_strand_with_progress(
        self,
        chrom_name: str,
        chrom_seq: str,
        guide: str,
        pam_patterns: Sequence[str],
        max_mismatches: int,
        strand: str,
        pbar,
    ):
        guide_len = len(guide)
        pam_len = len(pam_patterns[0])

        if strand == "+":
            seq_to_scan = chrom_seq
            chrom_len = len(chrom_seq)
        else:
            seq_to_scan = reverse_complement(chrom_seq)
            chrom_len = len(chrom_seq)

        max_i = len(seq_to_scan) - guide_len - pam_len + 1
        step = 500000

        for block_start in range(0, max_i, step):
            block_end = min(block_start + step, max_i)

            for i in range(block_start, block_end):
                protospacer = seq_to_scan[i:i + guide_len]
                pam = seq_to_scan[i + guide_len:i + guide_len + pam_len]

                if "N" in protospacer or "N" in pam:
                    continue
                if not self._pam_matches(pam, pam_patterns):
                    continue

                mismatches = self._count_mismatches(protospacer, guide, max_mismatches)
                if mismatches > max_mismatches:
                    continue

                if strand == "+":
                    guide_start = i + 1
                    guide_end = guide_start + guide_len - 1
                    pam_start = guide_end + 1
                    pam_end = pam_start + pam_len - 1
                else:
                    orig_window_end = chrom_len - i
                    guide_end = orig_window_end - pam_len
                    guide_start = guide_end - guide_len + 1
                    pam_start = guide_end + 1
                    pam_end = pam_start + pam_len - 1

                yield {
                    "chrom": chrom_name,
                    "strand": strand,
                    "guide_start": guide_start,
                    "guide_end": guide_end,
                    "pam_start": pam_start,
                    "pam_end": pam_end,
                    "protospacer": protospacer,
                    "pam": pam,
                    "mismatches": mismatches,
                }

            pbar.update(block_end - block_start)

    def _scan_pam_aware_offtargets(
        self,
        guide: str,
        pam_patterns: Sequence[str],
        max_mismatches: int,
    ) -> dict:
        guide = guide.upper()
        pam_patterns = tuple(p.upper() for p in pam_patterns)

        exact_hits = []
        mismatch_hits = []

        guide_len = len(guide)
        pam_len = len(pam_patterns[0])

        total_windows = 0
        for chrom_seq in self.sequences.values():
            windows = len(chrom_seq) - guide_len - pam_len + 1
            if windows > 0:
                total_windows += windows * 2

        with tqdm(
            total=total_windows,
            desc=f"[offtarget] scanning {guide}",
            unit="win",
            leave=False,
        ) as pbar:
            for chrom_name, chrom_seq in self.sequences.items():
                for hit in self._scan_one_strand_with_progress(
                    chrom_name=chrom_name,
                    chrom_seq=chrom_seq,
                    guide=guide,
                    pam_patterns=pam_patterns,
                    max_mismatches=max_mismatches,
                    strand="+",
                    pbar=pbar,
                ):
                    if hit["mismatches"] == 0:
                        exact_hits.append(hit)
                    else:
                        mismatch_hits.append(hit)

                for hit in self._scan_one_strand_with_progress(
                    chrom_name=chrom_name,
                    chrom_seq=chrom_seq,
                    guide=guide,
                    pam_patterns=pam_patterns,
                    max_mismatches=max_mismatches,
                    strand="-",
                    pbar=pbar,
                ):
                    if hit["mismatches"] == 0:
                        exact_hits.append(hit)
                    else:
                        mismatch_hits.append(hit)

        exact_count = len(exact_hits)
        one_mismatch_total = exact_count + len(mismatch_hits)

        return {
            "guide": guide,
            "pam_patterns": list(pam_patterns),
            "max_mismatches": max_mismatches,
            "exact_pam_match_count": exact_count,
            "genome_match_count": exact_count,
            "off_target_count": max(0, exact_count - 1),
            "is_unique": exact_count == 1,
            "one_mismatch_match_count": one_mismatch_total,
            "one_mismatch_offtarget_count": max(0, one_mismatch_total - 1),
            "one_mismatch_pam_match_count": max(0, one_mismatch_total - 1),
            "exact_hits_preview": exact_hits[:10],
            "mismatch_hits_preview": mismatch_hits[:10],
        }

    def count_guide_offtargets(
        self,
        guide: str,
        pam_patterns: Sequence[str] | None = None,
        guide_length: int | None = None,
        max_mismatches: int = 1,
    ) -> dict:
        guide = guide.upper()

        if guide_length is not None and len(guide) != guide_length:
            raise ValueError(
                f"Guide length mismatch: expected {guide_length}, got {len(guide)}"
            )

        if pam_patterns is None:
            pam_patterns = ("NGG",)

        pam_patterns = tuple(p.upper() for p in pam_patterns)
        if not pam_patterns:
            raise ValueError("pam_patterns must contain at least one PAM pattern")

        pam_len_set = {len(p) for p in pam_patterns}
        if len(pam_len_set) != 1:
            raise ValueError("All PAM patterns must have the same length")

        cache_key = (guide, pam_patterns, max_mismatches)
        cached = self._pam_offtarget_cache.get(cache_key)
        if cached is not None:
            return dict(cached)

        result = self._scan_pam_aware_offtargets(
            guide=guide,
            pam_patterns=pam_patterns,
            max_mismatches=max_mismatches,
        )

        self._pam_offtarget_cache[cache_key] = dict(result)
        return dict(result)

    def preload_guides_offtargets(
        self,
        guides: Iterable,
        count_1mm: bool = False,
        include_1mm: bool = False,
        verbose: bool = False,
        *args,
        **kwargs,
    ) -> None:
        unique_guides = []
        seen = set()

        for item in guides:
            if isinstance(item, dict):
                guide_seq = item.get("guide")
            else:
                guide_seq = str(item)

            if not guide_seq:
                continue

            guide_seq = guide_seq.upper()
            if guide_seq not in seen:
                seen.add(guide_seq)
                unique_guides.append(guide_seq)

        pam_patterns = kwargs.get("pam_patterns", ("NGG",))
        guide_length = kwargs.get("guide_length", None)
        max_mismatches = kwargs.get("max_mismatches", 1)

        do_1mm = count_1mm or include_1mm or (max_mismatches >= 1)
        mm = 1 if do_1mm else 0

        guide_iter = unique_guides
        if verbose:
            guide_iter = tqdm(unique_guides, desc="[offtarget] guides", unit="guide")

        for guide in guide_iter:
            self.count_guide_offtargets(
                guide=guide,
                pam_patterns=pam_patterns,
                guide_length=guide_length,
                max_mismatches=mm,
            )

    def annotate_guide_offtargets(
        self,
        guide_record: dict,
        pam_patterns: Sequence[str] | None = None,
        guide_length: int | None = None,
        max_mismatches: int = 1,
    ) -> dict:
        off = self.count_guide_offtargets(
            guide=guide_record["guide"],
            pam_patterns=pam_patterns or ("NGG",),
            guide_length=guide_length,
            max_mismatches=max_mismatches,
        )
        guide_record.update(off)
        return guide_record