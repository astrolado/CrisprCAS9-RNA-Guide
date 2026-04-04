from __future__ import annotations

from typing import Iterable


class GenomeStore:
    def __init__(self, fasta_path: str):
        self.fasta_path = fasta_path
        self.sequences: dict[str, str] = {}
        self.loaded = False

        # caches for simple off-target counting
        self._exact_match_cache: dict[str, int] = {}
        self._one_mismatch_cache: dict[str, int] = {}

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

    # ------------------------------------------------------------------
    # simple off-target counting
    # ------------------------------------------------------------------

    def _count_overlapping_exact(self, seq: str, guide: str) -> int:
        count = 0
        guide_len = len(guide)

        for i in range(len(seq) - guide_len + 1):
            if seq[i:i + guide_len] == guide:
                count += 1

        return count

    def _count_exact_matches_genomewide(self, guide: str) -> int:
        guide = guide.upper()
        if guide in self._exact_match_cache:
            return self._exact_match_cache[guide]

        total = 0
        for chrom_seq in self.sequences.values():
            total += self._count_overlapping_exact(chrom_seq, guide)

        self._exact_match_cache[guide] = total
        return total

    def _hamming_distance_leq_one(self, a: str, b: str) -> bool:
        mismatches = 0
        for x, y in zip(a, b):
            if x != y:
                mismatches += 1
                if mismatches > 1:
                    return False
        return True

    def _count_one_mismatch_or_less_genomewide(self, guide: str) -> int:
        guide = guide.upper()
        if guide in self._one_mismatch_cache:
            return self._one_mismatch_cache[guide]

        total = 0
        guide_len = len(guide)

        for chrom_seq in self.sequences.values():
            for i in range(len(chrom_seq) - guide_len + 1):
                window = chrom_seq[i:i + guide_len]
                if self._hamming_distance_leq_one(window, guide):
                    total += 1

        self._one_mismatch_cache[guide] = total
        return total

    def count_exact_matches(self, guide: str) -> int:
        return self._count_exact_matches_genomewide(guide)

    def count_1mm_matches(self, guide: str) -> int:
        return self._count_one_mismatch_or_less_genomewide(guide)

    def get_exact_match_count(self, guide: str) -> int:
        return self._count_exact_matches_genomewide(guide)

    def get_1mm_match_count(self, guide: str) -> int:
        return self._count_one_mismatch_or_less_genomewide(guide)

    def get_offtarget_count(self, guide: str) -> int:
        exact = self._count_exact_matches_genomewide(guide)
        return max(0, exact - 1)

    def get_1mm_offtarget_count(self, guide: str) -> int:
        total_1mm = self._count_one_mismatch_or_less_genomewide(guide)
        return max(0, total_1mm - 1)

    def is_unique_exact(self, guide: str) -> bool:
        return self._count_exact_matches_genomewide(guide) == 1

    def preload_guides_offtargets(
        self,
        guides: Iterable,
        count_1mm: bool = False,
        include_1mm: bool = False,
        verbose: bool = False,
        *args,
        **kwargs,
    ) -> None:
        """
        Precompute simple genome-wide exact and optional <=1 mismatch counts.

        Accepts either:
        - a list of guide strings
        - a list of guide dicts with a 'guide' key

        Extra args/kwargs are accepted so this stays compatible with whatever
        GuideFinder passes.
        """
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

        do_1mm = count_1mm or include_1mm or ("max_1mm_offtargets" in kwargs)

        if verbose:
            print(f"[offtarget] preloading off-target counts for {len(unique_guides)} unique guides")

        for idx, guide in enumerate(unique_guides, start=1):
            if verbose and (idx <= 5 or idx % 100 == 0):
                print(f"[offtarget] guide {idx}/{len(unique_guides)}: {guide}")

            self._count_exact_matches_genomewide(guide)

            if do_1mm:
                self._count_one_mismatch_or_less_genomewide(guide)

    def annotate_guide_offtargets(self, guide_record: dict) -> dict:
        """
        Convenience helper if GuideFinder wants a filled-in record.
        """
        guide = guide_record["guide"].upper()
        genome_match_count = self._count_exact_matches_genomewide(guide)
        off_target_count = max(0, genome_match_count - 1)

        guide_record["genome_match_count"] = genome_match_count
        guide_record["off_target_count"] = off_target_count
        guide_record["is_unique"] = genome_match_count == 1

        if guide in self._one_mismatch_cache:
            guide_record["one_mismatch_match_count"] = self._one_mismatch_cache[guide]
            guide_record["one_mismatch_offtarget_count"] = max(
                0, self._one_mismatch_cache[guide] - 1
            )

        return guide_record
