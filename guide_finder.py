from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


class GuideFinder:
    """
    Finds SpCas9 guides in a region and annotates them against a provided annotation store.

    Required annotation store interface:
      - get_features(chrom, start, end) -> list[dict]
    """

    def __init__(
        self,
        genome_store,
        annotation_store,
        guide_length: int = 20,
        pam_patterns: Sequence[str] = ("NGG",),
        relaxed_pam_patterns: Sequence[str] = ("NGG",),
    ) -> None:
        self.genome_store = genome_store
        self.annotation_store = annotation_store
        self.guide_length = guide_length
        self.pam_patterns = tuple(p.upper() for p in pam_patterns)
        self.relaxed_pam_patterns = tuple(p.upper() for p in relaxed_pam_patterns)

    @staticmethod
    def _pam_matches(pam_seq: str, pam_patterns: Sequence[str]) -> bool:
        pam_seq = pam_seq.upper()
        for pattern in pam_patterns:
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
    def _intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
        return not (a_end < b_start or b_end < a_start)

    @staticmethod
    def _sequence_distance(a: str, b: str) -> int:
        return sum(1 for x, y in zip(a, b) if x != y)

    @staticmethod
    def _gc_percent(seq: str) -> float:
        if not seq:
            return 0.0
        gc = sum(1 for b in seq if b in {"G", "C"})
        return round(100.0 * gc / len(seq), 1)

    @staticmethod
    def _longest_homopolymer(seq: str) -> int:
        best = 0
        current = 0
        prev = None
        for ch in seq:
            if ch == prev:
                current += 1
            else:
                current = 1
                prev = ch
            if current > best:
                best = current
        return best

    def _extract_relevant_features(self, chrom: str, start: int, end: int) -> List[Dict]:
        return self.annotation_store.get_features(chrom, start, end)

    def _collect_overlaps(
        self,
        chrom: str,
        guide_start: int,
        guide_end: int,
        features: Iterable[Dict],
    ) -> Dict:
        overlapping_genes = set()
        overlapping_transcripts = set()
        overlapping_exons: List[Dict] = []
        overlapping_cds: List[Dict] = []

        for feat in features:
            feat_start = int(feat["start"])
            feat_end = int(feat["end"])
            if not self._intervals_overlap(guide_start, guide_end, feat_start, feat_end):
                continue

            gene_name = feat.get("gene_name")
            transcript_id = feat.get("transcript_id")
            feature_type = str(feat.get("feature_type", "")).lower()

            if gene_name:
                overlapping_genes.add(gene_name)
            if transcript_id:
                overlapping_transcripts.add(transcript_id)

            if feature_type == "exon":
                overlapping_exons.append(
                    {
                        "gene_name": gene_name,
                        "transcript_id": transcript_id,
                        "exon_number": feat.get("exon_number"),
                        "start": feat_start,
                        "end": feat_end,
                        "strand": feat.get("strand"),
                    }
                )
            elif feature_type == "cds":
                overlapping_cds.append(
                    {
                        "gene_name": gene_name,
                        "transcript_id": transcript_id,
                        "exon_number": feat.get("exon_number"),
                        "start": feat_start,
                        "end": feat_end,
                        "strand": feat.get("strand"),
                    }
                )

        overlapping_exons.sort(key=lambda x: (x["start"], x["end"], str(x.get("transcript_id"))))
        overlapping_cds.sort(key=lambda x: (x["start"], x["end"], str(x.get("transcript_id"))))

        return {
            "overlapping_genes": sorted(overlapping_genes),
            "overlapping_transcripts": sorted(overlapping_transcripts),
            "overlapping_exons": overlapping_exons,
            "overlapping_cds": overlapping_cds,
            "overlapping_exon_count": len(overlapping_exons),
            "overlapping_cds_count": len(overlapping_cds),
            "overlaps_exon": len(overlapping_exons) > 0,
            "overlaps_cds": len(overlapping_cds) > 0,
        }

    def _best_exon_number_for_gene(self, overlaps: List[Dict], target_gene: Optional[str]) -> Optional[int]:
        if not target_gene:
            return None
        values: List[int] = []
        for entry in overlaps:
            if entry.get("gene_name") != target_gene:
                continue
            exon_number = entry.get("exon_number")
            try:
                values.append(int(exon_number))
            except (TypeError, ValueError):
                continue
        return min(values) if values else None

    def _annotate_target_specific_fields(self, guide: Dict, target_gene: Optional[str]) -> None:
        overlapping_genes = guide["overlapping_genes"]
        overlapping_exons = guide["overlapping_exons"]
        overlapping_cds = guide["overlapping_cds"]

        guide["target_gene"] = target_gene
        guide["hits_target_gene"] = bool(target_gene and target_gene in overlapping_genes)
        guide["target_gene_only"] = bool(
            target_gene and guide["hits_target_gene"] and set(overlapping_genes) == {target_gene}
        )
        guide["non_target_gene_count"] = 0
        if target_gene:
            guide["non_target_gene_count"] = sum(1 for g in overlapping_genes if g != target_gene)

        target_transcripts = {
            x["transcript_id"]
            for x in overlapping_exons
            if x.get("gene_name") == target_gene and x.get("transcript_id")
        }
        target_cds_transcripts = {
            x["transcript_id"]
            for x in overlapping_cds
            if x.get("gene_name") == target_gene and x.get("transcript_id")
        }

        guide["target_transcript_count"] = len(target_transcripts)
        guide["target_exon_count"] = sum(1 for x in overlapping_exons if x.get("gene_name") == target_gene)
        guide["target_cds_count"] = sum(1 for x in overlapping_cds if x.get("gene_name") == target_gene)
        guide["hits_target_cds"] = bool(target_gene and guide["target_cds_count"] > 0)

        guide["best_target_exon_number"] = self._best_exon_number_for_gene(overlapping_exons, target_gene)
        guide["best_target_cds_exon_number"] = self._best_exon_number_for_gene(overlapping_cds, target_gene)

    def _score_guide(self, guide: Dict) -> int:
        score = 0

        if guide["overlaps_exon"]:
            score += 25
        if guide["overlaps_cds"]:
            score += 50
        if guide["hits_target_gene"]:
            score += 25
        if guide["target_gene_only"]:
            score += 20
        if guide["hits_target_cds"]:
            score += 25

        score += min(guide["target_transcript_count"], 30)
        score += min(guide["target_cds_count"], 30)

        gc = guide["gc_percent"]
        if 40.0 <= gc <= 60.0:
            score += 20
        elif 35.0 <= gc <= 65.0:
            score += 10

        if guide["longest_homopolymer"] <= 2:
            score += 10
        elif guide["longest_homopolymer"] == 3:
            score += 5

        if not guide["has_tttt"]:
            score += 5

        exact_count = int(guide.get("exact_pam_match_count", guide.get("genome_match_count", 0)))
        one_mm_count = int(guide.get("one_mismatch_pam_match_count", 0))

        if exact_count == 1:
            score += 25
        else:
            score -= (exact_count - 1) * 40

        score -= one_mm_count * 10

        return int(score)

    def _group_key(self, guide: Dict) -> tuple:
        """
        Group near-duplicate genomic hits into a local bucket.

        Using exon/CDS-aware coarse locality:
        - chromosome position bucket
        - strand
        - target gene
        """
        bucket_size = 50
        bucket = guide["guide_genomic_start"] // bucket_size
        return (
            bucket,
            guide["strand"],
            guide.get("target_gene"),
        )

    def _deduplicate_by_sequence_best(self, guides: List[Dict]) -> List[Dict]:
        """
        Keep best-scoring site for each exact guide sequence.
        """
        best_by_seq: Dict[str, Dict] = {}
        for guide in guides:
            seq = guide["guide"]
            prev = best_by_seq.get(seq)
            if prev is None:
                best_by_seq[seq] = guide
                continue

            prev_key = (
                prev["score"],
                int(prev.get("hits_target_cds", False)),
                int(prev.get("target_gene_only", False)),
                prev.get("target_cds_count", 0),
                -prev.get("one_mismatch_pam_match_count", 0),
                -prev.get("off_target_count", 0),
            )
            curr_key = (
                guide["score"],
                int(guide.get("hits_target_cds", False)),
                int(guide.get("target_gene_only", False)),
                guide.get("target_cds_count", 0),
                -guide.get("one_mismatch_pam_match_count", 0),
                -guide.get("off_target_count", 0),
            )
            if curr_key > prev_key:
                best_by_seq[seq] = guide

        return list(best_by_seq.values())

    def _apply_diversity_filter(
        self,
        guides: List[Dict],
        max_guides: Optional[int] = None,
        min_position_spacing: int = 25,
        max_sequence_distance_for_redundancy: int = 2,
        one_per_local_bucket: bool = True,
    ) -> List[Dict]:
        """
        Build a non-redundant diverse list from an already-ranked guide list.

        Rules:
        - suppress guides too close on same strand
        - suppress guides differing by only a couple bases nearby
        - optionally allow only one guide per local genomic bucket
        """
        selected: List[Dict] = []
        used_buckets = set()

        for guide in guides:
            if one_per_local_bucket:
                gkey = self._group_key(guide)
                if gkey in used_buckets:
                    continue

            redundant = False
            for kept in selected:
                if guide["strand"] != kept["strand"]:
                    continue

                pos_close = abs(guide["guide_genomic_start"] - kept["guide_genomic_start"]) < min_position_spacing
                seq_close = self._sequence_distance(guide["guide"], kept["guide"]) <= max_sequence_distance_for_redundancy

                # same target neighborhood and very similar sequence
                if pos_close and seq_close:
                    redundant = True
                    break

                # also suppress exact-sequence duplicate if one slipped through
                if guide["guide"] == kept["guide"]:
                    redundant = True
                    break

            if redundant:
                continue

            selected.append(guide)
            if one_per_local_bucket:
                used_buckets.add(self._group_key(guide))

            if max_guides is not None and len(selected) >= max_guides:
                break

        return selected

    def find_guides_in_region(
        self,
        chrom: str,
        region_start: int,
        region_end: int,
        target_gene: Optional[str] = None,
        require_target_gene: bool = False,
        require_exon_overlap: bool = True,
        require_cds_overlap: bool = False,
        max_one_mismatch_offtargets: Optional[int] = None,
        require_unique_exact: bool = False,
        deduplicate_exact_sequences: bool = False,
    ) -> List[Dict]:
        seq = self.genome_store.get_sequence(chrom, region_start, region_end).upper()
        features = self._extract_relevant_features(chrom, region_start, region_end)

        pam_len = len(self.pam_patterns[0])
        guides: List[Dict] = []

        max_i = len(seq) - self.guide_length - pam_len + 1

        # plus-strand guides
        for i in range(max_i):
            protospacer = seq[i:i + self.guide_length]
            pam = seq[i + self.guide_length:i + self.guide_length + pam_len]
            if "N" in protospacer or "N" in pam:
                continue
            if not self._pam_matches(pam, self.pam_patterns):
                continue

            guide_start = region_start + i
            guide_end = guide_start + self.guide_length - 1

            guide = {
                "guide": protospacer,
                "pam": pam,
                "pam_genomic_start": guide_end + 1,
                "guide_genomic_start": guide_start,
                "guide_genomic_end": guide_end,
                "strand": "+",
            }
            guide.update(self._collect_overlaps(chrom, guide_start, guide_end, features))
            guide["gc_percent"] = self._gc_percent(protospacer)
            guide["longest_homopolymer"] = self._longest_homopolymer(protospacer)
            guide["has_tttt"] = "TTTT" in protospacer
            self._annotate_target_specific_fields(guide, target_gene)
            guides.append(guide)

        # minus-strand guides
        rc_seq = reverse_complement(seq)
        for i in range(max_i):
            protospacer_rc = rc_seq[i:i + self.guide_length]
            pam_rc = rc_seq[i + self.guide_length:i + self.guide_length + pam_len]
            if "N" in protospacer_rc or "N" in pam_rc:
                continue
            if not self._pam_matches(pam_rc, self.pam_patterns):
                continue

            guide_seq = protospacer_rc
            orig_window_end = region_end - i
            guide_end = orig_window_end - pam_len
            guide_start = guide_end - self.guide_length + 1
            pam_start = guide_end + 1

            guide = {
                "guide": guide_seq,
                "pam": pam_rc,
                "pam_genomic_start": pam_start,
                "guide_genomic_start": guide_start,
                "guide_genomic_end": guide_end,
                "strand": "-",
            }
            guide.update(self._collect_overlaps(chrom, guide_start, guide_end, features))
            guide["gc_percent"] = self._gc_percent(guide_seq)
            guide["longest_homopolymer"] = self._longest_homopolymer(guide_seq)
            guide["has_tttt"] = "TTTT" in guide_seq
            self._annotate_target_specific_fields(guide, target_gene)
            guides.append(guide)

        unique_guides = {g["guide"] for g in guides}
        self.genome_store.preload_guides_offtargets(
            unique_guides,
            pam_patterns=self.relaxed_pam_patterns,
            guide_length=self.guide_length,
            max_mismatches=1,
        )

        for guide in guides:
            off = self.genome_store.count_guide_offtargets(
                guide["guide"],
                pam_patterns=self.relaxed_pam_patterns,
                guide_length=self.guide_length,
                max_mismatches=1,
            )
            guide.update(off)
            guide["score"] = self._score_guide(guide)

        filtered: List[Dict] = []
        for guide in guides:
            if require_target_gene and not guide["hits_target_gene"]:
                continue
            if require_exon_overlap and not guide["overlaps_exon"]:
                continue
            if require_cds_overlap and not guide["overlaps_cds"]:
                continue
            if require_unique_exact and not guide.get("is_unique", False):
                continue
            if max_one_mismatch_offtargets is not None:
                if int(guide.get("one_mismatch_pam_match_count", 0)) > max_one_mismatch_offtargets:
                    continue
            filtered.append(guide)

        filtered.sort(
            key=lambda g: (
                -g["score"],
                -int(g.get("hits_target_cds", False)),
                -int(g.get("target_gene_only", False)),
                -int(g.get("target_cds_count", 0)),
                -int(g.get("target_transcript_count", 0)),
                int(g.get("one_mismatch_pam_match_count", 0)),
                int(g.get("off_target_count", 0)),
                abs(g["gc_percent"] - 50.0),
                g["guide_genomic_start"],
                g["guide"],
            )
        )

        if deduplicate_exact_sequences:
            filtered = self._deduplicate_by_sequence_best(filtered)
            filtered.sort(
                key=lambda g: (
                    -g["score"],
                    -int(g.get("hits_target_cds", False)),
                    -int(g.get("target_gene_only", False)),
                    -int(g.get("target_cds_count", 0)),
                    -int(g.get("target_transcript_count", 0)),
                    int(g.get("one_mismatch_pam_match_count", 0)),
                    int(g.get("off_target_count", 0)),
                    abs(g["gc_percent"] - 50.0),
                    g["guide_genomic_start"],
                    g["guide"],
                )
            )

        return filtered

    def select_diverse_top_guides(
        self,
        ranked_guides: List[Dict],
        top_n: int = 10,
        min_position_spacing: int = 25,
        max_sequence_distance_for_redundancy: int = 2,
        one_per_local_bucket: bool = True,
    ) -> List[Dict]:
        return self._apply_diversity_filter(
            guides=ranked_guides,
            max_guides=top_n,
            min_position_spacing=min_position_spacing,
            max_sequence_distance_for_redundancy=max_sequence_distance_for_redundancy,
            one_per_local_bucket=one_per_local_bucket,
        )