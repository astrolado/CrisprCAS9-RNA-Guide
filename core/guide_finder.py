def reverse_complement(seq):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    return "".join(complement[base] for base in reversed(seq))


def ranges_overlap(start1, end1, start2, end2):
    return start1 <= end2 and end1 >= start2


def gc_content(seq):
    gc_count = seq.count("G") + seq.count("C")
    return round((gc_count / len(seq)) * 100, 2)


def longest_homopolymer_run(seq):
    longest = 1
    current = 1

    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            if current > longest:
                longest = current
        else:
            current = 1

    return longest


def find_ngg_guides(seq, chrom, region_start):
    candidates = []

    for i in range(len(seq) - 2):
        pam = seq[i:i + 3]

        if pam[1:] == "GG" and i >= 20:
            guide = seq[i - 20:i]

            candidates.append({
                "guide": guide,
                "pam": pam,
                "pam_genomic_start": region_start + i,
                "guide_genomic_start": region_start + i - 20,
                "guide_genomic_end": region_start + i - 1,
                "strand": "+"
            })

    for i in range(len(seq) - 2):
        pam_forward = seq[i:i + 3]

        if pam_forward[:2] == "CC" and i + 3 + 20 <= len(seq):
            target_seq = seq[i + 3:i + 23]
            guide = reverse_complement(target_seq)
            pam = reverse_complement(pam_forward)

            candidates.append({
                "guide": guide,
                "pam": pam,
                "pam_genomic_start": region_start + i,
                "guide_genomic_start": region_start + i + 3,
                "guide_genomic_end": region_start + i + 22,
                "strand": "-"
            })

    return candidates


def annotate_guides(guides, region_features):
    for guide in guides:
        overlapping_genes = set()
        overlapping_transcripts = set()
        overlapping_exons = []
        overlapping_cds = []

        overlaps_exon = False
        overlaps_cds = False

        for feature in region_features:
            if ranges_overlap(
                guide["guide_genomic_start"],
                guide["guide_genomic_end"],
                feature["start"],
                feature["end"]
            ):
                if feature["gene_name"]:
                    overlapping_genes.add(feature["gene_name"])

                if feature["transcript_id"]:
                    overlapping_transcripts.add(feature["transcript_id"])

                if feature["feature"] == "exon":
                    overlaps_exon = True
                    overlapping_exons.append({
                        "gene_name": feature["gene_name"],
                        "transcript_id": feature["transcript_id"],
                        "exon_number": feature["exon_number"],
                        "start": feature["start"],
                        "end": feature["end"],
                        "strand": feature["strand"],
                    })

                if feature["feature"] == "CDS":
                    overlaps_cds = True
                    overlapping_cds.append({
                        "gene_name": feature["gene_name"],
                        "transcript_id": feature["transcript_id"],
                        "exon_number": feature["exon_number"],
                        "start": feature["start"],
                        "end": feature["end"],
                        "strand": feature["strand"],
                    })

        unique_exons = []
        seen_exons = set()
        for exon in overlapping_exons:
            key = (
                exon["gene_name"],
                exon["transcript_id"],
                exon["exon_number"],
                exon["start"],
                exon["end"]
            )
            if key not in seen_exons:
                seen_exons.add(key)
                unique_exons.append(exon)

        unique_cds = []
        seen_cds = set()
        for cds in overlapping_cds:
            key = (
                cds["gene_name"],
                cds["transcript_id"],
                cds["exon_number"],
                cds["start"],
                cds["end"]
            )
            if key not in seen_cds:
                seen_cds.add(key)
                unique_cds.append(cds)

        guide["overlaps_exon"] = overlaps_exon
        guide["overlaps_cds"] = overlaps_cds
        guide["overlapping_genes"] = sorted(overlapping_genes)
        guide["overlapping_transcripts"] = sorted(overlapping_transcripts)
        guide["overlapping_exons"] = unique_exons
        guide["overlapping_exon_count"] = len(unique_exons)
        guide["overlapping_cds"] = unique_cds
        guide["overlapping_cds_count"] = len(unique_cds)
        guide["gc_percent"] = gc_content(guide["guide"])
        guide["longest_homopolymer"] = longest_homopolymer_run(guide["guide"])
        guide["has_tttt"] = "TTTT" in guide["guide"]

    return guides


def add_simple_off_target_counts(guides, genome_store):
    guide_set = {g["guide"] for g in guides}
    counts = genome_store.count_selected_guides(guide_set)

    for guide in guides:
        match_count = counts.get(guide["guide"], 0)
        guide["genome_match_count"] = match_count
        guide["off_target_count"] = max(0, match_count - 1)
        guide["is_unique"] = match_count == 1

    return guides


def target_gene_metrics(guide, target_gene):
    if target_gene is None:
        guide["target_gene"] = None
        guide["hits_target_gene"] = False
        guide["target_gene_only"] = False
        guide["non_target_gene_count"] = max(0, len(guide["overlapping_genes"]))
        guide["target_transcript_count"] = 0
        guide["target_exon_count"] = 0
        guide["target_cds_count"] = 0
        guide["hits_target_cds"] = False
        guide["best_target_exon_number"] = None
        guide["best_target_cds_exon_number"] = None
        return guide

    target_exons = [
        exon for exon in guide["overlapping_exons"]
        if exon["gene_name"] == target_gene
    ]

    target_cds = [
        cds for cds in guide["overlapping_cds"]
        if cds["gene_name"] == target_gene
    ]

    target_transcripts = sorted({
        exon["transcript_id"]
        for exon in target_exons
        if exon["transcript_id"] is not None
    })

    exon_numbers = []
    for exon in target_exons:
        exon_number = exon["exon_number"]
        if exon_number is None:
            continue
        try:
            exon_numbers.append(int(exon_number))
        except ValueError:
            continue

    cds_exon_numbers = []
    for cds in target_cds:
        exon_number = cds["exon_number"]
        if exon_number is None:
            continue
        try:
            cds_exon_numbers.append(int(exon_number))
        except ValueError:
            continue

    hits_target_gene = target_gene in guide["overlapping_genes"]
    non_target_gene_count = len([
        gene for gene in guide["overlapping_genes"]
        if gene != target_gene
    ])

    guide["target_gene"] = target_gene
    guide["hits_target_gene"] = hits_target_gene
    guide["target_gene_only"] = hits_target_gene and non_target_gene_count == 0
    guide["non_target_gene_count"] = non_target_gene_count
    guide["target_transcript_count"] = len(target_transcripts)
    guide["target_exon_count"] = len(target_exons)
    guide["target_cds_count"] = len(target_cds)
    guide["hits_target_cds"] = len(target_cds) > 0
    guide["best_target_exon_number"] = min(exon_numbers) if exon_numbers else None
    guide["best_target_cds_exon_number"] = min(cds_exon_numbers) if cds_exon_numbers else None

    return guide


def basic_guide_score(guide, target_gene=None):
    score = 0
    gc = guide["gc_percent"]
    homopolymer = guide["longest_homopolymer"]

    if guide["overlaps_exon"]:
        score += 20

    if guide["overlaps_cds"]:
        score += 20

    if 40 <= gc <= 60:
        score += 20
    elif 30 <= gc <= 70:
        score += 10

    if guide["has_tttt"]:
        score -= 20

    if homopolymer >= 5:
        score -= 20
    elif homopolymer == 4:
        score -= 10

    if guide.get("is_unique") is True:
        score += 25

    score -= 15 * guide.get("off_target_count", 0)

    if target_gene is not None:
        if guide["hits_target_gene"]:
            score += 20

        if guide["target_gene_only"]:
            score += 20

        if guide["hits_target_cds"]:
            score += 25

        score += min(guide["target_transcript_count"], 20)
        score += min(guide["target_cds_count"], 20)

        if guide["best_target_cds_exon_number"] is not None:
            exon_num = guide["best_target_cds_exon_number"]
            if exon_num <= 3:
                score += 20
            elif exon_num <= 6:
                score += 12
            elif exon_num <= 10:
                score += 6
        elif guide["best_target_exon_number"] is not None:
            exon_num = guide["best_target_exon_number"]
            if exon_num <= 3:
                score += 8
            elif exon_num <= 6:
                score += 4

        score -= 15 * guide["non_target_gene_count"]

    return score


def add_target_gene_ranking(guides, target_gene=None):
    for guide in guides:
        target_gene_metrics(guide, target_gene)
        guide["score"] = basic_guide_score(guide, target_gene=target_gene)

    return guides


def filter_guides(guides, target_gene=None):
    filtered = []

    for guide in guides:
        if not guide["overlaps_exon"]:
            continue
        if guide["has_tttt"]:
            continue
        if guide["longest_homopolymer"] >= 5:
            continue
        if target_gene is not None and not guide["hits_target_gene"]:
            continue

        filtered.append(guide)

    return filtered


def rank_guides(guides, target_gene=None):
    return sorted(
        guides,
        key=lambda g: (
            -g["score"],
            g["off_target_count"],
            0 if g["hits_target_cds"] else 1,
            0 if g["target_gene_only"] else 1,
            -(g["target_cds_count"]),
            -(g["target_transcript_count"]),
            g["best_target_cds_exon_number"] if g["best_target_cds_exon_number"] is not None else 999,
            g["best_target_exon_number"] if g["best_target_exon_number"] is not None else 999,
            abs(g["gc_percent"] - 50),
            g["guide_genomic_start"]
        )
    )