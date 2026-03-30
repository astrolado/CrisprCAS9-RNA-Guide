import sys

from core.genome_store import GenomeStore
from core.annotation_store import AnnotationStore
from core.guide_finder import (
    find_ngg_guides,
    annotate_guides,
    add_simple_off_target_counts,
    add_target_gene_ranking,
    filter_guides,
    rank_guides,
)


def parse_region(region_string):
    chrom, coords = region_string.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def main():
    if len(sys.argv) not in {2, 3}:
        print("Usage:")
        print("  python -m cli.main chr:start-end")
        print("  python -m cli.main chr:start-end TARGET_GENE")
        return

    region_string = sys.argv[1]
    target_gene = sys.argv[2] if len(sys.argv) == 3 else None

    chrom, start, end = parse_region(region_string)

    print("Loading genome...")
    genome = GenomeStore("data/genome.fa")

    print("Loading annotations...")
    annot = AnnotationStore("data/annotation.gtf")

    print("Querying region...")
    results = annot.query_region(chrom, start, end)
    summary = annot.summarize_region(chrom, start, end)

    print(f"Region: {chrom}:{start}-{end}")
    print("Summary:", summary)

    if target_gene:
        print("Target gene:", target_gene)

    print("Getting sequence...")
    seq = genome.get_sequence(chrom, start, end)

    print("Sequence length:", len(seq))
    print("Sequence sample:", seq[:50])

    guides = find_ngg_guides(seq, chrom, start)
    guides = annotate_guides(guides, results)
    guides = add_simple_off_target_counts(guides, genome)
    guides = add_target_gene_ranking(guides, target_gene=target_gene)

    raw_count = len(guides)
    plus_count = sum(1 for g in guides if g["strand"] == "+")
    minus_count = sum(1 for g in guides if g["strand"] == "-")
    exon_count = sum(1 for g in guides if g["overlaps_exon"])
    cds_count = sum(1 for g in guides if g["overlaps_cds"])
    unique_count = sum(1 for g in guides if g["is_unique"])

    if target_gene:
        target_gene_hit_count = sum(1 for g in guides if g["hits_target_gene"])
        target_gene_only_count = sum(1 for g in guides if g["target_gene_only"])
        target_cds_hit_count = sum(1 for g in guides if g["hits_target_cds"])
        print("Guides hitting target gene:", target_gene_hit_count)
        print("Guides hitting only target gene:", target_gene_only_count)
        print("Guides hitting target CDS:", target_cds_hit_count)

    print("Unique exact-match guides:", unique_count)

    filtered_guides = filter_guides(guides, target_gene=target_gene)
    ranked_guides = rank_guides(filtered_guides, target_gene=target_gene)

    print("Number of raw guides found:", raw_count)
    print("Plus strand guides:", plus_count)
    print("Minus strand guides:", minus_count)
    print("Exon-overlapping guides:", exon_count)
    print("CDS-overlapping guides:", cds_count)
    print("Filtered guides:", len(filtered_guides))

    print("Top 10 ranked guides:")
    for g in ranked_guides[:10]:
        print(g)

    print("Detailed view of the top 3 guides:")
    for g in ranked_guides[:3]:
        print("-" * 60)
        print("Guide:", g["guide"])
        print("Strand:", g["strand"])
        print("PAM:", g["pam"])
        print("Guide coordinates:", f'{g["guide_genomic_start"]}-{g["guide_genomic_end"]}')
        print("Score:", g["score"])
        print("GC%:", g["gc_percent"])
        print("Genome match count:", g["genome_match_count"])
        print("Off-target count:", g["off_target_count"])
        print("Is unique:", g["is_unique"])
        print("Overlapping genes:", g["overlapping_genes"])
        print("Hits target gene:", g["hits_target_gene"])
        print("Target gene only:", g["target_gene_only"])
        print("Non-target gene count:", g["non_target_gene_count"])
        print("Hits target CDS:", g["hits_target_cds"])
        print("Target transcript count:", g["target_transcript_count"])
        print("Target CDS count:", g["target_cds_count"])
        print("Best target exon number:", g["best_target_exon_number"])
        print("Best target CDS exon number:", g["best_target_cds_exon_number"])
        print("Overlapping transcripts:", g["overlapping_transcripts"])
        print("Overlapping exon count:", g["overlapping_exon_count"])
        print("Overlapping CDS count:", g["overlapping_cds_count"])
        print("Overlapping exons:", g["overlapping_exons"])
        print("Overlapping CDS:", g["overlapping_cds"])


if __name__ == "__main__":
    main()