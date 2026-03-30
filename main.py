from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from cli.annotation_store import AnnotationStore
    from cli.genome_store import GenomeStore
    from cli.guide_finder import GuideFinder
except ModuleNotFoundError:
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from cli.annotation_store import AnnotationStore
    from cli.genome_store import GenomeStore
    from cli.guide_finder import GuideFinder


def resolve_existing_path(path_str: str, base_dir: Path) -> Path:
    p = Path(path_str)
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()


def parse_region(region: str) -> tuple[str, int, int]:
    try:
        chrom, coords = region.split(":")
        start_text, end_text = coords.split("-")
        start = int(start_text.replace(",", ""))
        end = int(end_text.replace(",", ""))
    except ValueError as exc:
        raise ValueError(
            f"Invalid region format: {region}. Expected format like chr17:43044295-43125483"
        ) from exc

    if start > end:
        raise ValueError(f"Invalid region: start ({start}) is greater than end ({end})")

    return chrom, start, end


def print_top_guides(guides: list[dict], top_n: int = 10) -> None:
    print(f"Top {min(top_n, len(guides))} ranked guides:")
    for g in guides[:top_n]:
        print(g)


def print_detailed_guides(guides: list[dict], count: int = 3) -> None:
    print(f"Detailed view of the top {min(count, len(guides))} guides:")
    for g in guides[:count]:
        print("-" * 60)
        print(f"Guide: {g.get('guide')}")
        print(f"Strand: {g.get('strand')}")
        print(f"PAM: {g.get('pam')}")
        print(
            f"Guide coordinates: "
            f"{g.get('guide_genomic_start')}-{g.get('guide_genomic_end')}"
        )
        print(f"Score: {g.get('score')}")
        print(f"GC%: {g.get('gc_percent')}")
        if "genome_match_count" in g:
            print(f"Genome match count: {g.get('genome_match_count')}")
        if "off_target_count" in g:
            print(f"Off-target count: {g.get('off_target_count')}")
        if "is_unique" in g:
            print(f"Is unique: {g.get('is_unique')}")
        print(f"Overlapping genes: {g.get('overlapping_genes')}")
        print(f"Hits target gene: {g.get('hits_target_gene')}")
        print(f"Target gene only: {g.get('target_gene_only')}")
        print(f"Non-target gene count: {g.get('non_target_gene_count')}")
        if "hits_target_cds" in g:
            print(f"Hits target CDS: {g.get('hits_target_cds')}")
        print(f"Target transcript count: {g.get('target_transcript_count')}")
        print(f"Target exon count: {g.get('target_exon_count')}")
        if "target_cds_count" in g:
            print(f"Target CDS count: {g.get('target_cds_count')}")
        print(f"Best target exon number: {g.get('best_target_exon_number')}")
        if "best_target_cds_exon_number" in g:
            print(f"Best target CDS exon number: {g.get('best_target_cds_exon_number')}")
        print(f"Overlapping transcripts: {g.get('overlapping_transcripts')}")
        print(f"Overlapping exon count: {g.get('overlapping_exon_count')}")
        if "overlapping_cds_count" in g:
            print(f"Overlapping CDS count: {g.get('overlapping_cds_count')}")
        print(f"Overlapping exons: {g.get('overlapping_exons')}")
        if "overlapping_cds" in g:
            print(f"Overlapping CDS: {g.get('overlapping_cds')}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Find CRISPR guides in a genomic region")
    parser.add_argument("region", help="Genomic region, e.g. chr17:43044295-43125483")
    parser.add_argument(
        "target_gene",
        nargs="?",
        default=None,
        help="Optional target gene name, e.g. BRCA1",
    )

    parser.add_argument("--genome-fa", default="data/genome.fa")
    parser.add_argument("--annotation-gtf", default="data/annotations.gtf")

    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--detail-n", type=int, default=3)

    parser.add_argument("--require-cds", action="store_true")
    parser.add_argument("--require-unique-exact", action="store_true")
    parser.add_argument("--max-1mm-offtargets", type=int, default=None)
    parser.add_argument("--deduplicate-exact-guides", action="store_true")
    parser.add_argument("--diverse-top", action="store_true")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parent.parent
    genome_path = resolve_existing_path(args.genome_fa, project_root)
    annotation_path = resolve_existing_path(args.annotation_gtf, project_root)

    if not genome_path.exists():
        raise FileNotFoundError(
            f"Genome FASTA not found: {genome_path}\n"
            f"Pass --genome-fa with the correct file path."
        )

    if not annotation_path.exists():
        raise FileNotFoundError(
            f"Annotation GTF not found: {annotation_path}\n"
            f"Pass --annotation-gtf with the correct file path."
        )

    chrom, start, end = parse_region(args.region)

    print("Loading genome...")
    genome_store = GenomeStore(str(genome_path))
    genome_store.load()
    print(f"Loaded chromosomes: {len(genome_store.chromosome_names())}")
    print(f"Chromosome examples: {genome_store.chromosome_names()[:10]}")

    print("Loading annotations...")
    annotation_store = AnnotationStore(str(annotation_path))

    print("Querying region...")
    print(f"Region: {args.region}")
    summary = annotation_store.summarize_region(chrom, start, end)
    print(f"Summary: {summary}")

    if args.target_gene:
        print(f"Target gene: {args.target_gene}")

    print("Getting sequence...")
    sequence = genome_store.get_sequence(chrom, start, end)
    print(f"Sequence length: {len(sequence)}")
    print(f"Sequence sample: {sequence[:60]}")

    finder = GuideFinder(genome_store, annotation_store)

    guides = finder.find_guides_in_region(
        chrom=chrom,
        region_start=start,
        region_end=end,
        target_gene=args.target_gene,
        require_target_gene=bool(args.target_gene),
        require_exon_overlap=True,
        require_cds_overlap=args.require_cds,
        max_one_mismatch_offtargets=args.max_1mm_offtargets,
        require_unique_exact=args.require_unique_exact,
        deduplicate_exact_sequences=args.deduplicate_exact_guides,
    )

    if args.diverse_top:
        guides = finder.select_diverse_top_guides(
            ranked_guides=guides,
            top_n=args.top_n,
        )

    print(f"Filtered guides: {len(guides)}")
    print_top_guides(guides, args.top_n)
    print_detailed_guides(guides, args.detail_n)


if __name__ == "__main__":
    main()