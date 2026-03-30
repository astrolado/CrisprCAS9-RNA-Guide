import sys

from core.genome_store import GenomeStore
from core.annotation_store import AnnotationStore


def parse_region(region_string):
    chrom, coords = region_string.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def main():
    if len(sys.argv) != 2:
        print("Usage: python -m cli.main chr:start-end")
        return

    region_string = sys.argv[1]
    chrom, start, end = parse_region(region_string)

    print("Loading genome...")
    genome = GenomeStore("data/genome.fa")

    print("Loading annotations...")
    annot = AnnotationStore("data/annotation.gtf")

    print("Querying region...")
    results = annot.query_region(chrom, start, end)

    print(f"Region: {chrom}:{start}-{end}")
    print(f"Found features: {len(results)}")

    for r in results[:10]:
        print(r)

    print("Getting sequence...")
    seq = genome.get_sequence(chrom, start, end)
    print("Sequence length:", len(seq))
    print("Sequence sample:", seq[:50])

    summary = annot.summarize_region(chrom, start, end)
    print("Summary:", summary)

if __name__ == "__main__":
    main()