from core.genome_store import GenomeStore
from core.annotation_store import AnnotationStore

print("Loading genome...")
genome = GenomeStore("data/genome.fa")

print("Loading annotations...")
annot = AnnotationStore("data/annotation.gtf")

print("Querying region...")
results = annot.query_region("chr17", 43044295, 43125483)

print("Found features:", len(results))

for r in results[:5]:
    print(r)

print("Getting sequence...")
seq = genome.get_sequence("chr17", 43044295, 43044350)
print("Sequence sample:", seq[:50])