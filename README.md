##WIP##

# Local Python CLI for CRISPR guide RNA design: scans genomic regions, annotates gene/CDS overlap, evaluates off-targets, and ranks candidates.

To utilize, access 

https://www.gencodegenes.org/human/

Then download the Basic Gene Annotation (PRI) and Genome sequence, primary assembly (GRCh38) (PRI).
Add these two files to a file named Data and run in cmd prompt CLI/main.py

## Why GenomeEngine?

CRISPR guide RNA design is not just about finding PAM sites. Effective guides must:
- target functionally relevant regions (e.g. coding sequence)
- minimize off-target effects
- maintain strong sequence properties for Cas9 activity

GenomeEngine provides a local, reproducible workflow to systematically evaluate and rank guide candidates using genomic data and annotation-aware scoring.


## Features

### Guide Discovery
- Scans genomic regions for SpCas9 (NGG) guides
- Supports both forward and reverse strands

### Annotation Awareness
- Identifies overlap with genes, transcripts, exons, and CDS
- Supports target-gene-specific filtering (e.g. BRCA1)

### Scoring System
- Prioritizes CDS-targeting guides
- Considers transcript coverage
- Evaluates GC content and sequence quality
- Penalizes homopolymers and poly-T sequences
- Includes off-target penalties (exact + 1 mismatch)

### Off-Target Evaluation
- Genome-wide exact match counting
- 1-mismatch off-target detection
- Uniqueness filtering

### Filtering & Ranking
- Require CDS overlap
- Require unique guides
- Limit mismatch off-targets
- Deduplicate identical sequences
- Diversity selection across genomic positions

### CLI Workflow
- Fully local execution
- Flexible parameter control





## Example: BRCA1 Guide Design

```bash
python -m genome_engine.cli.main \
  --genome data/genome.fa \
  --annotation data/annotations.gtf \
  --region chr17:43044295-43125483 \
  --target-gene BRCA1 \
  --require-cds \
  --require-unique-exact \
  --max-1mm-offtargets 0 \
  --diverse-top 10


```md
### Sample Output

| Rank | Guide Sequence | Strand | CDS Overlap | Off-targets (1mm) | Score |
|------|----------------|--------|-------------|-------------------|-------|
| 1    | GAGTCC...      | +      | Yes         | 0                 | 212   |
| 2    | CCTGGA...      | -      | Yes         | 0                 | 205   |
