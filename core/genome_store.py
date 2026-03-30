class GenomeStore:
    def __init__(self, fasta_path):
        self.sequences = {}
        self.load_fasta(fasta_path)

    def load_fasta(self, path):
        with open(path, "r") as f:
            chrom = None
            seq = []

            for line in f:
                if line.startswith(">"):
                    if chrom:
                        self.sequences[chrom] = "".join(seq)
                    chrom = line.strip().split()[0][1:]
                    seq = []
                else:
                    seq.append(line.strip())

            if chrom:
                self.sequences[chrom] = "".join(seq)

    def get_sequence(self, chrom, start, end):
        return self.sequences[chrom][start - 1:end]