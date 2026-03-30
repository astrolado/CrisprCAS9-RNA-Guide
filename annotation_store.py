class AnnotationStore:
    def __init__(self, gtf_path):
        self.features = {}
        self.load_gtf(gtf_path)

    def normalize_chrom_name(self, chrom):
        chrom = chrom.strip()
        if not chrom:
            return chrom

        if chrom.startswith(">"):
            chrom = chrom[1:]

        if chrom.lower().startswith("chr"):
            base = chrom[3:]
        else:
            base = chrom

        if base.upper() in {"M", "MT"}:
            return "chrM"

        return f"chr{base}"

    def parse_attributes(self, attr_string):
        attrs = {}
        fields = attr_string.strip().split(";")

        for field in fields:
            field = field.strip()
            if not field:
                continue

            parts = field.split(" ", 1)
            if len(parts) != 2:
                continue

            key, value = parts
            attrs[key] = value.strip('"')

        return attrs

    def load_gtf(self, path):
        with open(path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) != 9:
                    continue

                chrom, source, feature, start, end, score, strand, frame, attributes = parts

                if feature not in {"gene", "transcript", "exon", "CDS"}:
                    continue

                chrom = self.normalize_chrom_name(chrom)
                attr_dict = self.parse_attributes(attributes)

                record = {
                    "feature": feature,
                    "feature_type": feature,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_id": attr_dict.get("gene_id"),
                    "gene_name": attr_dict.get("gene_name"),
                    "transcript_id": attr_dict.get("transcript_id"),
                    "exon_number": attr_dict.get("exon_number"),
                }

                if chrom not in self.features:
                    self.features[chrom] = []

                self.features[chrom].append(record)

    def query_region(self, chrom, start, end):
        chrom = self.normalize_chrom_name(chrom)
        matches = []

        for feature in self.features.get(chrom, []):
            if feature["start"] <= end and feature["end"] >= start:
                matches.append(feature)

        return matches

    def get_features(self, chrom, start, end):
        return self.query_region(chrom, start, end)

    def summarize_region(self, chrom, start, end):
        results = self.query_region(chrom, start, end)

        genes = set()
        transcripts = set()
        exon_count = 0
        cds_count = 0

        for feature in results:
            if feature["gene_name"]:
                genes.add(feature["gene_name"])
            if feature["transcript_id"]:
                transcripts.add(feature["transcript_id"])
            if feature["feature"] == "exon":
                exon_count += 1
            if feature["feature"] == "CDS":
                cds_count += 1

        return {
            "genes": sorted(genes),
            "transcript_count": len(transcripts),
            "exon_count": exon_count,
            "cds_count": cds_count,
            "feature_count": len(results),
        }