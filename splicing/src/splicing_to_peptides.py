#!/usr/bin/env python3
import argparse
import gffutils
from pyfaidx import Fasta
from Bio.Seq import Seq
import csv
import sys

def build_gtf_db(gtf_path, db_path="hg.db"):
    return gffutils.create_db(
            gtf_path, dbfn=db_path, force=True,
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True
    )

def load_genome(fasta_path):
    return Fasta(fasta_path)

def extract_cdna(genome, chrom, strand, coords):
    seqs = [ genome[chrom][start:end].seq for start,end in coords ]
    cdna = "".join(seqs)
    if strand == "-":
        cdna = str(Seq(cdna).reverse_complement())
    return cdna

def best_orf(cdna):
    best = Seq("").translate()  # empty
    for frame in (0,1,2):
        pep = Seq(cdna[frame:]).translate(to_stop=True)
        if len(pep) > len(best):
            best = pep
    return str(best)

def tile_junction_peptides(prot, exonA_len, kmin=8, kmax=12):
    junction_pos = exonA_len
    peptides = []
    for k in range(kmin, kmax+1):
        for i in range(max(0, junction_pos-k+1), min(len(prot)-k+1, junction_pos+1)):
            if i < junction_pos < i+k:
                peptides.append(prot[i:i+k])
    return peptides

def main():
    p = argparse.ArgumentParser(
            description="Map peptides to rMATS skipped-exon variants"
    )
    p.add_argument("--gtf", required=True, help="GTF annotation")
    p.add_argument("--fasta", required=True, help="Genome FASTA")
    p.add_argument("--variants", required=True,
                   help="tab-delimited file with chr,strand,exonStart_0base,exonEnd,upstreamES,upstreamEE,downstreamES,downstreamEE,geneID,geneSymbol")
    p.add_argument("--output", required=True, help="output TSV of peptides")
    args = p.parse_args()

    # 1) load resources
    db = build_gtf_db(args.gtf)
    genome = load_genome(args.fasta)

    # 2) process each event
    writer = csv.writer(open(args.output, "w"), delimiter="\t", lineterminator="\n")
    writer.writerow([
        "geneID","geneSymbol","eventID","peptide","peptide_length"
    ])

    with open(args.variants) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for eid, row in enumerate(reader, start=1):
            chrom = row["chr"]
            strand = row["strand"]
            # 0-based exon coordinates → convert to Python slicing as-is
            exA = (int(row["upstreamES"]), int(row["upstreamEE"]))
            exB = (int(row["exonStart_0base"]), int(row["exonEnd"]))
            exC = (int(row["downstreamES"]), int(row["downstreamEE"]))

            # --- 1) Build the skipped‐exon (A→C) isoform only ---
            coords_skip = [ exA, exC ]
            cdna_skip = extract_cdna(genome, chrom, strand, coords_skip)

            # 2) Translate to protein (longest ORF)
            prot_skip = best_orf(cdna_skip)

            # 3) Output the full‐length skipped‐isoform protein to FASTA
            fasta_name = f"{row['geneSymbol']}_SE.fasta"
            header = f">{row['geneSymbol']}_{row['geneID']}_{row['exonStart_0base']}\n"
            with open(fasta_name, "w") as fa:
                fa.write(header)
                fa.write(prot_skip + "\n")

            # 4) Tile 8–12mer peptides across the novel A→C junction
            # exonA_len = length of exon A in aa
            exonA_len = len(Seq(
                extract_cdna(genome, chrom, strand, [exA])
            ).translate())
            peptides_skip = tile_junction_peptides(prot_skip, exonA_len, 8, 12)

            # 5) Write those junction‐spanning peptides to a second FASTA
            pep_fname = f"{row['geneSymbol']}_SE_peptides.fasta"
            with open(pep_fname, "w") as pf:
                for i, pep in enumerate(sorted(set(peptides_skip)), start=1):
                    pf.write(f">Peptide{i}\n{pep}\n")
                
    print(f"Done — peptides written to {args.output}")

if __name__=="__main__":
    main()
