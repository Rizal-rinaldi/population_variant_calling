import vcf
import pandas as pd
import csv
import argparse
import gzip


def parse_gff(gff_file):
    """
    Args:
        gff_file (str): Path to the GFF file.

    Returns:
        dict: Gene transcript count mapping.
        dict: Gene total transcripts mapping.
    """
    gff_gene_trans = {}
    transcript_geneid = {}
    with gzip.open(gff_file, 'rt') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'mRNA':
                attributes = fields[8].split(';')
                gene_id = None
                transcript_id = None
                for attr in attributes:
                    if attr.strip().startswith('Parent=gene:'):
                        gene_id = attr.strip()[12:]
                    elif attr.strip().startswith('ID=transcript:'):
                        transcript_id = attr.strip()[14:]
                if gene_id and transcript_id:
                    if gene_id in gff_gene_trans:
                        gff_gene_trans[gene_id].add(transcript_id)
                    else:
                        gff_gene_trans[gene_id] = {transcript_id}
                    transcript_geneid[transcript_id] = gene_id
    return gff_gene_trans, transcript_geneid

def calculate_frequency(record, samples_table):
    """
    Calculate frequency for each population in a VCF record.

    Args:
        record (vcf.Record): VCF record.
        samples_table (pd.DataFrame): DataFrame containing sample information.

    Returns:
        list: List of population frequencies.
    """
    pop_gt_dict = {}
    for sample in record.samples:
        population = samples_table.loc[samples_table["sample"].str.startswith(sample.sample), "population"].item()
        if population in pop_gt_dict:
            pop_gt_dict[population].append(sample["GT"])
        else:
            pop_gt_dict[population] = [sample["GT"]]
    frequencies = []
    for key, value in pop_gt_dict.items():
        points = 0
        for el in value:
            if el == "1/1":
                points += 2
            elif el == "1/0" or el == "0/1":
                points += 1
            elif el == "0/0":
                continue
        freq = points / (len(value) * 2) * 100
        frequencies.append(freq)
    return frequencies

def process_vcf_record(record, gff_gene_trans, transcript_geneid):
    """
    Process a VCF record and extract relevant information.

    Args:
        record (vcf.Record): VCF record (line in vcf file).
        gene_transcript_mapping (dict): Gene transcript count mapping.

    Returns:
        list: List of extracted information.
    """
    row_to_write = [record.CHROM, record.POS, record.INFO.get("TYPE"), record.QUAL]

    frequencies = calculate_frequency(record, samples_table)
    row_to_write.extend(frequencies)

    hom_samples, het_samples = [], []
    for hom in record.get_hom_alts():
        hom_samples.append(hom.sample)
    for het in record.get_hets():
        het_samples.append(het.sample)
    row_to_write.append(",".join(hom_samples))
    row_to_write.append(",".join(het_samples))

    ann, symbol, pc, str, tran, sub, transcripts = [], [], [], [], [], [], []
    csqs = record.INFO.get("BCSQ")
    gene_transcript_dict = {}
    if csqs is not None:
        for csq in csqs:
            elem = csq.split("|")
            for i in elem:
                if i.startswith('ENSMGAT'):
                    if elem[1] in gene_transcript_dict:
                        gene_transcript_dict[elem[1]].append(i)
                    else:
                        gene_transcript_dict[elem[1]] = [i]

            ann.append(elem[0])
            genes = elem[1] + ";" + elem[2]
            symbol.append(genes)
            pc.append(elem[3])
        for gene_id, transcripts_list in gene_transcript_dict.items():
            ensembl_id = transcript_geneid[transcripts_list[0]]
            gff_nr_trans = len(gff_gene_trans.get(ensembl_id, []))
            output = "{};{}".format(len(set(transcripts_list)), gff_nr_trans)
            transcripts.append(output)
        if len(elem) >= 5:
            str.append(elem[4])
        if len(elem) >= 6:
            tran.append(elem[5])
        if len(elem) >= 7:
            sub.append(elem[6])

    row_to_write.append(",".join(ann))
    row_to_write.append(",".join(symbol))
    row_to_write.append(",".join(pc))
    row_to_write.append(",".join(str))
    row_to_write.append(",".join(tran))
    row_to_write.append(",".join(sub))
    row_to_write.append(",".join(transcripts))
    return row_to_write



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF file")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--gff", required=True, help="Input GFF file")
    parser.add_argument("--sample", required=True, help="Input sample list")
    parser.add_argument("--out", required=True, help="Output TSV file")
    args = parser.parse_args()
    # Parse the GFF file 
    gene_id_to_transcript, transcript_geneid = parse_gff(args.gff)

    samples_table = pd.read_csv(args.sample, names=["sample", "bam", "population"])

    csvfile = open(args.out, "w")
    outfile = csv.writer(csvfile, delimiter='\t')

    populations = samples_table.population.unique()

    header = ["CHR", "POS", "TYPE", "QUAL"]
    population_headers = [p + "_freq" for p in populations]

    complete_header = header + population_headers + ["Homozygotes", "Heterozygotes", "Consequence", \
                                                     "Genes;transcripts", "Protein_coding", "Str", "AA_sub", "Nuc_sub",
                                                     "Nr transcripts vcf,total nr transcripts gene"]
    outfile.writerow(complete_header)
    reader = vcf.Reader(filename=args.vcf)
    for record in reader:
        row_to_write = process_vcf_record(record, gene_id_to_transcript, transcript_geneid)
        outfile.writerow(row_to_write)

    csvfile.close()
