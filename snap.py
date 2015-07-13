from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import pprint
import sys



## checks the original base call for each snp against the genome to make sure
## that the calls are the same. currently rejects any call that isn't correct
## though we can modify later to allow a % threshold of errors or soemthing
## along those lines.
def validate_genome(genome, snps) :
       status = True
       failed_snps = {}
       for key, val in snps.items() :
           if genome[int(key)] == val.split(":")[0] :
               status = True
           else :
               status = False
               failed_snps[key] = "Genome:" + genome[int(key)] + ",VCF_Call:" + val.split(":")[0]

       return status, failed_snps

## this function reads in a VCF file and parses out the needed information
## dummy function currently, just parses a simple tsv file
## will look for a module to parse VCF, otherwise we will have to write one
def vcf_read(vcf_file) :
    snps = {}
    print "Reading in VCF file %s..." %vcf_file
    handle = open(vcf_file, 'r')
    for line in handle :
        rec = line.strip().split()
        snps[rec[0]] = rec[1] + ":" + rec[2]

    return snps

## this function accepts a local genome file from genbank.
def genome_read(genome_file) :
    print "Reading in Genbank genome %s..." %genome_file
    genome = SeqIO.read(genome_file, 'genbank')

    return genome

def fetch_coding(genome, key, val) :
    result = []
    print "#"*100
    for g in genome.features :
        if g.type == "CDS" :
            start_codon = g.qualifiers['codon_start']
            base_positions = [ x for x in g.location ]
            gene_bases = []
            ### working. now need ot find out if the translation aspect with SNP exchange is working
            ### properly.
            if g.strand == -1 :
                for b in base_positions :
                    if genome.seq[b].upper() == 'A' :
                        gene_bases.append('T')
                    elif genome.seq[b].upper() == 'G' :
                        gene_bases.append('G')
                    elif genome.seq[b].upper() == 'T' :
                        gene_bases.append('A')
                    elif genome.seq[b].upper() == 'G' :
                        gene_bases.append('C')
                    else :
                        gene_bases.append(genome.seq[b])
            else :
                for b in base_positions :
                    gene_bases.append(genome.seq[b])
            if int(key) in base_positions :
                gene = g.qualifiers['gene']
                loc = key
                orig_dna = g.extract(genome.seq)
                orig_trans = orig_dna.translate()

                aa_translation = {}
                n = max(1, 3)
                aa_positions = [base_positions[i:i + n] for i in range(0, len(base_positions), n)]
                aa_index = 0
                for aa in orig_trans :
                    recs = aa_positions[aa_index]
                    for r in recs :
                        aa_translation[int(r)] = aa
                    aa_index += 1

                orig_aa = aa_translation[int(loc)]
                print orig_aa
                loc_index = base_positions.index(int(loc))
                new_dna_bases = gene_bases
                new_dna_bases[loc_index] = val.split(":")[-1]
                new_dna = Seq("".join(new_dna_bases))
                new_trans = new_dna.translate()
                
                print orig_trans
                print "\n"
                print new_trans

    print "#"*100

## main function. primarily used to parse command line ArgumentParser
## and call the worker methods
def main(args) :
   snps = vcf_read(args.vcf)
   genome = genome_read(args.genome)
   status, failed_snps = validate_genome(genome, snps)
   if status :
       print "All SNPS passed. Carry on."
       for key, val in snps.items() :
           fetch_coding(genome, key, val)

   else :
       print "Some SNPs failed validation. Will dump a file for further checking."
       print "%i / %i SNPs failed validation." %(len(failed_snps), len(snps))
       print "Please make sure you are using the correct genome."



if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description='SNAP: A tool to find changes in genes driven by SNP variation')
    parser.add_argument('-v', '--vcf', help='VCF input file.', required=True)
    parser.add_argument('-g', '--genome', help='Genome input file.', required=True)
    parser.add_argument('-o', '--out', help='Output file name.', required=True)
    args = parser.parse_args()
    main(args)
