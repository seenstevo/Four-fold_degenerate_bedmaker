#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 18:24:00 2017

@author: Sean Stevenson
"""

'''This script will take a fasta file created from the genome fasta and gff3 files described below. 
Codons are scanned and compared to the list of four-fold degenerate codon list. 
Matches are reported in bed format which point to the third position only but which also reports the full codon it is found in.
To generate the nessesary fasta file you can run:
awk '($3=="CDS") {print}' genes.gff | awk '{OFS="\t"; split($9, a, "="); split(a[2], b, "."); print $1, ($4-1), $5,b[1]"."b[2]"_"$1"_"$4}' \
| sort -k1,1 -k2,2n | bedtools getfasta -name -fi genome.fa -bed - -fo CDS.fasta
Note, the splitting may not work for all gff file formats depending on how the column of locus IDs is structured.
To run this script on command line use:
python FFDS_bedmaker.py [input.fasta]
This will print to STDOUT so can be piped to further commands (such as bedtools) or to file. 
'''     

from Bio.Seq import Seq
from Bio import SeqIO
import sys
        
fasta_in = sys.argv[1]
        
fourfold_deg = ["CTT", "CTA", "CTG", "CTC", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "GGT", "GGC", "GGA", "GGG"]

def codon_flavour_bed_maker (fasta):         
    old_left = ""
    old_x = 0
    old_header = ""
    for seq_record in SeqIO.parse(fasta, "fasta"):
        header_s = seq_record.name.split("::")[0].split("_")  # this splitting changes depending on the naming format of bedtools get fasts
        seq = str(seq_record.seq)
        Chr = header_s[1]
        strand = header_s[3]
        
        # Check whether we are still handling the same isomer to add previous trimmed bases to start
        if old_header == header_s[0]:
            main = old_left+seq
        else:
            main = seq
            old_x = 0
            
        # With main defined we can now trim the end to add to start of next cds
        x = len(main)%3
        if x!=0:
            left = main[-x:]
            main = main[:-x]
        else:
            main = main
            left = ""
            
        # Now we want to scan across the sequence to identify 4-fold degenerate sites
        for i in range(0, len(main), 3):
            # Adjust start depending on carry-over from previous
            adjusted_cds_start = int(header_s[2]) - int(old_x)
            if strand == "+":
                codon = main[i:i+3]
                fourfold_nuc_pos = adjusted_cds_start + (i+3)
            elif strand =="-":
                codon = Seq(main[i:i+3]).reverse_complement()
                codon = str(codon)
                fourfold_nuc_pos = adjusted_cds_start + i
            if codon in fourfold_deg:
                print("\t".join([Chr, str(fourfold_nuc_pos-1), str(fourfold_nuc_pos), header_s[0], codon]))
        # Set the current variables as 'old' for use with next seq
        old_header = header_s[0]
        old_left = left
        old_x = x
    
codon_flavour_bed_maker(fasta_in)


    
    
