# Four-fold_degenerate_bedmaker
A python script to identify four-fold degenerate sites from a CDS fasta

To generate the required [CDS.fasta] file the following unix command should be run:

grep CDS [genes.gff] | awk '{OFS="\t"; split($9, a, "="); split(a[2], b, "."); print $1, ($4-1), $5,b[1]"."b[2]"_"$1"_"$4}' | sort -k1,1 -k4,4 | bedtools getfasta -name -fi [genome.fa] -bed - -fo [CDS.fasta]

where [genes.gff] is the gene annotation file in gff3 format build from the same genome assembly as [genome.fa]. 

Once [CDS.fasta] is made, the program should be run from the command line from where it can be piped to further commands or to file:

python FFDS_bedmaker.py [CDS.fasta]
