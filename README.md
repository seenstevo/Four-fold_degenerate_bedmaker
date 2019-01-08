# Four-fold_degenerate_bedmaker
A python script to identify four-fold degenerate sites from a CDS fasta

To generate the required [CDS.fasta] file the following unix command should be run:

awk '($3=="CDS") {print}' [genes.gff] | awk '{OFS="\t"; split($9, a, ";"); split(a[1], b, ":"); split(b[2], c, "_"); print $1, ($4-1), $5,c[1]"."c[2]"_"$1"_"$4"_"$7}' | sort -k1,1 -k4,4 -k2,2 | bedtools getfasta -name -fi [genome.fa] -bed - -fo [CDS.fasta]

Where [genes.gff] is the gene annotation file in gff3/gft format build from the same genome assembly as [genome.fa]. NOTE that the splitting may vary depending on the structure of the gff3/gft file column 9 attribute column. These specific steps relate to the B73_RefGen_v4 Zea mays assembly.

Once [CDS.fasta] is made, the program should be run from the command line from where it can be piped to further commands or to file:

python FFDS_bedmaker.py [CDS.fasta]
