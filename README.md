# Four-fold_degenerate_bedmaker
A python script to identify four-fold degenerate sites from a CDS fasta

To generate the required [CDS.fasta] file the following unix command should be run:

awk '($3=="CDS") {print}' [genes.gff] | awk '{OFS="\t"; split($9, a, "="); split(a[2], b, "."); print $1, ($4-1), $5,b[1]"."b[2]"_"$1"_"$4}' | sort -k1,1 -k4,4 | bedtools getfasta -name -fi [genome.fa] -bed - -fo [CDS.fasta]

Where [genes.gff] is the gene annotation file in gff3/gft format build from the same genome assembly as [genome.fa]. NOTE that the splitting may vary depending on the structure of the gff3/gft file column 9 attribute column. This template is extracting the ID tag and retaining the full isoform value as based on the Arabidopsis thaliana TAIR10 gff3 file.

Once [CDS.fasta] is made, the program should be run from the command line from where it can be piped to further commands or to file:

python FFDS_bedmaker.py [CDS.fasta]
