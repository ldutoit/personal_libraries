#/bin/sh
infile=I.gatk.allsitesRM_sorted_chromcoord.vcf
awk -F "\t" '{print > $1".tosort"}' $infile
for i in *.tosort; do  echo $i ; sort -k2,2n $i > $i.s & done
sort -m -k1,1 -k2,2n *s > $infile_sorted.vcf
rm *.tosort*