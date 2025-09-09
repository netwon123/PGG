###extact genomoc feature, you need to annotate the version of gtf

1. used edited gtftools to extact the intro, exon, intergenic.(0-base)

2. used TBtools to extract the all information of genome 

3. awk '$3==""three_utr' tbtools.txt > 3utr.bed
   awk '$5==""five_utr' tbtools.txt >  5utr.bed




### annotate the feature in vcf
#merge and sort all bed file
for f in *.bed; do
    prefix=$(basename "$f" .bed)
    awk -v p="$prefix" '{print $0"\t"p}' "$f"
done > merged.bed
bedtools sort merged.bed -faidx genome.faidx > merge.sort.bed

bgzip merge.sort.bed
tabix -p bed merge.sort.bed

#https://github.com/brentp/vcfanno
vcfanno -ends conf.toml query.vcf.gz > query_anno.vcf


#extact the info
python extract_svlen.py panel_sv_anno.vcf panel_sv_anno.txt --extra-info-keys gene,feature
