







#split to small size file
split -d -l 1000 neardups-clusters_seq.tsv ./split/output --additional-suffix=.tsv

#using '--bastch' 
snakemake -s snkemake_dup --cores 4 --unlock
for BATCH in `seq 1 10`
do
    echo "Batch " $BATCH
#ls *cram | while read id; do samtools index $id ;done
    snakemake --scheduler greedy -s snkemake_dup --cores 4 --batch all=$BATCH/10
done
