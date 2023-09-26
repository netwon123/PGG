







#split to small size file
<code>split -d -l 1000 neardups-clusters_seq.tsv ./split/output --additional-suffix=.tsv<code>

#using '--bastch' 
<code>snakemake -s snkemake_dup --cores 4 --unlock<code>
<code>for BATCH in `seq 1 10`
do
    echo "Batch " $BATCH
#ls *cram | while read id; do samtools index $id ;done
    snakemake --scheduler greedy -s snkemake_dup --cores 4 --batch all=$BATCH/10
done<code>


<code>singularity pull docker://jmcdani20/hap.py:v0.3.12<code>

##evaluated the performance of 
<code>ls *pan.vcf.gz | while read id; do singularity exec --network=none -n ~/hap.py_latest.sif /opt/hap.py/bin/hap.py 1_ms.vcf.gz $id -r ../../Duroc.fa --threads 4 -o ${id} --unhappy --no-roc --no-json --engine=vcfeval;done<code>
