<strong>simulated different coverage {5x,10x,15x,20,25x}</strong>
<p>
<code>   vg paths -x PGG_all.giraffe.gbz -F -R -Q 1 > 1.fa<code>
<code> mutation-simulator ../1.fa args -sn 0.005 -in 0.001 -inmin 2 -inmax 49 -de 0.001 -demin 2 -demax 49 -n benchmark -s pig<code>
mutation-simulator  ../1.fa args -sn 0.005 -in 0.001 -inmin 50 -inmax 100000 -de 0.001 -demin 50 -demax 100000 -n benchmark -s pig  -iv 0.0001913997 -ivmin 50 -ivmax 100000 -du 0.0001913997 -dumin 50 -dumax 100000 -tl 0.0001913997 -tlmin 50 -tlmax 100000

randomreads.sh ref=1_ms.fa paired=true q=40 length=150 coverage=5 maxq=40 midq=35 minq=25 out1=5x_R1.fastq.gz out2=5x_R2.fastq.gz seed=0 illuminanames=t replacenoref=t addpairnum=t
randomreads.sh ref=1_ms.fa paired=true q=40 length=150 coverage=10 maxq=40 midq=35 minq=25 out1=10x_R1.fastq.gz out2=10x_R2.fastq.gz seed=0 illuminanames=t replacenoref=t addpairnum=t
randomreads.sh ref=1_ms.fa paired=true q=40 length=150 coverage=15 maxq=40 midq=35 minq=25 out1=15x_R1.fastq.gz out2=15x_R2.fastq.gz seed=0 illuminanames=t replacenoref=t addpairnum=t
randomreads.sh ref=1_ms.fa paired=true q=40 length=150 coverage=20 maxq=40 midq=35 minq=25 out1=20x_R1.fastq.gz out2=20x_R2.fastq.gz seed=0 illuminanames=t replacenoref=t addpairnum=t
randomreads.sh ref=1_ms.fa paired=true q=40 length=150 coverage=25 maxq=40 midq=35 minq=25 out1=25x_R1.fastq.gz out2=25x_R2.fastq.gz seed=0 illuminanames=t replacenoref=t addpairnum=t


