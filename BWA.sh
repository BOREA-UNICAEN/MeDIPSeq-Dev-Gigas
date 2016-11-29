####### BWA alignment script

bwa index -a is c.gigas_genome.fa

bwa aln -t 16 -f sample1_forward_read.sai c.gigas_genome.fa sample1_forward_read.fq

bwa aln -t 16 -f sample1_reverse_read.sai c.gigas_genome.fa sample1_reverse_read.fq

bwa sampe -f sample1.sam c.gigas_genome.fa sample1_forward.sai sample1_reverse.sai sample1_forward.fq sample_1_reverse.fq

samtools view -bS sample1.sam > sample1.bam

samtools flagstat sample1.bam

 