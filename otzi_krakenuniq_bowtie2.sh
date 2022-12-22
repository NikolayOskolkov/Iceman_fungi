#BUILDING KRAKENUNIQ FULL NT DATABASE ON VERY FAT NODE
module load bioinfo-tools perl
module load jellyfish/1.1.11
time krakenuniq-build --db DBDIR_KrakenUniq_Full_NT --kmer-len 31 --threads 80 --taxids-for-genomes --taxids-for-sequences

#RUNNING KRAKENUNIQ ON FULL NT DATABASE
ml bioinfo-tools perl
time krakenuniq --db DBDIR_KrakenUniq_Full_NT --fastq-input otzi.trimmed.fastq.gz --threads 80 --output sequences.krakenuniq_Full_NT --report-file krakenuniq.output_Full_NT --gzip-compressed --only-classified-out



#INDEX BOWTIE2 FULL NT DATABASE (was possible to do on a 1TB node)
ml bioinfo-tools bowtie2 perl
cd BOWTIE2_Full_NT_Database/20210127-003704_nt/library/nt
bowtie2-build --large-index library.fna library.fna --threads 20

#RUN BOWTIE ON FULL NT DATABASE
ml bioinfo-tools bowtie2 perl
time bowtie2 --large-index -x BOWTIE2_Full_NT_Database/20210127-003704_nt/library/nt/library.fna --end-to-end --threads 20 --very-sensitive -U otzi.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 10 - > otzi.trimmed_BOWTIE2_FULL_NT.bam

