#AUTHENTICATION WORKFLOW FOR BOWTIE2
#Example command line: 
#./authentic_bowtie.sh otzi.trimmed_BOWTIE2_FULL_NT.sam

TAXID=$1
SAM=$2

ml bioinfo-tools samtools mapDamage seqtk
mkdir $PWD/valid_${TAXID}
printf "\n"

echo "Grepping NCBI accession IDs for given TaxID = $TAXID"
grep -w $TAXID BOWTIE2_Full_NT_Database/20210127-003704_nt/seqid2taxid.map | cut -f1 > $PWD/valid_${TAXID}/${TAXID}.seq.ids

echo "Using NCBI accession IDs in order to grep reads from given SAM-file"
grep -wFf $PWD/valid_${TAXID}/${TAXID}.seq.ids $SAM > $PWD/valid_${TAXID}/${TAXID}.sam

echo "Select only reads with MAPQ >= 1 shorter than 100nt and make BAM-file"
awk 'length($10)<100' $PWD/valid_${TAXID}/${TAXID}.sam > $PWD/valid_${TAXID}/${TAXID}.pruned.sam
samtools view -bS $PWD/valid_${TAXID}/${TAXID}.pruned.sam -q 1 -h > $PWD/valid_${TAXID}/${TAXID}.bam

echo "Running PMDtools on BAM-file for given TaxID"
samtools view $PWD/valid_${TAXID}/${TAXID}.bam | python2 PMDtools/pmdtools.0.60.py --platypus --requirebaseq 30 --requiremapq 1 --maxlength 100 --minlength 30 --number 100000 > $PWD/valid_${TAXID}/PMD_temp.txt
cd $PWD/valid_${TAXID}
R CMD BATCH PMDtools/plotPMD.v2.R
cd ..

if [ "$(samtools view $PWD/valid_${TAXID}/${TAXID}.bam | wc -l)" -lt "100" ]; 
then
    echo "There are fewer than 100 reads with MAPQ >= 20 and shorter than 100 nt"
    exit
fi

echo "Extracting reference fasta-sequences for NCBI accession IDs"
seqtk subseq BOWTIE2_Full_NT_Database/20210127-003704_nt/library/nt/library.fna $PWD/valid_${TAXID}/${TAXID}.seq.ids > $PWD/valid_${TAXID}/${TAXID}.fasta

echo "Running mapDamage on BAM-file for given TaxID = $TAXID"
mapDamage -i $PWD/valid_${TAXID}/${TAXID}.bam -r $PWD/valid_${TAXID}/${TAXID}.fasta --merge-reference-sequences --no-stats -d $PWD/valid_${TAXID}/results_${TAXID}

printf "\n"
