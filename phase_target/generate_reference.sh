javapath=`tail -n+1 phasing_settings | head -n1`;
sequencing=`tail -n+4 phasing_settings | head -n1`;
numberofcores=`tail -n+8 phasing_settings | head -n1`;
reference=`tail -n+9 phasing_settings | head -n1`;

nosamples=`wc -l samples.txt | awk '{print $1}'`;

sed -i 's/\?/N/g' $reference;
sed -i 's/-//g' $reference;
bwa index -a is $reference;
samtools faidx $reference;

# going to run samples in parallel by defining function bysample below
bysample () {
name=`tail -n+$i samples.txt | head -n1`;
forward_proto=`tail -n+5 phasing_settings | head -n1`;
forward=`eval "echo $forward_proto"`;

if [ $sequencing == paired ]
then
reverse_proto=`tail -n+6 phasing_settings | head -n1`;
reverse=`eval "echo $reverse_proto"`;
bwa mem -V -t $numberofcores $reference $forward $reverse > $name.temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem -V -t $numberofcores $reference $forward > $name.temp.sam;
fi

# Convert temp.sam to temp.bam, outputting only reads that have mapping quality
samtools view -@ $numberofcores -b -F 5 -T $reference $name.temp.sam > $name.temp.bam
# Sorting bam and indexing resulting file
samtools sort -@ $numberofcores $name.temp.bam > $name.tempsorted.bam
samtools index -@ $numberofcores $name.tempsorted.bam
# Printing out pileup
samtools mpileup -f $reference $name.tempsorted.bam > $name.pileup
rm $name.temp*

}

# for-loop allowing bysample to be parallelized
for i in `seq 1 $nosamples`;
do bysample "$i" & done
wait

# Generating the fasta files based on the pileup file
Rscript filtering_fasta_on_pileup.R
