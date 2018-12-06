javapath=`tail -n+1 phasing_settings | head -n1`;
gatk=`tail -n+2 phasing_settings | head -n1`;
picard=`tail -n+3 phasing_settings | head -n1`;
sequencing=`tail -n+4 phasing_settings | head -n1`;
gatk38=`tail -n+7 phasing_settings | head -n1`
numberofcores=`tail -n+8 phasing_settings | head -n1`;

nosamples=`wc -l samples.txt | awk '{print $1}'`;

for j in fasta_files/*;
do sed -i 's/\?/N/g' $j;
sed -i 's/-//g' $j;
bwa index -a is $j;
samtools faidx $j;
$javapath -jar $picard CreateSequenceDictionary R=$j O=$j.dict;

samplearray=( `seq 1 $nosamples` )

for i in "${samplearray[@]}";
do name=`tail -n+$i samples.txt | head -n1`;
forward_proto=`tail -n+5 phasing_settings | head -n1`;
forward=`eval "echo $forward_proto"`;

if [ $sequencing == paired ]
then
reverse_proto=`tail -n+6 phasing_settings | head -n1`;
reverse=`eval "echo $reverse_proto"`;
bwa mem -V -t $numberofcores $j $forward $reverse > $name.temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem -V -t $numberofcores $j $forward > $name.temp.sam;
fi

# Convert temp.sam to temp.bam, outputting only reads that have mapping quality
samtools view -@ $numberofcores -b -F 5 -T $j $name.temp.sam > $name.temp.bam
# Sorting bam and indexing resulting file
samtools sort -@ $numberofcores $name.temp.bam > $name.tempsorted.bam
samtools index -@ $numberofcores $name.tempsorted.bam
# Printing out pileup
samtools mpileup -f $j $name.tempsorted.bam > $name.pileup
rm $name.temp*

# Generating the fasta files based on the pileup file
Rscript filtering_fasta_on_pileup.R
rm $name.pileup

# Aligning the sequences in the output file
nolines=`wc -l ${name}_pileup.fasta | awk '{print $1}'`
if [ "$nolines" -gt "2" ]
then
echo $i
echo $j
mafft --thread $numberofcores ${name}_pileup.fasta > ${name}_pileup_aligned.fasta
break 2
else
locusname=`echo $j | sed 's"fasta_files/""g' | sed 's/\.fasta//g'`
echo ">"$locusname >> $name.reference.fa
tail -n+2 ${name}_pileup.fasta >> $name.reference.fa
fi

rm ${name}_pileup.fasta

done

done

