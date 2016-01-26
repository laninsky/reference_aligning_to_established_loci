javapath=`tail -n+1 phasing_settings | head -n1`
gatk=`tail -n+2 phasing_settings | head -n1`
picard=`tail -n+3 phasing_settings | head -n1`
sequencing=`tail -n+4 phasing_settings | head -n1`

sed -i 's/\?/N/g' reference.fa
sed -i 's/-//g' reference.fa
bwa index -a is reference.fa
samtools faidx reference.fa
$javapath -jar $picard CreateSequenceDictionary R=reference.fa O=reference.dict

nosamples=`wc -l samples.txt | awk '{print $1}'`

for i in `seq 1 $nosamples`;
do name=`tail -n+$i samples.txt | head -n1`;
forward_proto=`tail -n+5 phasing_settings | head -n1`;
forward=`eval "echo $forward_proto"`;

if [ $sequencing == paired ]
then
reverse_proto=`tail -n+6 phasing_settings | head -n1`;
reverse=`eval "echo $reverse_proto"`;
bwa mem reference.fa $forward $reverse > temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem $i $forward > temp.sam;
fi

$javapath -jar $picard AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=phase SM=everyone;
$javapath -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE;
$javapath -jar $picard SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam;
samtools index tempsortmarked.bam;
$javapath -jar $gatk -T RealignerTargetCreator -R reference.fa -I tempsortmarked.bam -o tempintervals.list;
$javapath -jar $gatk -T IndelRealigner -R reference.fa -I tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
$javapath -jar $gatk -T DepthOfCoverage -R reference.fa -I temp_realigned_reads.bam -o temp.coverage
rm -rf temp.coverage.sample_*
echo $name > name
Rscript coverage.R


$javapath -jar $gatk -T HaplotypeCaller -R reference.fa -I temp_realigned_reads.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 -o temp_raw_variants.vcf;
$javapath -jar $gatk -T ReadBackedPhasing -R reference.fa -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -o temp_phased_SNPs.vcf;
$javapath -jar $gatk -T FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R reference.fa -o temp_alt.fa;

Rscript onelining.R;

rm -rf $name.*;
mv temp_alt2.fa $name.fa;
rm -rf temp*;

sed -i 's/\?/N/g' $i;
sed -i 's/-//g' $i;
bwa index -a is $i;
samtools faidx $i;
java -jar $picard CreateSequenceDictionary R=$i O=$name.dict;

if [ $sequencing == paired ]
then
bwa mem $i $forward $reverse > temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem $i $forward > temp.sam;
fi

java -jar $picard AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=phase SM=everyone;
java -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE;
java -jar $picard SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam;
samtools index tempsortmarked.bam;
java -jar $gatk -T RealignerTargetCreator -R $i -I tempsortmarked.bam -o tempintervals.list;
java -jar $gatk -T IndelRealigner -R $i -I  tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
java -jar $gatk -T HaplotypeCaller -R $i -I temp_realigned_reads.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 -o temp_raw_variants.vcf;
java -jar $gatk -T ReadBackedPhasing -R $i -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -o temp_phased_SNPs.vcf;
java -jar $gatk -T FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R $i -o temp_alt.fa;

Rscript onelining.R;

mv $i safe.$i.ref.fa
rm -rf $name.*;
mv temp_alt2.fa $name.1.fa;
mv safe.$i.ref.fa $name.2.fa
rm -rf temp*;

done

for i in `ls *1.fa`;
do sed -i 's/-//g' $i;
done
