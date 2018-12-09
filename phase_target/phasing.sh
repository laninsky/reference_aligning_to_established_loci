javapath=`tail -n+1 phasing_settings | head -n1`;
gatk=`tail -n+2 phasing_settings | head -n1`;
picard=`tail -n+3 phasing_settings | head -n1`;
sequencing=`tail -n+4 phasing_settings | head -n1`;
gatk38=`tail -n+7 phasing_settings | head -n1`
numberofcores=`tail -n+8 phasing_settings | head -n1`;

echo "samplename" "locus" "ref_length(bp)" "bp_covered_by_seq" "min_cov" "max_cov" "mean_inc_0" "mean_exc_0" > coverage_summary.txt;

nosamples=`wc -l samples.txt | awk '{print $1}'`;

bysample () {
name=`tail -n+$i samples.txt | head -n1`;

sed -i 's/\?/N/g' ${name}.reference.fa;
sed -i 's/-//g' ${name}.reference.fa;
bwa index -a is ${name}.reference.fa;
samtools faidx ${name}.reference.fa;
$javapath -jar $picard CreateSequenceDictionary R=${name}.reference.fa; O=${name}.reference.fa.dict;

forward_proto=`tail -n+5 phasing_settings | head -n1`;
forward=`eval "echo $forward_proto"`;

if [ $sequencing == paired ]
then
reverse_proto=`tail -n+6 phasing_settings | head -n1`;
reverse=`eval "echo $reverse_proto"`;
bwa mem -t $numberofcores ${name}.reference.fa $forward $reverse > ${name}.temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem -t $numberofcores ${name}.reference.fa $forward > ${name}.temp.sam;
fi

samtools view -@ $numberofcores -b -F 5 -T ${name}.reference.fa $name.temp.sam > $name.temp.bam
# Sorting bam and indexing resulting file
samtools sort -@ $numberofcores $name.temp.bam > $name.tempsorted.bam
$javapath -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=$name.tempsorted.bam O=$name.tempsorteddups.bam M=temp.metrics AS=TRUE;
$javapath -jar $picard AddOrReplaceReadGroups I=$name.tempsorteddups.bam O=$name.tempsorteddupsrg.bam LB=rglib PL=illumina PU=phase SM=everyone;

samtools index -@ $numberofcores $name.tempsorteddupsrg.bam;

# nt/nct option doesn't work for this one
$javapath -jar $gatk38 -T DepthOfCoverage -R ${name}.reference.fa -I $name.tempsorteddupsrg.bam -o ${name}.temp.coverage;
rm -rf ${name}.temp.coverage.sample_*;

# Getting coverage for the sample
Rscript coverage.R ${name}.temp.coverage

# Phasing
$gatk/gatk HaplotypeCaller -R ${name}.reference.fa -I $name.tempsorteddupsrg.bam -stand-call-conf 30 -O $name.temp_raw_variants.vcf;

# nt/nct option doesn't work for these commands
$javapath -jar $gatk38 -T FindCoveredIntervals -R ${name}.reference.fa -I $name.tempsorteddupsrg.bam -cov 1 -o $name.temp_covered.list;
$javapath -jar $gatk38 -T FastaAlternateReferenceMaker -V $name.temp_raw_variants.vcf -R ${name}.reference.fa -L $name.temp_covered.list -o $name.temp_alt.fa;

Rscript modref.R;

mv temp_alt2.fa $name.fa;
rm -rf temp*;

sed -i 's/\?/N/g' $name.fa;
sed -i 's/-//g' $name.fa;
bwa index -a is $name.fa;
samtools faidx $name.fa;
$javapath -jar $picard CreateSequenceDictionary R=$name.fa O=$name.dict;

if [ $sequencing == paired ]
then
bwa mem -t $numberofcores $name.fa $forward $reverse > temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem -t $numberofcores $name.fa $forward > temp.sam;
fi

$javapath -jar $picard AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=phase SM=everyone;
$javapath -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE;
$javapath -jar $picard SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam;
samtools index tempsortmarked.bam;
# These tools were depreciated in gatk version 4
#$gatk/gatk -T RealignerTargetCreator -R $name.fa -I tempsortmarked.bam -o tempintervals.list;
#$gatk/gatk -T IndelRealigner -R $name.fa -I  tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
# The -stand_emit_conf 30 option is deprecated in GATK v 3.7 and was removed from this code on the 5-June-2017
$gatk/gatk HaplotypeCaller -R $name.fa -I tempsortmarked.bam -stand-call-conf 30 -O temp_raw_variants.vcf;
$javapath -jar $gatk38 -T FindCoveredIntervals -R $name.fa -I tempsortmarked.bam -cov 1 -o temp_covered.list;
$javapath -jar $gatk38 -T FastaAlternateReferenceMaker -V temp_raw_variants.vcf -R $name.fa -o temp_alt.fa;

Rscript modref.R;

mv $name.fa safe.$name.fa.ref.fa;
rm -rf $name.*;
mv temp_alt2.fa $name.1.fa;
mv safe.$name.fa.ref.fa $name.2.fa;
rm -rf temp*;

done

for i in `seq 1 $nosamples`;
do 

for i in `ls *1.fa`;
do sed -i 's/-//g' $i;
done
