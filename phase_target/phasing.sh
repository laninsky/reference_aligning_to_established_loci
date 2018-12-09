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
$javapath -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=$name.tempsorted.bam O=$name.tempsorteddups.bam M=$name.temp.metrics AS=TRUE;
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

Rscript modref.R $name.temp_alt.fa

mv $name.temp.reference.fa $name.2.reference.fa;
rm -rf $name.temp*;
rm -rf $name.reference.fa.*
rm -rf $name.reference.dict

sed -i 's/\?/N/g' $name.2.reference.fa;
sed -i 's/-//g' $name.2.reference.fa;
bwa index -a is $name.2.reference.fa;
samtools faidx $name.2.reference.fa;
$javapath -jar $picard CreateSequenceDictionary R=$name.2.reference.fa O=$name.2.reference.fa.dict;

if [ $sequencing == paired ]
then
bwa mem -t $numberofcores ${name}.2.reference.fa $forward $reverse > ${name}.temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem -t $numberofcores ${name}.2.reference.fa $forward > ${name}.temp.sam;
fi

samtools view -@ $numberofcores -b -F 5 -T ${name}.2.reference.fa $name.temp.sam > $name.temp.bam
# Sorting bam and indexing resulting file
samtools sort -@ $numberofcores $name.temp.bam > $name.tempsorted.bam
$javapath -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=$name.tempsorted.bam O=$name.tempsorteddups.bam M=$name.temp.metrics AS=TRUE;
$javapath -jar $picard AddOrReplaceReadGroups I=$name.tempsorteddups.bam O=$name.tempsorteddupsrg.bam LB=rglib PL=illumina PU=phase SM=everyone;

samtools index -@ $numberofcores $name.tempsorteddupsrg.bam;

mv $name.2.reference.fa.dict $name.2.reference.dict

# Phasing
$gatk/gatk HaplotypeCaller -R ${name}.2.reference.fa -I $name.tempsorteddupsrg.bam -stand-call-conf 30 -O $name.temp_raw_variants.vcf;

# nt/nct option doesn't work for these commands
$javapath -jar $gatk38 -T FindCoveredIntervals -R ${name}.2.reference.fa -I $name.tempsorteddupsrg.bam -cov 1 -o $name.temp_covered.list;
$javapath -jar $gatk38 -T FastaAlternateReferenceMaker -V $name.temp_raw_variants.vcf -R ${name}.2.reference.fa -L $name.temp_covered.list -o $name.temp_alt.fa;

Rscript modref.R $name.temp_alt.fa

mv $name.temp.reference.fa $name.1.reference.fa;
rm -rf $name.temp*;
rm -rf $name.2.reference.fa.*
rm -rf $name.2.reference.dict

sed -i 's/\?/N/g' $name.1.reference.fa;
sed -i 's/-//g' $name.1.reference.fa;
}

# for-loop allowing bysample to be parallelized
for i in `seq 1 $nosamples`;
do bysample "$i" & done
wait

# feeding the per-sample coverage summary into the total table
for i in `seq 1 $nosamples`;
do name=`tail -n+$i samples.txt | head -n1`;
cat $name.coverage_summary.txt >> coverage_summary.txt;
rm $name.coverage_summary.txt;
done

