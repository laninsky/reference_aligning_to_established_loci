javapath=`tail -n+1 phasing_settings | head -n1`;
gatk=`tail -n+2 phasing_settings | head -n1`;
picard=`tail -n+3 phasing_settings | head -n1`;
sequencing=`tail -n+4 phasing_settings | head -n1`;

nosamples=`wc -l samplenames.txt | awk '{print $1}'`;

echo "samplename" "locus" "ref_length(bp)" "bp_covered_by_seq" "min_cov" "max_cov" "mean_inc_0" "mean_exc_0" > coverage_summary.txt;

for i in `seq 1 $nosamples`;
do name=`tail -n+$i samplenames.txt | head -n1`;
cp $name/reference.fa ./;
sed -i 's/\?/N/g' reference.fa;
sed -i 's/-//g' reference.fa;
bwa index -a is reference.fa;
samtools faidx reference.fa;
$javapath -jar $picard CreateSequenceDictionary R=reference.fa O=reference.dict;
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
bwa mem reference.fa $forward > temp.sam;
fi

$javapath -jar $picard AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=phase SM=everyone;
$javapath -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE;
$javapath -jar $picard SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam;
samtools index tempsortmarked.bam;
# These tools were depreciated in gatk version 4
#$gatk/gatk -T RealignerTargetCreator -R reference.fa -I tempsortmarked.bam -o tempintervals.list;
#$gatk/gatk -T IndelRealigner -R reference.fa -I tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
$gatk/gatk DepthOfCoverage -R reference.fa -I temp_realigned_reads.bam -O temp.coverage;
rm -rf temp.coverage.sample_*;
echo $name > name;
Rscript coverage.R;

# The -stand_emit_conf 30 option is deprecated in GATK v 3.7 and was removed from this code on the 5-June-2017
$gatk/gatk HaplotypeCaller -R reference.fa -I temp_realigned_reads.bam -stand-call-conf 30 -O temp_raw_variants.vcf;
$gatk/gatk ReadBackedPhasing -R reference.fa -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -O temp_phased_SNPs.vcf;
$gatk/gatk FindCoveredIntervals -R reference.fa -I temp_realigned_reads.bam -cov 1 -O temp_covered.list;
$gatk/gatk FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R reference.fa -L temp_covered.list -O temp_alt.fa;

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
bwa mem $name.fa $forward $reverse > temp.sam;
fi

if [ $sequencing == single ]
then
bwa mem $name.fa $forward > temp.sam;
fi

$javapath -jar $picard AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=phase SM=everyone;
$javapath -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE;
$javapath -jar $picard SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam;
samtools index tempsortmarked.bam;
# These tools were depreciated in gatk version 4
#$gatk/gatk -T RealignerTargetCreator -R $name.fa -I tempsortmarked.bam -o tempintervals.list;
#$gatk/gatk -T IndelRealigner -R $name.fa -I  tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
# The -stand_emit_conf 30 option is deprecated in GATK v 3.7 and was removed from this code on the 5-June-2017
$gatk/gatk HaplotypeCaller -R $name.fa -I temp_realigned_reads.bam -stand-call-conf 30 -O temp_raw_variants.vcf;
$gatk/gatk ReadBackedPhasing -R $name.fa -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -O temp_phased_SNPs.vcf;
$gatk/gatk FindCoveredIntervals -R $name.fa -I temp_realigned_reads.bam -cov 1 -O temp_covered.list;
$gatk/gatk FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R $name.fa -O temp_alt.fa;

Rscript modref.R;

mv $name.fa safe.$name.fa.ref.fa;
rm -rf $name.*;
mv temp_alt2.fa $name.1.fa;
mv safe.$name.fa.ref.fa $name.2.fa;
rm -rf temp*;
rm -rf reference.*;
done

for i in `ls *1.fa`;
do sed -i 's/-//g' $i;
done
