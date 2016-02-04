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
$javapath -jar $gatk -T RealignerTargetCreator -R reference.fa -I tempsortmarked.bam -o tempintervals.list;
$javapath -jar $gatk -T IndelRealigner -R reference.fa -I tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
$javapath -jar $gatk -T DepthOfCoverage -R reference.fa -I temp_realigned_reads.bam -o temp.coverage;
rm -rf temp.coverage.sample_*;
echo $name > name;
Rscript coverage.R;
rm reference.*
done;
