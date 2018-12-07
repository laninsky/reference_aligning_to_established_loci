

# Aligning the sequences in the output file
nolines=`wc -l ${name}_pileup.fasta | awk '{print $1}'`
if [ "$nolines" -gt "2" ]
then
echo $i
mafft --thread $numberofcores ${name}_pileup.fasta > ${name}_pileup_aligned.fasta
break 2
else
locusname=`echo $j | sed 's"fasta_files/""g' | sed 's/\.fasta//g'`
echo ">"$locusname >> $name.reference.fa
tail -n+2 ${name}_pileup.fasta >> $name.reference.fa
fi

rm ${name}_pileup.fasta


These steps are the same regardless of whether your target fasta has the same locus more than once or not
```
nolines=`wc -l $allelefilename | awk '{print $1}'`
mkdir fasta_files
```
For files where the target locus is present just once
```
headercontents=`tail -n+1 $allelefilename | head -n1`
locus_name=`echo $headercontents | sed 's/>//g'`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g'`
echo $linecontents >> fasta_files/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta
```
For files where the target locus is present more than once you'll need to define how the sequence header combines the 
the sample name and the gene/locus name e.g. what delimits them, and what is the first field
```
delimiter="-"
first_field="sample"
# or alternatively first_field="locus"

if [ $first_field = "sample" ]; then

headercontents=`tail -n+1 $allelefilename | head -n1`
locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta

else

headercontents=`tail -n+1 $allelefilename | head -n1`
locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j $allelefilename | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
echo $headercontents >> fasta_files/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

fi
```
