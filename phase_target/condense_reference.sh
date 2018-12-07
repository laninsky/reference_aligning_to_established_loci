numberofcores=`tail -n+8 phasing_settings | head -n1`;
nosamples=`wc -l samples.txt | awk '{print $1}'`;
delimiter=`tail -n+10 phasing_settings | head -n1`;
locusfield=`tail -n+11 phasing_settings | head -n1`;

# going to run samples in parallel by defining function bysample below
bysample () {
name=`tail -n+$i samples.txt | head -n1`;
mkdir $name
nolines=`wc -l ${name}_pileup.fasta | awk '{print $1}'`
headercontents=`tail -n+1 ${name}_pileup.fasta | head -n1`

if [ $locusfield == "second" ]; then
locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
echo $headercontents >> $name/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j ${name}_pileup.fasta | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> $name/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\2/g"`
headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\1/g"`
echo $headercontents >> $name/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

echo $fasta_sequence | awk '{print $1}' >> $name/$locus_name.fasta

else
locus_name=`echo $headercontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
headercontents=`echo $headercontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
echo ">"$headercontents >> $name/$locus_name.fasta
fasta_sequence=""

for j in `seq 2 $nolines`;
do linecontents=`tail -n+$j ${name}_pileup.fasta  | head -n1`;
if echo $linecontents | grep --quiet ">";
then echo $fasta_sequence | awk '{print $1}' >> $name/$locus_name.fasta;
locus_name=`echo $linecontents | sed 's/>//g' | sed -E "s/(.*)$delimiter(.*)/\1/g"`
headercontents=`echo $linecontents | sed -E "s/(.*)$delimiter(.*)/\2/g"`
echo ">"$headercontents >> $name/$locus_name.fasta
fasta_sequence=""
else
fasta_sequence="$fasta_sequence$linecontents"
fi
done;

fi

# Aligning the sequences in the output file
for j in ${name}/*.fasta;

nolines=`wc -l $j | awk '{print $1}'`
if [ "$nolines" -gt "2" ]
then
mafft --thread $numberofcores $j > $j.aligned
break 2
else
locusname=`echo $j | sed "s'$name/''g" | sed 's/\.fasta//g'`
echo ">"$locusname >> $name.reference.fa
tail -n+2 $j >> $name.reference.fa
fi

rm ${name}_pileup.fasta





}

# for-loop allowing bysample to be parallelized
for i in `seq 1 $nosamples`;
do bysample "$i" & done
wait

