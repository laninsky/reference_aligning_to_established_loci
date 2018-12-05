# Phase_target

You'll have a fasta file with the target loci that you are interested in pulling down from your sample. This might be a file that doesn't have multiple taxa (i.e. each target locus is present just once), or it may have multiple taxa per locus. This first little script is going to split out the fasta file into locus specific files.

### Getting your locus-specific fasta files together
First, define your file name
```
fastafilename=target_Zhang_et_al.fa
```
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
i=$(expr $i + 1);
fi
done;

echo $fasta_sequence | awk '{print $1}' >> fasta_files/$locus_name.fasta
```
