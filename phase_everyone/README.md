The following instructions are for phasing and getting the coverage of UCEs from the phyluce pipeline. They are a little tweaked from the main instructions, but the coverage.R file and phasing_settings file should be the same as in the main instructions, and being familiar with the main instructions would be a good idea. 

1) After converting the 75% incomplete dataset (WITHOUT designators inserted for missing taxa) into fasta (separate file per locus) using the convert formats tool from phyluce, getting rid of line breaks in fasta with the onelining.R code at: https://github.com/laninsky/ambigoos_into_structure/blob/master/onelining.R
```
for i in `ls *.fa*`; do mv $i temp; Rscript onelining.R; mv tempout $i; done;
```
2) Get a list of sample names and put them in file called "samplenames.txt" (one name per line) in the folder with a copy of your fasta files. Run the following (cribbed from step 2 of : https://github.com/laninsky/reference_aligning_to_established_loci). This will also give us a list of the loci which we can use in step 4.
```
nosamples=`wc -l samplenames.txt | awk '{print $1}'`;

for i in `seq 1 $nosamples`;
do refname=`tail -n+$i samplenames.txt | head -n1`;
mkdir $refname;
for j in `ls *.fasta`;
do if grep --quiet $refname $j
then echo ">"$j >> $refname/reference.fa;
temp=`grep -A1 $refname $j`;
echo $temp | awk '{print $2}' >> $refname/reference.fa;
fi
done;
done;
ls *.fasta > fasta_names
```
The *.fasta files for each locus (but not reference.fa files) can theb be removed.

3) Make sure the F and R reads are named the same within the READ1 and READ2 files within the cleaned reads folder - you only have to do this step if the reads have been prefaced with a 1: in the READ1 folder and a 2: in the READ2 folder (running this from inside the cleaned reads folder):
```
for i in `\ls`; do gunzip $i/split-adapter-quality-trimmed/*READ2*; done;
for i in `\ls`; do sed -i 's/ 2:/ 1:/g' $i/split-adapter-quality-trimmed/*READ2*; done;
for i in `\ls`; do gzip $i/split-adapter-quality-trimmed/*READ2*; done;
```

4) Run the UCE_phasing.sh script in this folder (make sure the coverage.R from the main folder is also in here, and phasing_settings etc has been set up as for the main instructions).

5) Run summarize_coverage.R from the main folder to obtain a locus-by-locus summary of coverage and sample size. The output file will give the locus name in the first column, how many samples that locus was present in, and then the average coverage across these samples. Rejoin the main set of scripts, resuming at step 5, at:
https://github.com/laninsky/reference_aligning_to_established_loci

Follow any instructions for UCE-related points.