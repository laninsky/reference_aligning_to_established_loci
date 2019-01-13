# phase_everyone v0.0
The following instructions were originally designed for phasing samples and getting the coverage of UCEs from the phyluce pipeline. They are a little tweaked from the main instructions, but the coverage.R file and phasing_settings file should be the same as in the main instructions, and being familiar with the main instructions would be a good idea. These can be adapted to any source of loci as long as you have a separate fasta file for each locus, with each sample present a maximum of once (i.e. not broken down by alleles). If you are coming from pyRAD, you can use the scripts in helper_scripts to get your fasta files together from the *.loci file.

1) After converting the incomplete UCE dataset (following gblocks etc etc) into fasta (separate file per locus, maximum one sequence per sample e.g. not phased alleles) using the convert formats tool from phyluce, get rid of line breaks in the fasta files with the onelining_firststep.R code within this repository (taken from https://github.com/laninsky/ambigoos_into_structure/blob/master/onelining.R). Note, this onelining code is different to that in the main folder.
```
for i in `ls *.fa*`; do mv $i temp; Rscript onelining_firststep.R; mv tempout $i; done;
rm temp
```
2) Get a list of sample names and put them in file called "samplenames.txt" (one name per line - make sure there are no blank lines at the end of the file) in the folder with a copy of your fasta files. Run the following (modified from step 2 of : https://github.com/laninsky/reference_aligning_to_established_loci). This will also give us a list of the loci which we can use in step 4 (excluding loci for samples which are made up of nothing but missing data).
```
nosamples=`wc -l samplenames.txt | awk '{print $1}'`;

for i in `seq 1 $nosamples`;
do refname=`tail -n+$i samplenames.txt | head -n1`;
mkdir $refname;
for j in `ls *.fasta`;
do temp=`grep -A1 $refname $j`
echo $temp | awk '{print $2}' > temp
if grep --quiet [^\?] temp;
then echo ">"$j >> $refname/reference.fa;
cat temp >> $refname/reference.fa;
fi;
done;
done;
ls *.fasta > fasta_names
rm temp
```
The *.fasta files for each locus (but not reference.fa files) can then be removed.

3) Make sure the F and R reads are named the same within the READ1 and READ2 files within the cleaned reads folder - you only have to do this step if the reads have been prefaced with a 1: in the READ1 folder and a 2: in the READ2 folder (running this from inside the cleaned reads folder):
```
for i in `\ls`; do gunzip $i/split-adapter-quality-trimmed/*READ2*; done;
for i in `\ls`; do sed -i 's/ 2:/ 1:/g' $i/split-adapter-quality-trimmed/*READ2*; done;
for i in `\ls`; do gzip $i/split-adapter-quality-trimmed/*READ2*; done;
```

4) Back in the main folder where your locus.fasta files used to be, and where you now have your subfolders consisting of each of your samples (each containing a reference.fa file): run the phase_everyone.sh script in this folder (make sure coverage.R from the main folder and modref.R from this folder are also in here, and phasing_settings etc has been set up as for the main instructions - note the comment about ${name} in the main instructions).
```
bash phasing_everyone.sh
```

5) Run summarize_coverage.R from the main repository to obtain a locus-by-locus summary of coverage and sample size. The output file will give the locus name in the first column, how many samples that locus was present in, and then the average coverage across these samples. Rejoin the main set of scripts, resuming at step 5, at:
https://github.com/laninsky/reference_aligning_to_established_loci
```
sed -i 's/-//g' fasta_names
mv samplenames.txt samples.txt
Rscript summarize_coverage.R
```

Follow any instructions for UCE-related points in the main scripts.

### Programs/packages necessary for the pipeline:
```
stringr
bwa
samtools
java
picard
gatk
R
```

Along with the programs above, to cite this pipeline:
```
Alexander, A. 2018. phase_everyone v0.0. Avaiable from: https://github.com/laninsky/reference_aligning_to_established_loci/edit/master/phase_everyone

Baca, S.M., Alexander, A., Gustafson, G.T. and Short, A.E., 2017. Ultraconserved elements show utility in phylogenetic inference of A dephaga (C oleoptera) and suggest paraphyly of ‘Hydradephaga’. Systematic Entomology, 42(4), pp.786-795.
```

### Version history
0.0: Version used in Baca et al. (2017)

