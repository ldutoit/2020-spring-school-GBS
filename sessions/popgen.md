# Assembly with a reference genome

**Developed by:** Julian Catchen, Nicolas Rochette

**Adapted by:** Ludovic Dutoit

### Goals
Learn to run basic population genetics analyses
Compare a dataset with and without a reference genome


### Introduction

...

Each of you does one of the two...

Run Structure

Then all into R into the script [run_popgen.R](run_popgen.R) USE BULLY TUTO
 
The code is pre-made.

### Structure
Our goal now is to export a subset of loci for analysis in Structure, which analyzes the distribution of multi-locus genotypes within and among populations in a Bayesian framework to make predictions about the most probable population of origin for each individual. The assignment of each individual to a population is quantified in terms of Bayesian posterior probabilities, and visualized via a plot of posterior probabilities for each individual and population.

A key user defined parameter is the hypothesized number of populations of origin which is represented by K. Sometimes the value of K is clear from from the biology, but more often a range of potential K-values must be explored and evaluated using a variety of likelihood based approaches to decide upon the ultimate K. In the interest of time we won’t be exploring different values of K here, but this will be a key step for your own data sets. In addition, at the moment Structure will take a long time to run on the number of loci generated in a typical RAD data set because of the MCMC algorithms involved in the Bayesian computations. We therefore want to randomly choose a random subset of loci that are well represented in our three populations. Nonetheless, this random subset contains more than enough information to define population structure:

The final stage of the denovo_map.pl pipeline is to run the populations program. We want to execute just populations, rather than the full denovo_map.pl pipeline, to specify filters that will give us only the most well represented loci. Populations is a very useful pice of software both for filtering and for outputting population genetics. It can work with non-stacks generated data too.

• Since we won't be able to use all loci for our quick downstream analysis today, we will run populations again, specifying that loci must be present in at least 80% of individuals in all three populations to cut down on the total number of loci. You will have to tell populations where to find the output of the Stacks pipeline (this should be in the stacks output directory).

• Make sure to output a structure file! Output a .vcf file as well. We will be combing back that .vcf file later today later.

• One final detail: Structure assumes that each SNP locus is independent, so we don’t want to output multiple SNPs from the same RAD locus, since they are not independent but are in linkage blocks within each RAD tag. We can achieve this behavior by specifying the --write_single_snp parameter to populations.

Now we want to select 1,000 loci randomly from the results and save these loci into a file. We can easily do this using the shell given a list of catalog IDs output in the previous step. The populations.sumstats.tsv file gives a list of all polymorphic loci. Use the cat, grep, cut, sort, uniq, shuf, and head commands to generate a list of 1000 random loci. Save this list of loci as whitelist.txt, that we can feed back into populations. This operation can be done in a single shell command. That sounds challenging, but the instructions below should help you create one command with several pipes to create that whitelist.txt. Create that command step by step:

• First, use cat to concatenante stacks/populations.sumstats.tsv.

• Then, use grep with -v to exclude all headers (i.e. excluding "#")

• Select the first column with cut -f 1

• Then use sort before using uniq. uniq will collapse succesive identical lines into single lines, so that you have one line per locus. Lucky us, those two commands don't require any arguments.

• Now try adding shuf which will mix all these lines all over.

• Select the first one thousand lines with head and put it all into whitelist.txt.

• You got this! If you are new to bash, I am sure that seemed impossible yesterday, so take a minute to congratulate yourself on the progress made even if you required a little help!

Now execute populations again, this time feeding back in the whitelist you just generated. This will cause populations to only process the loci in the whitelist.txt. Specify that a Structure output file be included this time and again insist that only a single SNP is output from each RAD locus. Finally, you will need to again specify the population map that you generated above to populations so that this information is passed into the Structure output file. You might also want to output a vcf file: this format is handy for all kind of population genetics applications.

• We've run commands to generate the structure file several times, but how many structure files are there in the stacks directory? If you wanted to save several different vcf and structure files generated using different populations options, what would you have to do?

Create a new directory called structure within the denovo folder and copy the Structure output file that Stacks generated to this directory. cd into your new structure directory.

• Edit the Structure output file to remove the comment line (first line in the file, starts with “#”).
• The parameters to run Structure (with value of K=3) have already been prepared, you can find them here: /nesi/nobackup/nesi02659/source_data/denovo/mainparams and /nesi/nobackup/nesi02659/source_data/denovo/extraparams
• Copy them into your structure directory as well.

So far, when we've gone to run programs, we've been able to use module spider to figure out the program name, and then module load program_name to get access to the program and run it. However, structure is not an available module on Mahuika. Instead, we've done a local installation of the progam into /nesi/nobackup/nesi02659/source_data/structure. Run Structure:

 /nesi/nobackup/nesi02659/source_data/structure
The program should ive you some information as it runs. If the program immediately finishes, something has gone wrong! Do a less on populations.structure.console. Do you see WARNING! Probable error in the input file.? In our mainparams file it says that we have 1000 loci, but due to filters, it is possible that the populations module of Stacks actually output less than the 1000 loci we requested in whitelist.txt. In the output of populations.log in your stacks directory, how many variant sites remained after filtering? This is the number of loci actually contained in your structure file. You will need to adjust the number of loci in the mainparams file to match this exact Stacks output.

You will need to download Structure with the graphical front end and use scp to download the populations.structure.out_f file from the cluster. You can then load this file into the graphical interface for Structure on your local computer. Select the File menu and then Load structure results to load the Structure output. Choose the Barplot menu and then Show.

• Are the three Oregon threespine stickleback populations related to one another? How can you tell?

Congrats, you just finished our tutorial for denovo RAD-Seq. If you have plenty of time, you could try different parameters for the populations module of Stacks or familiarise yourself with the idea of a reference based approach: ref_map.pl. You could also play with your own data. Finally, you could use the vcf you generated into the stacks folder in the R package adegenet to create other visualisations if you are comfortable with R
