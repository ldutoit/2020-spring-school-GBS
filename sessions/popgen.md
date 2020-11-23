# Population genetics analyses

**Developed by:** Julian Catchen, Nicolas Rochette

**Adapted by:** Ludovic Dutoit

### Goals
Learn to run basic population genetics analyses



### Structure

Our goal for now is to export a subset of loci for analysis in [Structure](https://web.stanford.edu/group/pritchardlab/structure.html), which analyzes the distribution of multi-locus genotypes within and among populations in a Bayesian framework to make predictions about the most probable population of origin for each individual. The assignment of each individual to a population is quantified in terms of Bayesian posterior probabilities, and visualized via a plot of posterior probabilities for each individual and population.

A key user defined parameter is the hypothesized number of populations of origin which is represented by K. Sometimes the value of K is clear from from the biology, but more often a range of potential K-values must be explored and evaluated using a variety of likelihood based approaches to decide upon the ultimate K. In the interest of time we won’t be exploring different values of K here, but this will be a key step for your own data sets. In addition, at the moment Structure will take a long time to run on the number of loci generated in a typical RAD data set because of the MCMC algorithms involved in the Bayesian computations. We therefore want to randomly choose a random subset of loci that are well represented in our three populations. Nonetheless, this random subset contains more than enough information to define population structure:

The final stage of the Stacks pipeline is to run the `populations` program. Now, we want to execute just populations, rather than the full Stacks pipeline, to specify filters that will give us only the most well represented loci. Populations is a very useful piece of software both for filtering and for outputting population genetics. It can work with non-stacks generated data too.

The help of populations is [here](https://catchenlab.life.illinois.edu/stacks/comp/populations.php)

 • Since we won't be able to use all loci for our quick downstream analysis today, we will run populations again, specifying that loci must be present in at least 80% of individuals in all three populations to cut down on the total number of loci in the directory).

 • Make sure to output a `structure file` and a `.vcf` file as well. We might be coming back to that `.vcf` file later today.
 
 • Specify `output_refmap` as the path to the directory containing the Stacks files 
  
 • Also specify `output_refmap` as the output folder

This command is now ready, assuming you are into GBS, run it! It is a quickie, so no need to put it in a job.
 
How many SNPs do you have left? We want to select 1,000 loci randomly from the results to run a relatively quick Structure analysis. We will save these loci into a file. We can easily do this using the shell given a list of catalog IDs output in the previous step. The populations.sumstats.tsv file gives a list of all polymorphic loci. 

With the help of the instructions below, use the cat, grep, cut, sort, uniq, shuf, and head commands to generate a list of 1000 random loci. Save this list of loci as whitelist.txt, that we can feed back into populations. This operation can be done in a single shell command. That sounds challenging, but the instructions below should help you create one command with several pipes to create that `whitelist.txt` file. The idea of a pipe is to connect commands using `|`. `command 1 | command 2` outputs the content of command 1 into command 2 instead of outputting it to the screen.

Create that command step by step:

• First, use `cat` to concatenante `output_refmap/populations.sumstats.tsv`.

• Then, add a `|` in your command and use the command `grep` with `-v` to exclude all headers (i.e. `-v` stands for exclude, we want to exclude all the lines with "#")

• Add another `|` and select the first column with `cut -f 1`

• Then using two more pipes, use `sort` before using `uniq`. `uniq` will collapse succesive identical lines into single lines, so that you have one line per locus. Lucky us, those two commands don't require any arguments.

• Now try adding `shuf` which will mix all these lines all over.

• Select the first one thousand lines with `head -n 1000` and put it all into whitelist.txt. (*hint*: use ">" to redirect into a file).

• Do not worry about the shuf: `write error: Broken pipe`, it is simply because head stops the command before the end.
 

You got this! If you are new to bash, I am sure that seemed impossible on monday, so take a minute to congratulate yourself on the progress made even if you required a little help!

Now we will execute `populations` again, this time feeding back in the whitelist you just generated, check out the help of populations to see how to use a white list. This will cause populations to only process the loci in the `whitelist.txt`. 

• The [help of Populations](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) will tell you how to pass a white list.

• Specify that a Structure output file be included this time.
 
• Specify `--write-single-snp` that should generate one snp per locus as we need genetically unlinked SNPs (i.e. statistically)

• Finally, you will need to again specify the population map so that this information is passed into the Structure output file. **Careful** This time we will specify the population map with the proper population information. That file is at `/nesi/project/nesi02659/obss_2020/resources/day3/complete_popmap.txt`. You can copy this file here, link it, or simply specify the full path, your call! Al of these 3 solutions should work.

• We've run commands to generate the structure file two times now, but how many structure files are there in the stacks directory? If you wanted to save several different vcf and structure files generated using different populations options, what would you have to do?

Create a new directory called structure within the `GBS` folder and copy the Structure output file that Stacks generated to this directory. `cd` into your new structure directory.

• Edit the Structure output file to remove the comment line (i.e. first line in the file, starts with “#”).

• The parameters to run Structure (with value of K=3) have already been prepared, you can find them here: `/nesi/project/nesi02659/obss_2020/resources/day3/mainparams` and `/nesi/project/nesi02659/obss_2020/resources/day3/extraparams`. Copy them into your structure directory as well.

• So far, when we've gone to run programs, we've been able to use `module spider` to figure out the program name, and then module load program_name to get access to the program and run it. Do it one more time for `structure`

•  run structure by simply typing `structure` 

 
Do you see `WARNING! Probable error in the input file.?` In our mainparams file it says that we have 1000 loci, but due to filters, it is possible that the populations module of Stacks actually output less than the 1000 loci we requested in whitelist.txt. In the output of populations.log in your `output_refmap` directory, how many variant sites remained after filtering? This is the number of loci actually contained in your structure file. You will need to adjust the number of loci in the mainparams file to match this exact Stacks output.

!!!You will need to download [Structure] with the graphical front end and use scp to download the populations.structure.out_f file from the cluster. You can then load this file into the graphical interface for Structure on your local computer. Select the File menu and then Load structure results to load the Structure output. Choose the Barplot menu and then Show.

• Are the three Oregon threespine stickleback populations related to one another? How can you tell?

Congrats, you just finished our tutorial for denovo RAD-Seq. Are you already done? It looks like it, here are a few things you coyld do:

• You could spend a bit of time going through what you have done today and make sure you understand the set of steps we did in those few exercises.

• You could try running this set of analyses on a denovo dataset. You can use the one you generated yourself or a -M 2 dataset that is in `/nesi/project/nesi02659/obss_2020/resources/day3/output_denovo_M2`.

• You could have a look at this [tutorial](populationstructure_tuto/populationstructure_tuto.md). It is a small tutorial I wrote once that go over a different set of population genetics analyses in R. You could even try reproducing it using the vcf you generated above. The vcf for that particular set of analysis can be downloaded [here on the original repos](https://github.com/ldutoit/bully_gbs/blob/master/output_files/populations.snps.vcf), should you want to download it.

[Back to homepage](../index.md)
