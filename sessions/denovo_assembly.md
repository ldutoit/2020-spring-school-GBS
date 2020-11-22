# De-novo assembly without a reference genome

**Developed by:** Julian Catchen, Nicolas Rochette

**Adapted by:** Ludovic Dutoit

### Goals
  
  - Learn to call SNPs without a reference genome and to optimise the procedure
  - Learn to run SLURM scripts

### Introduction

In this second exercise we will  be working on  threespine stickleback data sampled from throughout Oregon, on the west coast of the United States. These data consist of 30 samples in three populations: a coastal marine population, a coastal freshwater, and an inland river population. These stickleback can be found in a
number of habitats from costal
marine and freshwater habitats, to
inland river habitats, to high
mountain lakes in the interior of
Oregon. We want to understand how
these populations relate to one
another and in this exercise, you will
examine three of these populations:
a coastal marine population, a costal
freshwater, and an inland river
population. For more information on
the background of this study, see [Catchen et al, 2013](https://onlinelibrary.wiley.com/doi/10.1111/mec.12330).

Without access to a reference genome, we
want to assemble the RAD loci and
examine population structure. However, before we can do that, we want to explore
the *de novo* parameter space for the in order to be confident that we are assembling our data
in an optimal way. Stack (i.e. locus) formation is controlled by three main
parameters: 

-m :  the minimum amount of reads to create a locus)

**-M : the number of mismatches allowed between alleles of the same locus**

-n : The number of mismatches between between loci between individuals

If that does not make sense or you would like to know more, have a quick read of [this explanation from the manual](http://catchenlab.life.illinois.edu/stacks/param_tut.php).

Here, we will optimize the parameter M using the collaborative power of our wonderful team. After
We will be using the guidelines of parameter optimization as outlined in [Paris
et al. (2017)](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12775) to assess which parameters value for M recovers the highest number of new polymorphic loci found across 80% of the individuals (r80 loci).

This approach is described more in [Paris et al. (2017)](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12775)

*"After putative alleles are formed, stacks performs a search to match alleles together into putative loci. This search is governed by the M parameter, which controls for the maximum number of mismatches allowed between putative alleles [...] Correctly setting **M** requires a balance – set it too low and alleles from the same locus will not collapse, set it too high and paralogous or repetitive loci will incorrectly merge together."*

As a giant research team,  we will run the *denovo* pipeline with different parameters. The results from the different parameters will be shared using [The Google sheet](https://docs.google.com/spreadsheets/d/13qm_fFZ4yoegZ6Gyc_-wobHFb7HZxp27mrAHGPmnjRU/edit#gid=0). We'll be able to use the best dataset downstream for population genetics analyses and comparison with a pipeline that utilises a reference genome.

## Build your denovo_map.pl command

1. We will assemble loci and use the collaborative power of this classroom to determine the best parameters. 

2. In your ```GBS``` workspace, create a directory called ```output_denovo``` to contain
the assembled data for this exercise.

3. To avoid duplicating the raw data for each of us, we will use a link to the source data. This effectively creates a shortcut to another path without copying all the files. 
`ln -s /path/you/want/to/link` will create a shortcut to a given path right where you are! The raw data is in 
```/scale_wlg_persistent/filesets/project/nesi02659/obss_2020/resources/day3/oregon_stickleback/``` Using the above example, create a link to this folder right here!


4. Run Stacks’ denovo_map.pl pipeline program according to the following set of instructions. Following those instructions you will bit by bit create the complete `denovo_map.pl` command:
    
    • Make sure you load the ```Stacks``` module (you can check if you already loaded it using `module list`)
    
    • Get back in the ```GBS``` folder if you wandered away.
    
    • Information on denovo_map.pl and its parameters can be found [online](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php). You will use this below to build your command.
    
    • We want Stacks to understand which individuals in our study belong to which population. Stacks use a so-called population map. The file contains sample names as well as populations. The file should be formatted in 2 columns like [this](http://catchenlab.life.illinois.edu/stacks/manual/#popmap). All 30 samples are at the file path below. Copy it in the folder `GBS you` should currently be in.
    
    ```/scale_wlg_persistent/filesets/project/nesi02659/obss_2020/resources/day3/denovo_popmap.txt```

    • Make sure you specify this population map to the denovo_map.pl command.
    
    • There are three important parameters that must be specified to denovo_map.pl, the
        minimum stack/locus depth (`m`), the distance allowed between stacks/loci (`M`), and the distance allowed
        between catalog loci (`n`) **that should be M+2**. Use the values we determined for these parameters in the
        previous exercise, but do not restrict the loci to just those found in 80% like we did in the opt runs.
        Choose which values of M (M<8) you want to run, not overlapping with other people parameters and insert them in [The google sheet](https://docs.google.com/spreadsheets/d/13qm_fFZ4yoegZ6Gyc_-wobHFb7HZxp27mrAHGPmnjRU/edit#gid=0). You can vary M and n(which should be M+2) as well as `-r` anywhere between 50 and 100%.
    
    • You must set the stacks directory as the output, and use 6 threads (6 CPUs so your analysis finishes faster than 1!).
        
    • Specify the path to the directory containing your sample files (*hint* use your *oregon_stickleback/* link here!.       The denovo_map.pl program will read the sample names out of the population map, and
        look for them in the samples directory you specify.
       
    • Your command should be ready, try to execute the Stacks pipeline. 

    • Is it starting alright?  Good, now  **Use `control + c` to stop your command**

5.Running the commands directly on the screen is not common practice. You now are on ga-vl01 which is a reserved amount of resources for this workshop and this allows us to run pur command directly. On a day to day basis, you would be evolving on the login node (i.e. The place you reach when you login). All the resources are tucked away from the login node. You generally run your commands as jobs that are sent to this resources, not on the login node itself. We will use this denovo_map.pl command as a perfect example to run our first job.

  • copy an example jobfile into this directory. The example is at :                  ```/scale_wlg_persistent/filesets/project/nesi02659/obss_2020/resources/day3/denovojob.sh```

  • Open it with a text editor, have a look at what is there. The first bit are parameters for the job that starts with . system informing on who you are, which type of resources you need and for how long.

  • There are a number of <...> followed by a comment starting by `#` that tells you what should be there, fill in your job script.

  • Once you are done, save it. Then run it using:

    sbatch denovojob.sh


 • You can check what is the status of your job using

    squeue -u <yourusername>
 
 • We used a few basic options of sbatch, time, memory and ... . In reality, there are many many more options, have a quick look at sbatch --help out of interest. NeSI also has its own handy guide on how to submit a job [here](https://support.nesi.org.nz/hc/en-gb/articles/360000684396-Submitting-your-first-job).

• Once `squeue` is empty, your job ran and what would have printed to your screen is into denovo.out. That should take about 30mn to run, so in the meantime, sit back and relax, we are getting to lunch!


### Analysing the data from our collaborative optimisation

5. Examine the Stacks log and output files when execution is complete. It should be in `output_denovo`
    
    • After processing all the individual samples through ustacks and before creating the catalog with cstacks, denovo_map.pl   will print a [table containing the depth of coverage](http://catchenlab.life.illinois.edu/stacks/manual/#cov) of  each sample. Find this table in the log, what were the depths of coverage? 
    
    • Examine the output of the populations program in the file XXXXX.log inside your `output_denovo` folder. (*hint*: use the `less` command).
    
    • How many loci were identified?

    • How many SNPs were identified?

     Enter that information in the collaborative [Google Sheet](https://docs.google.com/spreadsheets/d/13qm_fFZ4yoegZ6Gyc_-wobHFb7HZxp27mrAHGPmnjRU/edit?usp=sharing)
    
    • How many were filtered and for what reasons?
    
    • Familiarize yourself with the population genetics statistics produced by the populations component of stacks populations.sumstats_summary.tsv
    
    • What is the mean value of nucleotide diversity (π) and FIS for each of the three
        populations? [HINT: The less -S command may help you view these files easily by avoiding the wrapping]


Congratulations, you obtained variants.
