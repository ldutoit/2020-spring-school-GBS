# De-novo assembly of RAD tags without a genome

**Developed by:** Julian Catchen, Nicolas Rochette

**Adapted by:** Ludovic Dutoit


## Exercise 2: Call variants in the absence of a reference genome

### Goals
  
  - Learn to call SNPs without a reference genome and to optimise the procedure
  - Learn to run SLURM scripts

### Introduction

1. In this second exercise we will  be working on  threespine stickleback data sampled from throughout Oregon, on the west coast of the United States. These data consist of three populations: a coastal marine population, a coastal freshwater, and an inland river population. These stickleback can be found in a
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
the *de novo* parameter space in order to be confident that we are assembling our data
in an optimal way. Stack formation is controlled by three main
parameters: m (the minimum read depth); M (the number of mismatches between
alleles) and n (the number of mismatches between loci in the catalog). Here, we will
optimize M for the stickleback data using a subset of the full dataset provided. After
this, we can use the optimal value we have found for M in the *de novo* exercise below. We will be using the guidelines of parameter optimization as outlined in [Paris
et al. (2017)](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12775) to assess which value for M recovers the highest number of new polymorphic loci found across 80% of the individuals (r80 loci).

This approach is described more in [Paris et al. (2017)](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12775)

*"After putative alleles are formed, stacks performs a search to match alleles together into putative loci. This search is governed by the M parameter, which controls for the maximum number of mismatches allowed between putative alleles [...] Correctly setting **M** requires a balance – set it too low and alleles from the same locus will not collapse, set it too high and paralogous or repetitive loci will incorrectly merge together."*

As a giant research team,  we will run the *denovo* pipeline with different parameters. The results from the different parameters will be shared using [Th Google sheet]. We'll be able to use the best dataset downstream for population genetics analyses and comparison with a pipeline that utilises a reference genome.


2. We will assemble loci and use the collaborative power of this classroom to determine the best parameters.  
3. In your ```GBS``` workspace, create a directory called ```output_denovo``` to contain
the assembled data for this exercise.

4. Run Stacks’ denovo_map.pl pipeline program according to the following set of instructions:
    
    • Load the ```Stacks``` module
    
    • Get back in the ```GBS/denovo``` folder.
    
    • As for the previous exercise, information on denovo_map.pl and its parameters can be found [online](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php)
    
    •We want Stacks to understand which individuals in our study belong to which population. To specify this, create a file ```complete_popmap.txt``` in the denovo directory called popmap, using
        an editor. The file should be formatted in 2 columns like [this](http://catchenlab.life.illinois.edu/stacks/manual/#popmap). Include all 30 samples in this file and specify which individuals            belong to which populations. You must supply the population map to denovo_map.pl when you execute it. You could for example use ```ls -1 *fa.gz``` to see all the samples in a list before adding the populations. Add the populations as simple integers (i.e. 1, 2 and 3) t
    
    • There are three important parameters that must be specified to denovo_map.pl, the
        minimum stack depth (`m`), the distance allowed between stacks (`M`), and the distance allowed
        between catalog loci (`n`). Use the values we determined for these parameters in the
        previous exercise, but do not restrict the loci to just those found in 80% like we did in the opt runs.
    
    Choose which parameters you want to run, GOOGLE SHEET
    • Also, you must set the stacks directory as the output, and use 6 threads (6 CPUs so your analysis finishes faster than 1!).
        
    • Finally, specify the path to the directory containing your sample files. The
        denovo_map.pl program will read the sample names out of the population map, and
        look for them in the samples directory you specify.
    • To optimize for r80 loci you will need to tell denovo_map.pl to use the '''-r''' parameter to filter for loci in 80% of the
        samples) program. We will keep ```m``` at 3. Initially, we will set M to 4. We will also follow the general rule of ```M = n``` and we will tell [denovo_map.pl](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) to output to the M4 folder.
    
        Put the command here : ....

    • Execute the Stacks pipeline. 

    let's look at it, is it starting alright? ...

    hat should take approximately 30min, ideal time for a break! WAIT NO, DON'T DO THAT THIS IS HORRIBLE PRACTICE

    Introduce SLURM

    4.Running the commands directly on the screen is not common practice. You now are on ga-vl01 which is a reserved amount of resources for this workshop and this iallows us to run pur command directly. On a day to day basis, you would be evolving on the login node (i.e. The place you reach when you login). All the resources are tucked away from the login node. You generally run your commands as jobs that are sent to this resources, not on the login node itself. We will use this process_radtags command as a perfect example to run our first job.

• copy an example jobfile into this directory. The example is at : '''/nesi/nobackup/nesi02659/source_data/example_job.sh'''

Describe the basic parameters in there, link it to the actual help of batch.


• Open it with nano, have a look at what is there. As you can see, the first bit is parameters for the job . system informing on who you are, which type of ressources you need and for how long

• Once you are done, save it and run it using

    sbatch examplejob.sh

• You can check what is the status of your job using

    squeue -u <yourusername>
• Once this place is empty, your job ran and what would have printed to your screen is into prcoessrads.out. You should also have recieved an email!
    
5. Examine the Stacks log and output files when execution is complete.
    
    • After processing all the individual samples through ustacks and before creating the catalog with cstacks, denovo_map.pl    will print a [table containing the depth of coverage](http://catchenlab.life.illinois.edu/stacks/manual/#cov) of  each sample. Find this table in the log, what were the depths of coverage? 
    
    • Examine the output of the populations program in the log.
    
    • How many loci were identified?

    • How many SNPs were identified?

     Enter that information in the collaborative [Google Sheet](https://docs.google.com/spreadsheets/d/13qm_fFZ4yoegZ6Gyc_-wobHFb7HZxp27mrAHGPmnjRU/edit?usp=sharing)
    
    • How many were filtered and for what reasons?
    
    • Familiarize yourself with the population genetics statistics produced by the populations component of stacks populations.sumstats_summary.tsv
    
    • What is the mean value of nucleotide diversity (π) and FIS for each of the three
        populations? [HINT: The less -S command may help you view these files easily by avoiding the wrapping]


FOLLOWING X, what are the best parameters? 


##  Running Stacks denovo

1. Go to your ```XXX``` workspace, create a directory called ```denovo``` to contain all the
data for this exercise. *Inside* that directory, create two additional directories:
```samples```, and ```opt```. 

   To save time, we have already cleaned and demultiplexed this
   data and will start from the cleaned samples stage. Inside the opt directory, create four
   additional directories: ```M4```, ```M5```, ```M6```   and ```M7```. 

   • Go to ```samples``` directory. Copy the dataset below in the ```samples``` directory: ```/nesi/nobackup/nesi02659/source_data/denovo/oregon_stickleback.tar```
    
   • Extract it. The unarchived dataset contains 30 stickleback
     samples (Catchen, et al., [2013](https://onlinelibrary.wiley.com/doi/10.1111/mec.12330)), but we will use only 3 of them (`cs_1335.01`,  `pcr_1211.04`, `stl_1274.33`) in this first part of the exercise as we will run denovo_map.pl just a few times for optimisation. 
     
2. We will run the Stacks’ ```denovo_map.pl``` pipeline program, each time changing the value for
```M```. `denovo_map.pl` will run ustacks, cstacks, and sstacks on the individuals in our
study as well as the ```population``` program. The set of instructions below should help you build your command.
    
   • Get back into the ```opt``` folder. (*Hint*:  Use```pwd``` if you don't know where you are anymore.
    

    
   • Information on denovo_map.pl and its parameters can be found [online]( http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php)
    
   • We want Stacks to only use the 3 individuals in our parameter optimization (cs_1335.01,  pcr_1211.04, stl_1274.33).
   
   • To specify this, create a file in the opt directory called ```opt_popmap.txt```, using an editor.
        The file should be formatted like [indicated in the manual](http://catchenlab.life.illinois.edu/stacks/manual/#popmap).        
         Note: do not include the extension ```.fa.gz``` in the sample name.
    
   • Include samples in this file and **specify that all individuals belong to one
        single population** (e.g. give them all the same population code). You will need to supply this        ```opt_popmap.txt```population map to [denovo_map.pl](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) when you
        execute it for each parameter run.
    
   • To optimize for r80 loci you will need to tell denovo_map.pl to use the '''-r''' parameter to filter for loci in 80% of the
        samples) program. We will keep ```m``` at 3. Initially, we will set M to 4. We will also follow the general rule of ```M = n``` and we will tell [denovo_map.pl](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) to output to the M4 folder.
        
   • With this information, you should be able to launch the M4 run now. It will take a couple of minutes.
       
3. Once done, you should now see the denovo_map.pl output files in the directory M4.

    • After M4 is completed, do the same process for M = 5 through M = 7 (this should take around 10min total).    
    
    • While you are running M5 through m7, open a *new command window*, login to Mahuika and re-access       the reserved machine ```ssh -Y ga-vl01```. Then go to your ```working/denovo``` folder so that you can keep working while M5 through M7 are running.
    
    • To see how many r80 loci were assembled for each parameter run you will want to start looking at the                        ```populations.hapstats.tsv``` file using the [Stacks manual](http://catchenlab.life.illinois.edu/stacks/manual/#files)
    to inform you on the data
    contained in each of the columns. 

    • What is the number of the first locus assembled for M4?

    • What is the number of the last locus assembled for M4?
    
    • How many loci in total? (*Hint:* count the lines)
        
    • Using this technique, how many loci were assembled for M5 to M7 once they finish running?
    
    • Which iteration of M provided the highest number of r80 loci?   

You should now be able to choose your optimised parameters according to the Paris et al. (2017) method! ("select those values which maximize the number of polymorphic loci found in 80% of the individuals in your study")


