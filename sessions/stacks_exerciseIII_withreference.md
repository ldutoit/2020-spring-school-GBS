# Assembly with a reference genome

**Developed by:** Julian Catchen, Nicolas Rochette

**Adapted by:** Ludovic Dutoit

### Goals
  
  - Learn to call SNPs with a reference genome

### Introduction

Obtaining an assembly without a reference genome is easy and possible. However, having some reference genome allow us to avoid several issues. We do not have to make assumptions selcting the -M parameters, we reduce the risk of collapsing as once different loci or to separate single loci into several "erroneous loci". Several studies  have demonstrated that having some kind of reference genome is the single best step to improve a GBS SNP calling (... , ... , ...) . In this exercise, we will use the publicly available Stickleback genome to extract variants on our 30 samples datasets. We will conclude today by comparing popuation structure between the two dataset with and without a reference genome.

The stacks pipeline for samples with a reference genome is [ref_map.pl](https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php), it skips the creation of loci and the catalog steps as determining whether two reads are belonging to the same position in the genome is not dependent on assumptions derived from their genetic distances (-M and -n of denovo_map.pl). They belong to the same stack/locus if they map to the same location of the reference genome. 


The 30 samples have been mapped for you using BWA to save you time and effort. If you're curious, the mapping script is [here](MAPPING SCRIPT). 
THE FILES HAVE BEEN MAPPED WITH YOU BUT YOU CAN SEE AN EXAMPLE MAPPING SCRIPT [HERE](LINKED TO MAPPING).

### Run Stacks with a reference genome

In the folder `GBS` create the output folder for this analysis, `output_refmap`.

we'll use the same samples as before. The mapped version is in:


```/scale_wlg_persistent/filesets/project/nesi02659/obss_2020/resources/day3/oregon_stickleback_mapped/```

create a link (*hint* ln -s /.../../ .)

refmap_map.pl has less options since the mapping take care of many of the steps from denovo_map.pl such as the creation of loci for each individuals before a comparison of all loci across all individuals. Use the help [HERE]() to build your refmap command.

• Like in the previous exercise, ask for 6 threads 

• Specify the path to the output folder output_refmap
Is your command ready? Run it briefly to check that it starts properly, once it does, **stop it using ctrl + c**


Time to put that into a job script. You can use the job-script from the last exercise. We will simply make a copy of it under a different name.

   ```cp denovojob.sh refmapjob.sh```

now open Adjust the time to ... , the job name to. ... and the memory.

Run the job. That should take about X minues.

Colelct the main results.

In the last exercise, we obtained Y ... and ..

How does this run compare? 

As a thought exercise Do you have any idea on how to check whether those are the same loci or not?

..

Well done, you now know how to call SNPs with or without a reference genome. It is worth nothing that even a poor-quality reference genome will improve the quality of your SNP calling. 




