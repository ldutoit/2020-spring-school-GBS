# Exercise I. Data preparation    

**Developed by:** Julian Catchen, Nicolas Rochette

**Adapted by:** Ludovic Dutoit


## Goal

The goal of this first exercise is to take you from raw reads data to individual sample files.

• This first data set comprises just a small proportion of a lane of single-end
standard RAD data.

• You will use this first dataset to become familiar with the structure of RAD sequences, as well as to become proficient with the pre-processing (i.e. cleaning and de-multiplexing) of data before alignment or assembly. 

## Part 1: Single-end reads

The first step in the analysis of all short-read sequencing data, including RAD-seq
data, is removing low quality sequences and separating out reads from different
samples that were individually barcoded. This **‘de-multiplexing’** serves to associate
reads with the different individuals or population samples from which they were
derived.

1. Let's organise our space, get comfortable moving around and copy our data :
    
    • Log into Jupyter at [https://jupyter.nesi.org.nz/hub/login](https://jupyter.nesi.org.nz/hub/login). Not sure how to do it? Just follow the instructions [jupyter_hub.md](jupyter_hub.md)
       
    
    • For each exercise, you will set up a directory structure on the remote server that will hold your data and the different steps of your analysis. We will start by making the directory ```GBS``` in your working space, so let's `cd` (change directory) to your working location:
       
       
       ```cd /scale_wlg_persistent/filesets/project/nesi02659/obss_2020/users/<yourusername>/```
       
       OR from the launch of Jupyter:
     
       ```cd users/<yourusername>/```
       *NOTE* If you get lost an any time today, you can always cd in your home following this upper link.
       
    • Once there, create the directory `GBS` and then change directory into `GBS`:
      
        ``` 
        mkdir GBS
        cd GBS
        ```
     • The exercise is hands-on, the instructions are there to guide you through the process but you will have to come up with the commands yourself. Fear not tho, the instructions in the text are there to help you andthe room is full of friendly faces here to help you get through it. 
      
     •   We will create more subdirectories to hold our analyses. Be careful that you are reading and writing files to the appropriate directories within your hierarchy. You’ll be making many directories, so stay organized!
    
    • Each step of your analysis goes into the hierarchy of the workspace, and each step of  
        the analysis takes its input from one directory and places it into another director. We will name the                   directories in a way that correspond to each stage and that allow us to remember where they are. A well
        organized workspace makes analyses easier and prevents data from being overwritten.
    
    • First let's make a few directories. In ```GBS```, create a directory called ```dataprep``` to contain all the data  for this exercise. Inside that directory we will create two additional directories: ```lane1``` and ```samples```. 
    
    • As a check that we've set up our workspace correctly, go back to your ```<username>``` directory (*hint*: `cd ..`) and use the `ls -R` (the `ls` command with the recursive flag). It should show you the following:
    
    ```
    GBS/:
    dataprep

    GBS/dataprep:
    lane1  samples

    GBS/dataprep/lane1:

    GBS/dataprep/samples:
    
    ```
    
    • Copy the data set 1 (DS1) to your ```lane1``` directory. The data set is in the file
       `/scale_wlg_persistent/filesets/project/nesi02659/obss_2020/resources/day3/lane1.tar` 
       (*hint*: `cp /path/to/what/you/want/to/copy /destination/you/want/it/to/go`)          
    
    • `cd` to your ```lane1``` folder to extract/unzip the content of this ```tar``` archive. this is a compressed folder. We realise that we have not told you how to do so! But a quick look to a friendly search engine will show you how easy it is to find this kind of information on basic bash commands (your instructors *still* spend a lot of time doing this themselves!).    
    *hint* : you might try searching for "How to extract a tar archive"

2. Have a look at what is there now. These gz-compressed fastq files have millions of reads in them, too many for you to examine in a spreadsheet or word processor. However, we can examine the contents of the set of files in the terminal
(the ```less``` command may be of use).
    
3. Let's have a closer look at this data. Over the last couple of days, you learnt to run FastQC to evaluate the quality of the data. We'll save you the trouble of running it here ... In reality, sequencing platform often provide you with the quality reports too. This link [](...) is the fastqc report for file ... . Download it and have a quick look at it.

Should you want to run this yourself:

```
module load FastQC
fastqc *gz
```

Let's look at this FastQC report together:

    • What is this weird thing in the base-pair content from base 7 to 12-13?

    • You probably noticed that not all of the data is high quality. In general, you will want
      to remove the lowest quality sequences from your data set before you proceed.
      However, the stringency of the filtering will depend on the final application. In
      general, higher stringency is needed for *de novo* assemblies as compared to
      alignments to a reference genome. However, low quality data can
      affect downstream analysis for *de novo* and reference-based approaches, producing false positives, such as errant SNP predictions.

4.We will use the Stacks’s program **process_radtags** to remove low quality sequences (also known as cleaning data) and to demultiplex our samples. [Here is the Stacks [anual](http://catchenlab.life.illinois.edu/stacks/manual/) as well as the specific [manual page for
process_radtags](http://catchenlab.life.illinois.edu/stacks/manual/#procrad) on the Stacks website to find information         and examples. 
    
  • Get back into your ```dataprep``` folder by running:
  
  cd ```dataprep```
    
  • It is time to load the ```stacks``` module to be able to access the ```process_radtags``` command. Find it load it.
  
  ```
  module spider Stacks
  module load Stacks
  ```
  
   • You will need to specify the set of barcodes used in the construction of the RAD library.
        Remember, each P1 adaptor in RAD has a particular DNA sequence (an inline
        barcode) that gets sequenced first, allowing data to be associated with samples such as
        individuals or populations.
    
   • To save you some time, the barcode file is  ```/scale_wlg_persistent/filesets/project/nesi02659/obss_2020/resources/day3/lane1_barcodes.x``` Copy it in `dataprep` Should you need a quick check, 
   

   •  Normally, these sample names would
        correspond to the individuals used in a particular experiment (e.g. individual ID etc), but for this exercise, 
         samples are named in a simple way. Have a look at the inside of this file using the `less` command.
        
    
   • Based on the barcode file, can you check how many samples were multiplexed together in this
        RAD library? (*Hint:* you can count the lines in the file using `wc -l lane1_barcodes.txt`)
        
      
   • Have a look at the [help online](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) to prepare          your `process_radtags` command.  You will need to specify the restriction enzyme used to construct the library          (SbfI), the directory of input files (the ```lane1``` directory), the list of barcodes, the output directory
        (```samples```), and specify that process_radtags ```clean, discard, and rescue reads.``` as options of                 `process_radtags`. 
        
   • You should now be able to run the ```process_radtags``` command from the ```dataprep``` directory. It will take a couple of minutes to run. 
   
      -   If you find that something is possibly missing from your process_radtags
                input, correct the error and give running process_radtags another try.

   • The process_radtags program will write a log file into the output directory. Have a look in there.
        Examine the log and answer the following questions:
    
    -   How many reads were retained?
    
    -   Of those discarded, what were the reasons? 
    
    -   In the process_radtags log file, what can the list of “sequences not recorded” tell
                you about the barcodes analyzed and about the sequencing quality in general?


Well done! Have a breath, sit back or help your neighbour, we will be back shortly!
