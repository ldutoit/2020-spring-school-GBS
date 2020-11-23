## Mapping

This small file goes through the mapping of the 30 stickleback samples to the reference genome.

First, we link to the unmapped files, and create a working directory called `oregon_stickleback_mapped`
```bash
mkdir oregon_stickleback_mapped/
ln -s /scale_wlg_persistent/filesets/project/nesi02659/obss_2020/resources/day3/oregon_stickleback .
cd oregon_stickleback_mapped
```

The reference genome is at [http://asia.ensembl.org/Gasterosteus_aculeatus/Info/Index](http://asia.ensembl.org/Gasterosteus_aculeatus/Info/Index)

We can download the repeat masked file that includes all chromosomes from this link. The download command making use of the `wget` command is below.

```bash
wget ftp://ftp.ensembl.org/pub/release-101/fasta/gasterosteus_aculeatus/dna/Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa.gz
gunzip Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa.gz
```

We then create the bwa index.

```
module load BWA
bwa index Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa
```

The loop below runs the alignment for each file. It should be run in a batch job file.

```bash
 module load BWA SAMtools
 for file in ../oregon_stickleback/*gz
  do echo $file
  sample=$(basename ${file} .fa.gz) #filename is file without the extension
    bwa mem Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa $file  | samtools sort  | samtools view -hb >  ${sample}.bam #output a sorted bam
  done	
```

Those BAM files are now ready to be processed by ref_map.pl!


