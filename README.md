# PoolHapX: De novo Haplotype Reconstruction from Pooled NGS Data

To tackle the high-dimensional, multifaceted problem of haplotype reconstruction, we have developed PoolHapX, a computational program unifying short-range genomic information in sequencing reads, and long-range evidence for linkage disequilibrium in pooled populations. PoolHapX harnesses both 1) between-host information such as shared haplotypes distributed across several within-host pathogen populations, in conjunction with 2) within-host information. Long-range haplotypes are resolved de novo from short NGS reads generated from multiple complex mixtures of individual particles, which will be referred to as pools hereafter. 

The union of multiple strategies combines the strengths and compensates for the individual weaknesses of genomics-based algorithms and population genetics-based data-sharing programs. PoolHapX is written in Java (version 1.8), which is operating system-agnostic, and is distributed as a JAR file, allowing users to efficiently use the program out-of-the-box.

<img src="https://github.com/theLongLab/PoolHapX/blob/master/Figure_3.png" width="900">

**Overview of the PoolHapX methodology:** The PoolHapX (PoolHapX) program approximates the genotypic resolution of single-cell sequencing using only pooled sequencing data by integrating population genetics models with genomics algorithms to reconstruct haplotypes. A) PoolHapX first determines locations of physical linkage uncertainty using short NGS reads, and then divides the full genome into shorter regions. B) Regional haplotypes are solved for and joined together for a parsimonious global distribution of haplotypes. C) The within-pool frequency of each haplotype is estimated by regression to solve for each within-pool haplotype distribution.

## Getting Started

### Installation
Create a PoolHapX directory and download the source code from repository. 

### Prerequisites
Java 1.8+: https://www.java.com/en/download/
If there is no data-specific variant-caller needed:
	* bwa: https://github.com/lh3/bwa
	* samtools 0.1.19+: http://www.htslib.org/download/
	* GATK 4+: https://software.broadinstitute.org/gatk/download/index
	* A job scheduler capable of running job arrays and job dependencies (ex. Slurm). Makes it easy to a) run PHX_Step1.pl in parallel for each sample and b) automate PHX_Step2.pl to run immediately after PHX_Step1.pl finishes successfully. 

## How to Use

**PoolHapX Input:** Allele-annotated reads and observed allele frequencies. For examples, see 0_0_p0.vef and 0_0_vars.intra_freq.txt in sample_io_files/. 

**PoolHapX Output:** Within-pool haplotype distributions (allelic compositions, within-pool and across-pool frequencies). For examples, see 0_0.inter_freq_vars.txt and 0_0.intra_freq.txt in sample_io_files/. 

### What is your input data? 
*Option 1:* Raw or aligned sequencing reads (FASTQ)
*Option 2:* Sorted and aligned sequencing reads (BAM) and segregating sites called from data-specific variant-caller (VCF).
* If you have two or fewer FASTQ or BAM files, it is possible to run the PoolHapX workflow on a local machine (ex. laptop) using the command line without running out of memory. Instead of submitting PHX_Submitter.pl, simply copy and paste the commands from PHX_Step1.pl and PHX_Step2.pl. 

### Running Option 1 or 2:
1. Create the following subdirectories in your main working directory: input, intermediate, output, programs. These will correspond to the organizational structures referred to in all of the scripts.  
2. Place PHX.properties and your reference sequence in input/. 
3. PHX.properties is a JSON-like file setting parameters for PoolHapX.jar. Update the absolute paths of the directories. For an explanation of each parameter, see below. 
4. Decide on a common prefix across all of your sample FASTQs and name them as: `{$prefix}\_p{$sample\_number}\_read\_{1_or_2}.fastq`. Place them in intermediate/.

### Running Option 1: 
1. Index the reference sequence file with bwa, samtools, and Picard.
2. Place PHX_Submitter.pl, PHX_Step1.pl and PHX_Step2.pl in your main working directory. 
3. If job arrays are not possible, alter PHX_Submitter.pl such that it submits PHX_Step1.pl as multiple jobs instead of a parallel job array, and that PHX_Step2.pl is not submitted dependent on the variable `$SLURM_ARRAY_TASK_ID` (i.e.: that all PHX_Step1.pl jobs finish).
4. In PHX_Step1.pl, replace the directories for bwa, Samtools, Java, and GATK. If you already have BAM files, comment out the bwa, samtools, and GATK AddOrReplaceReadGroups commands as necessary. Alter the GATK variant-calling parameters as necessary. Alter the memory requirements as necessary (current setting has been tested with GATK on up to 80 segregating sites).
	* Check the version of Java installed in the compute nodes specifically, not the head node. If it is 1.8+, simply replace `/path/to/java` with java. If not, replace with the absolute path to your personally installed copy of Java.
	* BAM files must be sorted and indexed to be processed by GATK. 
	* If your BAM files already have read-group (@RG) information, comment out step 1.1.2. 
5. In PHX_Step2.pl, replace the directories for Java, GATK, and Samtools. Replace the name of the reference sequence file at $ref. Alter the memory requirements as necessary (current setting has been tested with PoolHapX on up to 80 segregating sites).
6. Place BAMFormatter_GATK.jar, PairedReadLinker.jar, and PoolHapX.jar in programs/.
7. Make sure the FASTQ files are unzipped (if applicable).
8. Submit PHX_Submitter.pl to the job scheduler. The result files (allele composition and within/across-pool frequencies) will show up in output/. For example, if a Slurm job scheduler is available:
`cd main_work_directory
sbatch PHX_Submitter.pl main_work_dir prefix ref_seq num_pools`

### Running Option 2: 
1. In PHX_Step2_Alt.pl, replace the directories for Java and Samtools. Alter the memory requirements as necessary (current setting has been tested with PoolHapX on up to 80 segregating sites).
	* Check the version of Java installed in the compute nodes specifically, not the head node. If it is 1.8+, simply replace `/path/to/java` with java. If not, replace with the absolute path to your personally installed copy of Java.
3. Place BAMFormatter_Alt.jar, PairedReadLinker.jar, and PoolHapX.jar in programs/.
3. Submit PHX_Step2_Alt.pl to the job scheduler. The result files (allele composition and within/across-pool frequencies) will show up in output/. 
`cd main_work_directory
sbatch PHX_Step2_Alt.pl main_work_dir prefix num_pools`

### Troubleshooting:

**Q:** PoolHapX fails to complete the first graph-colouring step within a day. What's happening?
**A:** The linkage in allele-annoted reads may have reached saturation, and graph-colouring cannot handle a) that many reads with b) that many segregating sites. There are several ways to solve this problem. 
	* Reducing the number of segregating sites: If you suspect that some segregating sites are not relevant to your investigation (ex. mutations in the middle of an eukaryotic intron), then filtering them at the GATK SelectVariant step will be useful. 
	* Reducing the read coverage: If there are simply too many reads, there will be too many edges drawn and far, far too many possible connections for PoolHapX to handle. You can down-sample the SAM file with a blunt coverage threshold, thereby keeping proportionally more reads in areas with already-low coverage. Follow the instructions here: http://lindenb.github.io/jvarkit/Biostar154220.html, with further explanations here: https://www.biostars.org/p/154220/. 
	* The output of either of these two solutions will be a SAM file and a VCF file. Run the appropriate version of BAMFormatter, PairedReadLinker, and PoolHapX as normal afterwards.

## Authors

* **Lauren Mak:** Developed the PoolHapX algorithm, ran all tests and comparisons to other haplotype reconstruction programs on simulated and real data.
* **Chen Cai:** Implemented the initial version of the graph-colouring algorithm. 
* **Jingni He:** Wrote scripts and programs to assist with comparisons to other haplotype reconstruction programs. 
* **Quan Long:** Principal investigator. Developed the original PoolHap algorithm, implemented the initial versions of the divide-and-conquer and linkage-informed LASSO regression algorithms.

## MIT Open-Source License

This project is licensed under the MIT License. See the [license.md](license.md) file for details.

## Citation

Please cite (reference publication to be added). 

***TODO:***
* Describe parameters. 
