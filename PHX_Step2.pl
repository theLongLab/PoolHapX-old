#!/usr/bin/perl
#SBATCH --mem=40gb	 		# Job memory request

# $ARGV[0] = The main working directory
# $ARGV[1] = The user-chosen prefix
# $ARGV[2] = Number of pools

# Checklist: Have you...
# 1) ...installed and added the paths to samtools and GATK?
# 2) ...checked the version of Java available on the compute nodes? 
# 3) ...altered the GATK variant-calling parameters, if necessary?

# Step 2.0) Standardize filenames. Make sure all intermediate files are put in intermediate, and all programs can be accessed. 
my $master_dir = $ARGV[0];
my $java_cmd = "/path/to/java";
my $gatk_cmd = "/path/to/gatk";
my $samtools_cmd = "/path/to/samtools"; 
my $own_pdir = "$master_dir/programs";	# The directory within main containing PoolHapX and accessory programs.
my $input_dir = "$master_dir/input";	# The directory within main containing the reference sequence and PoolHapX properties file.
my $inter_dir = "$master_dir/intermediate";
my $prefix = $ARGV[1];	# Corresponds to prefix in the previous script. 
my $pools = $ARGV[2];

my $ref =  "$input_dir/HIV_HXB2.fa";
my $gVCFs; 
for (my $p = 0; $p < $pools; $p++) {
	$gVCFs .= "-V $inter_dir/$prefix\_p$p.raw.g.vcf "; 
}

# Step 2.1.1) Join all of the pool-specific gVCFs into a joint gVCF file usign Combine GVCFs.
print STDOUT "Step 2.1.1) Join all of the pool-specific gVCFs into a joint gVCF file usign Combine GVCFs.\n\n";
system("$java_cmd -jar $gatk_cmd CombineGVCFs -R $ref $gVCFs -O $inter_dir/$prefix.g.vcf");

# Step 2.1.2) Convert the joint gVCF into a VCF file using GenotypeGVCFs.
print STDOUT "\nStep 2.1.2) Convert the joint gVCF into a VCF file using GenotypeGVCFs.\n\n";
system("$java_cmd -jar $gatk_cmd GenotypeGVCFs -R $ref -V $inter_dir/$prefix.g.vcf -ploidy 150 -O $inter_dir/$prefix.raw.vcf"); 
system("$java_cmd -jar $gatk_cmd SelectVariants -R $ref -V $inter_dir/$prefix.raw.vcf -O $inter_dir/$prefix.vcf"); 

# Step 2.2) Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files. 
print STDOUT "\nStep 2.2) Convert each pool-specific BAM file to SAM (i.e.: text), then VEF files, then paired-read VEF files.\n\n";
for (my $p = 0; $p < $pools; $p++) { 
	my $pfile = "$inter_dir/$prefix\_p$p";
	system("$samtools_cmd view -ho $pfile.srt.sam $pfile.srt.bam");
}
system("$java_cmd -jar $own_pdir/BAMFormatter_GATK.jar $inter_dir/ $prefix $pools");
for (my $p = 0; $p < $pools; $p++) {
    my $pfile = "$inter_dir/$prefix\_p$p";
    system("$java_cmd -jar $own_pdir/PairedReadLinker.jar $pfile");
}

# Step 2.3) Reconstructing haplotypes from multi-pool VEF files.
print STDOUT "\nStep 2.3) Reconstructing haplotypes from multi-pool VEF files...\n\n";
system("$java_cmd -jar $own_pdir/PoolHapX.jar $input_dir/PHX.properties $prefix $pools");