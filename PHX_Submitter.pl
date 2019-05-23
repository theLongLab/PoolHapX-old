#!/usr/bin/perl

# $ARGV[0] = The main working directory
# $ARGV[1] = The user-chosen prefix
# $ARGV[2] = Filename of the reference sequence
# $ARGV[3] = Number of pools

# Checklist: Have you...
# 1) ...checked to see if your job scheduler can run job arrays?

my $master_dir = $ARGV[0];
my $pools = $ARGV[3];

my $inter_dir = "$master_dir/intermediate/";
my $prefix = $ARGV[1];
my $ref = "$master_dir/input/";
$ref .= $ARGV[2]; 

# Step 0) Submit the jobs for Steps 1 and 2. 
print STDOUT "\nStep 1) For each pool, align the simulated reads to a reference sequence, and call variants using GATK HaplotypeCaller in gVCF mode.\n";
my $step3 = 0;
my $array_end = $pools - 1;
my $sub3 = `sbatch --array=0-$array_end -J Step1_$ARGV[1] -o Step1_$ARGV[1] PHX_Step1.sh $inter_dir/$prefix $ref`;
$sub3 =~ /^Submitted batch job (\d+)/; 
$step3 = $1;
print STDOUT "The array job ID was $step3.";

print STDOUT "\nStep 2) Complete the joint variant-calling process and run PoolHapX on the output from the pre-processing pipeline.\n";
system("sbatch --dependency=afterok:$step3 -J Step2_$ARGV[1] -o Step2_$ARGV[1] PHX_Step2.pl $master_dir $prefix $pools");

