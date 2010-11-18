#!/usr/bin/perl

# Usage:runReleaseSanityCheck  [-sting <path to GATK jar>]  [-dry]

use Getopt::Long;

$dry;
$sting = "/humgen/gsa-scr1/ebanks/Sting_dev/dist/GenomeAnalysisTK.jar";
GetOptions( "dry!" => \$dry,
            "sting=s" => \$sting);

$command_prefix = "java -Xmx4096m -jar $sting -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta -l OFF"; 

$random_number = rand();
$tmp_bam = "/tmp/$random_number.bam";

print "Executing DepthOfCoverage...";
$command = "$command_prefix -T DepthOfCoverage -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.ESP.WEx.chr1.bam -L chr1:10000000-10100000 -o /dev/null";
run($command, $dry);

print "Executing CountCovariatesWholeExome...";
$command = "$command_prefix -T CountCovariates -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.ESP.WEx.chr1.bam -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -L /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/whole_exome_agilent_designed_120.targets.chr1.interval_list -standard -OQ -recalFile /dev/null -XL chr1:1,000,000-247179187";
run($command, $dry);

print "Executing CountCovariatesWholeGenome...";
$command = "$command_prefix -T CountCovariates -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.GAII.chr1.50MB.bam -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -L chr1:1,900,000-2,000,000 -standard -OQ -recalFile /dev/null";
run($command, $dry);

print "Executing TableRecalibratorWholeExome...";
$command = "$command_prefix -T TableRecalibration -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.ESP.WEx.chr1.bam -L /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/whole_exome_agilent_designed_120.targets.chr1.interval_list -OQ -recalFile /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.ESP.WEx.chr1.recal.csv --out $tmp_bam -XL chr1:1,000,000-247179187";
run($command, $dry);

print "Executing TableRecalibratorWholeGenome...";
$command = "$command_prefix -T TableRecalibration -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.GAII.chr1.50MB.bam -L chr1:1,950,000-2,000,000 -OQ -recalFile /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.GAII.chr1.50MB.recal.csv --out $tmp_bam";
run($command, $dry);

print "Executing IndelRealignerWholeExome...";
$command = "$command_prefix -T IndelRealigner -LOD 5 -maxConsensuses 100 -greedy 100 -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.ESP.WEx.chr1.bam -L chr1:900,000-1,000,000 -compress 1 -targetIntervals /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.ESP.WEx.chr1.realigner.intervals -o $tmp_bam";
run($command, $dry);

print "Executing IndelRealignerWholeGenome...";
$command = "$command_prefix -T IndelRealigner -LOD 5 -maxConsensuses 100 -greedy 100 -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.GAII.chr1.50MB.bam -L chr1:975,000-1,000,000 -compress 1 -targetIntervals /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.GAII.chr1.50MB.realigner.intervals -o $tmp_bam";
run($command, $dry);

print "Executing UnifiedGenotyperWholeExome...";
$command = "$command_prefix -T UnifiedGenotyper -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.ESP.WEx.chr1.bam -L /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/whole_exome_agilent_designed_120.targets.chr1.interval_list -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -o /dev/null -XL chr1:1,500,000-247179187";
run($command, $dry);

print "Executing UnifiedGenotyperWholeGenome...";
$command = "$command_prefix -T UnifiedGenotyper -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.GAII.chr1.50MB.bam -L chr1:750,000-1,000,000 -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -o /dev/null";
run($command, $dry);

print "Executing UnifiedGenotyperWholeGenomeMultithreaded...";
$command = "$command_prefix -T UnifiedGenotyper -I /humgen/gsa-hpprojects/GATK/data/Evaluation_Data/NA12878.GAII.chr1.50MB.bam -L chr1:500,000-1,000,000 -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -o /dev/null -nt 4";
run($command, $dry);

unlink $tmp_bam;

sub run {

    $command = $_[0];
    $dry = $_[1];

    local $start = time;

    if ($dry) {
	print "$command\n";
    } else {
	system($command);
    }

    $total_time = time - $start;
    print " [$total_time sec]\n";
}
