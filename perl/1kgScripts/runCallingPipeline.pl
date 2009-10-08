#!/usr/bin/perl -w

use strict;
use Getopt::Long;

sub usage {
    print "Usage: perl runCallingPipeline.pl\n\t-i <GATK input bam command>\n\t-o <output file head>\n\t[-sting Sting dir]\n\t[-frac multiplier for indel fractions]\n\t[-snps Should we call snps?]\n\t[-sample for writing vcf with -snps]\n\t[-badsnps bad snps file from cleaning with -snps]\n\t[-doc DepthOfCoverage for filtering with -snps]\n\t[-mq mapping quality zero for filtering with -snps ]\n\t[-q farm queue; default:gsa]\n\t[-wait farm wait id]\n\t[-dry]\n";
    exit(1);
}

my $inputBamStr = undef;
my $outputHead = undef;
my $wait = undef;
my $dry;
my $snps;
my $badsnps = undef;
my $sample;
my $doc = 100;
my $mq = 100;
my $indelFractionMultiplier = "";
my $queue = "gsa";
my $sting = "/humgen/gsa-scr1/ebanks/Sting";

GetOptions( "i=s" => \$inputBamStr,
	    "o=s" => \$outputHead,
	    "q:s" => \$queue,
	    "dry!" => \$dry,
	    "snps!" => \$snps,
	    "sample:s" => \$sample,
	    "doc:s" => \$doc,
	    "mq:s" => \$mq,
	    "frac:s" => \$indelFractionMultiplier,
	    "badsnps:s" => \$badsnps,
	    "wait:s" => \$wait,
	    "sting:s" => \$sting );

usage() if ( !$inputBamStr || !$outputHead );

my $indelsHigh = "$outputHead.indels.high.calls";
my $bsub = "bsub -q $queue -o $indelsHigh.sdout";
if ($wait) {
    $bsub .= " -w \"ended($wait)\"";
}
my $command = "java -Djava.io.tmpdir=/broad/hptmp/ -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T IndelGenotyper -R /broad/1KG/reference/human_b36_both.fasta $inputBamStr -minConsensusFraction 0.5 -minCnt 2 -1kg -minFraction 0.".$indelFractionMultiplier."3 -O $indelsHigh";
execute("$bsub $command", $dry);

my $indelsLow = "$outputHead.indels.low.calls";
$bsub = "bsub -q $queue -o $indelsLow.sdout -J $outputHead";
if ($wait) {
    $bsub .= " -w \"ended($wait)\"";
}
$command = "java -Djava.io.tmpdir=/broad/hptmp/ -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T IndelGenotyper -R /broad/1KG/reference/human_b36_both.fasta $inputBamStr -minConsensusFraction 0.5 -minCnt 2 -1kg -minFraction 0.".$indelFractionMultiplier."1 -O $indelsLow";
execute("$bsub $command", $dry);

if ($snps) {
    my $snpsFile = "$outputHead.snps.unfiltered.calls";
    $bsub = "bsub -q $queue -o $snpsFile.sdout -J $outputHead";
    if ($wait) {
	$bsub .= " -w \"ended($wait)\"";
    }
    $command = "java -Djava.io.tmpdir=/broad/hptmp/ -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T UnifiedGenotyper -R /broad/1KG/reference/human_b36_both.fasta $inputBamStr -varout $snpsFile -lod 0.0";
    execute("$bsub $command", $dry);

    my $filterFile = "$outputHead.snps.filtered.calls";
    my $vcfFile = "$outputHead.snps.vcf";
    $bsub = "bsub -q $queue -o $filterFile.sdout -w \"ended($outputHead)\"";
    $command = "java -Djava.io.tmpdir=/broad/hptmp/ -Xmx4096m -jar $sting/dist/GenomeAnalysisTK.jar -S SILENT -T VariantFiltration -R /broad/1KG/reference/human_b36_both.fasta $inputBamStr -vcf $vcfFile -included $filterFile -sample $sample -B dbsnp,dbsnp,/humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod,variant,Variants,$snpsFile,";
    if ($badsnps) {
	$command .= "cleaned,CleanedOutSNP,$badsnps,";
    }
    $command .= "indels,SimpleIndel,$indelsLow -X DepthOfCoverage:max=$doc -X AlleleBalance:low=0.25,high=0.75 -X FisherStrand:pvalue=0.00001 -X LodThreshold:lod=5.0 -X MappingQualityZero:max=$mq -X IndelArtifact -X ClusteredSnps:window=7,snps=3";
    execute("$bsub $command", $dry);
}

sub execute {

    my $cmd = $_[0];
    my $dry = $_[1];

    if ($dry) {
	print "$cmd\n";
    } else {
	system($cmd);
    }
}
