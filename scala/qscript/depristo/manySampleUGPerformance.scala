import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction

class ManySampleUGPerformanceTesting extends QScript {
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

val TARGET_INTERVAL = "my.intervals"
val FULL_BAM_LIST = new File("/humgen/1kg/processing/allPopulations_chr20_june_release/allPopulations.june.bam.list")
val MERGED_DIR = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/manySampleUGPerformance/")

trait UNIVERSAL_GATK_ARGS extends CommandLineGATK { 
    this.logging_level = "INFO"; 
    this.jarFile = gatkJarFile;
    this.intervals = new File(TARGET_INTERVAL);
    this.reference_sequence = referenceFile; 
    this.jobQueue = "gsa"; 
    this.et = Option(org.broadinstitute.sting.gatk.phonehome.GATKRunReport.PhoneHomeOption.STANDARD);
    this.dcov = Option(50);
    //this.commandDirectory = new File("results");
    }

def script = {
    for (nSamples <- List(1, 2, 5, 10, 50, 100, 200, 300, 400, 500)) {
        val sublist = new SliceList(nSamples)
        val mergeSublist = new MergeBAMs(sublist.list)
        
    	add(sublist)
    	add(mergeSublist)
    	add(new Index(mergeSublist.o) )
    	add(new Call(sublist.list, nSamples, "dynamic_merge"))
    	add(new Call(mergeSublist.o, nSamples, "pre_merge"))
    }
}

class Index(bamIn: File) extends SamtoolsIndexFunction {
    this.jobQueue = "gsa"
    bamFile = bamIn
}

class MergeBAMs(bamList: File) extends PrintReads with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = Some(3)
    this.input_file :+= bamList
    this.o = new File(MERGED_DIR + "/" + bamList.getName + ".bam")
  }

class Call(@Input(doc="foo") bamList: File, n: Int, name: String) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outVCF: File = new File("%s.%d.%s.vcf".format(bamList.getName, n, name))
    this.memoryLimit = Some(4)
    this.input_file :+= bamList
    this.jobQueue = "gsa"
    this.stand_call_conf = Option(10.0)
    this.o = outVCF
  }

class SliceList(n: Int) extends CommandLineFunction {
    @Output(doc="foo") var list: File = new File("bams.%d.list".format(n))
    this.jobQueue = "gsa"
    def commandLine = "head -n %d %s > %s".format(n, FULL_BAM_LIST, list)
  }
}

