import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.extensions.gatk._

class ManySampleUGPerformanceTesting extends QScript {
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  @Argument(shortName = "bams", doc="BAMs", required=true)
  val FULL_BAM_LIST: File = null;

  @Argument(shortName = "intervals", doc="intervals", required=true)
  val TARGET_INTERVAL: String = null;

  @Argument(shortName = "preMerge", doc="preMerge", required=false)
  val PRE_MERGE: Boolean = false;

  @Argument(shortName = "dcov", doc="dcov", required=false)
  val DCOV: Int = 50;

  @Argument(shortName = "exome", doc="exome ",required=false)
  val EXOME_NSAMPLES: Boolean = false;

  val MERGED_DIR = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/manySampleUGPerformance/")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.intervals = List(new File(TARGET_INTERVAL));
    this.reference_sequence = referenceFile;
    //this.jobQueue = "gsa";
    this.memoryLimit = 4
    //this.commandDirectory = new File("results");
  }

  def script = {
    for (nSamples <- if ( EXOME_NSAMPLES) List(1, 2, 5, 10, 25, 50, 100, 200, 300, 400, 500) else List(1, 2, 5, 10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900) ) {
//    for (nSamples <- List(10)) {
      val sublist = new SliceList(nSamples)
      val mergeSublist = new MergeBAMs(sublist.list)

      val name: String = if ( PRE_MERGE ) "pre_merge" else "dynamic_merge"
      val bams: File = if ( PRE_MERGE ) mergeSublist.o else sublist.list

      add(sublist)
      if ( PRE_MERGE ) {
        add(mergeSublist)
        add(new Index(mergeSublist.o) )
      }

      // SNP calling
      //add(new Call(sublist.list, nSamples, "dynamic_merge"))
      add(new Call(bams, nSamples, name));

      val gtWithBAQ = new Call(bams, nSamples, name + "_baq");
      gtWithBAQ.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE
      add(gtWithBAQ)

      // SNP calling -- no annotations
      //add(new Call(bams.list, nSamples, "dynamic_merge_no_annotations") { this.G :+= "None"; })

      // CountLoci
      //add(new MyCountLoci(sublist.list, nSamples, "dynamic_merge"))
      add(new MyCountLoci(bams, nSamples, name))
    }
  }

  class Index(bamIn: File) extends SamtoolsIndexFunction {
    //this.jobQueue = "gsa"
    bamFile = bamIn
  }

  class MergeBAMs(bamList: File) extends PrintReads with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = 3
    this.input_file :+= bamList
    this.memoryLimit = 16
    this.o = new File(MERGED_DIR + "/" + bamList.getName + ".bam")
  }

  class Call(@Input(doc="foo") bamList: File, n: Int, name: String) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outVCF: File = new File("%s.%d.%s.vcf".format(bamList.getName, n, name))
    this.input_file :+= bamList
    this.stand_call_conf = 10.0
    this.dcov = DCOV;
    this.o = outVCF
  }

  class MyCountLoci(@Input(doc="foo") bamList: File, n: Int, name: String) extends CountLoci with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outFile: File = new File("%s.%d.%s.txt".format(bamList.getName, n, name))
    this.input_file :+= bamList
    this.dcov = DCOV;
    this.o = outFile
  }

  class SliceList(n: Int) extends CommandLineFunction {
    @Output(doc="foo") var list: File = new File("bams.%d.list".format(n))
    def commandLine = "head -n %d %s > %s".format(n, FULL_BAM_LIST, list)
    //this.jobQueue = "gsa";
  }
}

