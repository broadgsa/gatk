import java.io.File
import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools._

class ManySampleUGPerformanceTesting extends QScript {
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("/home/radon01/depristo/dev/GenomeAnalysisTKStable/trunk/dist/GenomeAnalysisTK.jar")

  @Argument(shortName = "R", doc="ref", required=false)
  var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  val TARGET_INTERVAL = "my.intervals"
  val TEST_BAM_LIST = new File("ten.bam.list")
  val FULL_BAM_LIST = new File("/humgen/1kg/processing/allPopulations_chr20_june_release/allPopulations.june.bam.list")
  val BAM_LIST = FULL_BAM_LIST
  val HM3 = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr.hg19_fwd.vcf")

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.intervals :+= new File(TARGET_INTERVAL);
    this.reference_sequence = referenceFile;
    this.jobQueue = "gsa";
    this.et = Option(org.broadinstitute.sting.gatk.phonehome.GATKRunReport.PhoneHomeOption.STANDARD);
    this.dcov = Option(50);
  }

  def script = {
    // SNP calling
    add(new MyQSample(BAM_LIST));
  }

  class MyQSample(@Input(doc="foo") bamList: File) extends QSample with UNIVERSAL_GATK_ARGS {
    this.memoryLimit = Some(4)
    this.input_file :+= bamList
    //this.BTI = "genotypes"
    this.nt = Option(10)
    this.rodBind :+= RodBind("genotypes", "VCF", HM3)
    this.o = new File("%s.qsample".format(bamList.getName))
  }
}

