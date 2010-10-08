import java.io.{FileReader, File, BufferedReader}
import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType

class omni_qc extends QScript {
  qscript =>

  // NON-OMNI VCF FILES
  var pilot3_release_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/pilot3/merge_release/ALL.exon.2010_03.genotypes.vcf","vcf")
  var pilot1_ceu_vcf = new TaggedFile("/humgen/1kg/releases/pilot_project/2010_07/low_coverage/snps/CEU.low_coverage.2010_07.genotypes.vcf.gz","vcf")
  var pilot1_chb_vcf = new TaggedFile("/humgen/1kg/releases/pilot_project/2010_07/low_coverage/snps/CHBJPT.low_coverage.2010_07.genotypes.vcf.gz","vcf")
  var pilot1_yri_vcf = new TaggedFile("/humgen/1kg/releases/pilot_project/2010_07/low_coverage/snps/YRI.low_coverage.2010_07.genotypes.vcf.gz","vcf")
  var august_calls_vcf = null
  var hiseq_calls_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/omni/resources/NA12878.HiSeq.v9.b36.vcf.gz","vcf")
  var pilot1_with_na12878_vcf = new TaggedFile("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/calls/v2/N60/lowpass.N60.recal.mG6.retranche.vcf","vcf")

  // OMNI VCF FILES
  var OMNI_b36_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/omni/resources_oct7/Omni_2_5_795_b36.vcf","vcf")
  var OMNI_b37_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/omni/resources_oct7/Omni_2_5_795_b37.deduped.vcf.gz","vcf")
  var OMNI_hapmap_b36_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/omni/resources_oct7/Omni_2_5_pilot_b36.vcf.gz","vcf")

  // INTERVALS
  var pilot3_interval_list: String = "/humgen/gsa-hpprojects/1kg/1kg_pilot3/documents/CenterSpecificTargetLists/results/p3overlap.targets.b36.interval_list"
  var pilot1_interval_list: String = _
  var hiseq_interval_list: String = _
  var production_interval_list: String = "20"

  // OTHER
  var failed_qc_samples = new File("/humgen/gsa-scr1/chartl/projects/omni/resources_oct7/failed_qc_samples.txt")

  trait OmniArgs extends CommandLineGATK {
    this.reference_sequence = new File("/humgen/gsa-hpprojects/1kg/reference/human_b36_both.fasta")
    this.jarFile = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")
  }

  protected def sampleFileToString(samFile: File) : List[String] = {
    var reader = new BufferedReader(new FileReader(samFile))
    var line: String = ""
    var sampleList: List[String] = Nil
    line = reader.readLine
    while ( line != null ) {
      sampleList :+= line
      line = reader.readLine
    }

    return sampleList
  }

  class vcf2bed(b: String) extends CommandLineFunction {
    @Input(doc="A VCF file to be put into an interval list") var in_vcf: File = _
    @Output(doc="An interval list file to be used with -L") var out_list: File = _

    def commandLine = "python /humgen/gsa-scr1/chartl/projects/omni/scripts/vcf2bed.py %s %s".format(in_vcf.getAbsolutePath,out_list.getAbsolutePath)
  }

  protected def runme(samples: File, intervals: String, compVCF: TaggedFile, outDir: String, base: String) {
    var subset = new SelectVariants with OmniArgs
    subset.variantVCF = Omni_chip_vcf
    subset.analysisName = "Subset_"+base
    subset.sample = sampleFileToString(samples)
    if ( intervals != null ) {
      subset.intervalsString :+= intervals
    }
    subset.out = new TaggedFile(outDir+base+".subset.vcf","vcf")
    
    var makeOmniSiteList = new vcf2bed("foo")
    makeOmniSiteList.analysisName = "Omni_list_"+base
    makeOmniSiteList.in_vcf = subset.out
    makeOmniSiteList.out_list = new TaggedFile(outDir+base+".subset.omni.intervals.list","intervals.list")

    var eval = new VariantEval with OmniArgs
    eval.rodBind :+= new RodBind("evalOmni","vcf",subset.out)
    eval.rodBind :+= new RodBind("comp"+base,"vcf",compVCF)
    eval.intervals = makeOmniSiteList.out_list
    eval.evalModule :+= "GenotypeConcordance"
    eval.evalModule :+= "SimpleMetricsBySample"
    eval.reportLocation = new TaggedFile(outDir+base+".subset.omni.eval","eval")
    eval.reportType = Some(VE2TemplateType.R)
    eval.analysisName = base+"_Eval"
    eval.select ++= List("AC==1","AC==2","AC==3","AC==4","AC==5","AC>5")
    eval.selectName ++= List("AC1","AC2","AC3","AC4","AC5","AC6+")
    eval.disI = true
    eval.outputVCF = outDir+base+".eval.interesting.sites.vcf"

    if ( base.equalsIgnoreCase("production") ) {
      subset.variantVCF = Omni_b37_vcf
      subset.reference_sequence = new File("/humgen/gsa-hpprojects/1kg/reference/human_g1k_v37.fasta")
      eval.reference_sequence = new File("/humgen/gsa-hpprojects/1kg/reference/human_g1k_v37.fasta")
      eval.DBSNP = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_b37.rod")
    }



    add(subset,makeOmniSiteList,eval)
  }

  def script = {
    runme(pilot1_overlap_samples,pilot1_interval_list,pilot1_vcf,"/humgen/gsa-scr1/chartl/projects/omni/scratch/queue/pilot1/","pilot1")
    runme(pilot3_overlap_samples,pilot3_interval_list,pilot3_vcf_broad,"/humgen/gsa-scr1/chartl/projects/omni/scratch/queue/pilot3/","pilot3")
    runme(hiseq_overlap_samples,hiseq_interval_list,hiSeq_vcf,"/humgen/gsa-scr1/chartl/projects/omni/scratch/queue/hiSeq/","NA12878_HiSeq")
    runme(production_overlap_samples,production_interval_list,production_vcf,"/humgen/gsa-scr1/chartl/projects/omni/scratch/queue/production/","production")
    runme(pilot3_overlap_samples,pilot3_interval_list,pilot3_vcf,"/humgen/gsa-scr1/chartl/projects/omni/scratch/queue/pilot3/","pilot3_release")
  }
}
