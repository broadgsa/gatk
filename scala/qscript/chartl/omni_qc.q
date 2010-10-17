import java.io.{FileReader, File, BufferedReader}
import net.sf.picard.reference.FastaSequenceFile
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils
import org.broadinstitute.sting.gatk.DownsampleType
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel.Model
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType
import scala.collection.mutable.HashMap

class omni_qc extends QScript {
  qscript =>

  // NON-OMNI VCF FILES
  var pilot3_release_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/pilot3/merge_release/ALL.exon.2010_03.genotypes.vcf","vcf")
  var pilot1_ceu_vcf = new TaggedFile("/humgen/1kg/releases/pilot_project/2010_07/low_coverage/snps/CEU.low_coverage.2010_07.genotypes.vcf.gz","vcf")
  var pilot1_chb_vcf = new TaggedFile("/humgen/1kg/releases/pilot_project/2010_07/low_coverage/snps/CHBJPT.low_coverage.2010_07.genotypes.vcf.gz","vcf")
  var pilot1_yri_vcf = new TaggedFile("/humgen/1kg/releases/pilot_project/2010_07/low_coverage/snps/YRI.low_coverage.2010_07.genotypes.vcf.gz","vcf")
  var august_calls_EUR = new TaggedFile("/humgen/1kg/processing/release/august/EUR.vcf","vcf")
  var august_calls_ASN = new TaggedFile("/humgen/1kg/processing/release/august/ASN.vcf","vcf")
  var august_calls_AFR = new TaggedFile("/humgen/1kg/processing/release/august/AFR.vcf","vcf")
  var august_calls_EUR_refined = new TaggedFile("/humgen/1kg/processing/release/august/EUR.beagle.vcf.gz","vcf")
  var august_calls_ASN_refined = new TaggedFile("/humgen/1kg/processing/release/august/ASN.beagle.vcf.gz","vcf")
  var august_calls_AFR_refined = new TaggedFile("/humgen/1kg/processing/release/august/AFR.beagle.vcf.gz","vcf")
  var hiseq_calls_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/omni/resources/NA12878.HiSeq.v9.b36.vcf.gz","vcf")
  var pilot1_with_na12878_vcf = new TaggedFile("/humgen/1kg/analysis/bamsForDataProcessingPapers/lowpass_b36/calls/v2/N60/lowpass.N60.recal.mG6.retranche.vcf","vcf")
  //var august_calls_other_genotypes = _

  // OMNI VCF FILES
  var OMNI_b36_vcf = new TaggedFile("/humgen/illumina/1kg_seq_vcfs/Illumina_HapMap_Omni_2.5_764samples.vcf","vcf")
  var OMNI_b37_vcf = new TaggedFile("/broad/shptmp/chartl/Omni_2.5_764_samples.b37.deduped.vcf","vcf")
  var OMNI_hapmap_b36_vcf = new TaggedFile("/humgen/gsa-scr1/chartl/projects/omni/resources_oct7/Omni_2_5_pilot.b36.vcf","vcf")

  // INTERVALS
  var pilot3_interval_list: String = "/humgen/gsa-hpprojects/1kg/1kg_pilot3/documents/CenterSpecificTargetLists/results/p3overlap.targets.b36.interval_list"
  var pilot1_interval_list: String = "/broad/shptmp/chartl/omni/resources/Omni_b36_sites.interval.list"
  var hiseq_interval_list: String = "/broad/shptmp/chartl/omni/resources/Omni_b36_sites.interval.list"
  var production_interval_list: String = "/broad/shptmp/chartl/omni/resources/Omni_b37_sites.chr20.interval.list"

  // REFERENCES
  var b36_ref = new File("/humgen/1kg/reference/human_b36_both.fasta")
  var b37_ref = new File("/humgen/1kg/reference/human_g1k_v37.fasta")

  // OTHER
  val analysis_dir = "/broad/shptmp/chartl/omni/"
  val resources_dir = analysis_dir + "resources/"
  val scratch_dir = analysis_dir + "scratch/"
  val eval_dir = analysis_dir + "eval/"
  val vcf_dir = analysis_dir + "vcfs/"

  trait OmniArgs extends CommandLineGATK {
    this.jarFile = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")
  }

  class vcf2bed extends CommandLineFunction {
    @Input(doc="A VCF file to be put into an interval list") var in_vcf: File = _
    @Output(doc="An interval list file to be used with -L") var out_list: File = _

    def commandLine = "python /humgen/gsa-scr1/chartl/projects/omni/scripts/vcf2bed.py %s %s".format(in_vcf.getAbsolutePath,out_list.getAbsolutePath)
  }

  class GetSampleOverlap extends CommandLineFunction {
    @Input(doc="A list of VCF files for which to calculate the sample overlap") var in_vcfs: List[File] = Nil
    @Output(doc="A file to which to write the overlapping sample names") var outFile: File = _

    /*def commandLine = "grep #CHR %s | sed 's/.vcf:/\\t/g' | cut -f11- | tr '\\t' '\\n' | sort | uniq -c | awk '$1 == %d' | awk '{print $2}' > %s".format(
    in_vcfs.foldLeft[String]("")( (str,f) => if ( str.equals("") ) str + f.getAbsolutePath else str + " " + f.getAbsolutePath),
    in_vcfs.size,
    outFile.getAbsolutePath
    )*/
    def commandLine = "python /humgen/gsa-scr1/chartl/projects/omni/scripts/getOverlapSamples.py %s %s".format(
      in_vcfs.foldLeft[String]("")( (str,f) => if ( str.equals("") ) str + f.getAbsolutePath else str + " " + f.getAbsolutePath),
      outFile.getAbsolutePath
      )
  }

  class GunzipFile extends CommandLineFunction {
    @Input(doc="file to gunzip") var gunzipMe: File = _
    @Output(doc="file to gunzip to") var outFile: File = _

    def commandLine = "gunzip -c %s > %s".format(gunzipMe.getAbsolutePath,outFile.getAbsolutePath)
  }

  def script = {

    /** Convert other chips to merged VCFs **/

    //var august_call_other_chips: List[(String,File)] = processAuxiliaryChipData(august_calls_other_genotypes)


    /** Unzip the pilot 1 VCFs and dump them into the resources directory **/

    var gunzip_p1_ceu = new GunzipFile
    var gunzip_p1_chb = new GunzipFile
    var gunzip_p1_yri = new GunzipFile
    var gunzip_hiseq = new GunzipFile
    var gunzip_ag_eur = new GunzipFile
    var gunzip_ag_asn = new GunzipFile
    var gunzip_ag_afr = new GunzipFile

    gunzip_p1_ceu.gunzipMe = pilot1_ceu_vcf
    gunzip_p1_ceu.outFile = new File(resources_dir+"CEU.low_coverage.genotypes.vcf")
    gunzip_p1_chb.gunzipMe = pilot1_chb_vcf
    gunzip_p1_chb.outFile = new File(resources_dir+"CHB.low_coverage.genotypes.vcf")
    gunzip_p1_yri.gunzipMe = pilot1_yri_vcf
    gunzip_p1_yri.outFile = new File(resources_dir+"YRI.low_coverage.genotypes.vcf")
    gunzip_hiseq.gunzipMe = hiseq_calls_vcf
    gunzip_hiseq.outFile = new File(resources_dir+"HiSeq.b36.vcf")
    gunzip_ag_eur.gunzipMe = august_calls_EUR_refined
    gunzip_ag_eur.outFile = new File(resources_dir+"EUR.refined.vcf")
    gunzip_ag_asn.gunzipMe = august_calls_ASN_refined
    gunzip_ag_asn.outFile = new File(resources_dir+"ASN.refined.vcf")
    gunzip_ag_afr.gunzipMe = august_calls_AFR_refined
    gunzip_ag_afr.outFile = new File(resources_dir+"AFR.refined.vcf")

    add(gunzip_p1_ceu,gunzip_p1_yri,gunzip_p1_chb,gunzip_hiseq,gunzip_ag_eur,gunzip_ag_asn,gunzip_ag_afr)

    /** fix the omni ref bases **/
    var fix_421 = new FixRefBases with OmniArgs
    var fix_764 = new FixRefBases with OmniArgs
    var fix_764_b37 = new FixRefBases with OmniArgs

    fix_421.variantVCF = OMNI_hapmap_b36_vcf
    fix_421.reference_sequence = b36_ref
    fix_421.out = new File(vcf_dir+swapExt(OMNI_hapmap_b36_vcf.getName,".vcf",".ref_fixed.vcf"))
    fix_764.variantVCF = OMNI_b36_vcf
    fix_764.reference_sequence = b36_ref
    fix_764.out = new File(vcf_dir+swapExt(OMNI_b36_vcf.getName,".vcf",".ref_fixed.vcf"))
    fix_764_b37.variantVCF = OMNI_b37_vcf
    fix_764_b37.reference_sequence = b37_ref
    fix_764_b37.out = new File(vcf_dir+swapExt(OMNI_b37_vcf.getName,".vcf",".ref_fixed.vcf"))

    add(fix_421,fix_764,fix_764_b37)

    /** Propagate AC/AN annotations to Omni files via variant annotator **/
    var annotate_421 = new VariantAnnotator with OmniArgs
    var annotate_764 = new VariantAnnotator with OmniArgs
    var annotate_764_b37 = new VariantAnnotator with OmniArgs

    annotate_421.variantVCF = fix_421.out
    annotate_421.reference_sequence = b36_ref
    annotate_421.annotation :+= "ChromosomeCounts"
    annotate_421.out = new File(vcf_dir+swapExt(fix_421.out.getName,".vcf",".annot.vcf"))
    annotate_764.variantVCF = fix_764.out
    annotate_764.reference_sequence = b36_ref
    annotate_764.annotation :+= "ChromosomeCounts"
    annotate_764.out = new File(vcf_dir+swapExt(fix_764.out.getName,".vcf",".annot.vcf"))
    annotate_764_b37.variantVCF = fix_764_b37.out
    annotate_764_b37.reference_sequence = b37_ref
    annotate_764_b37.annotation :+= "ChromosomeCounts"
    annotate_764_b37.out = new File(vcf_dir+swapExt(fix_764_b37.out.getName,".vcf",".annot.vcf"))

    add(annotate_421,annotate_764,annotate_764_b37)

    /** Eval the omni chip against the various comps **/
    runEval(annotate_764.out,gunzip_p1_ceu.outFile,"OMNI_764","Pilot1_CEU",pilot1_interval_list, b36_ref)
    runEval(annotate_421.out,gunzip_p1_ceu.outFile,"OMNI_421","Pilot1_CEU",pilot1_interval_list, b36_ref)
    runEval(annotate_764.out,gunzip_p1_chb.outFile,"OMNI_764","Pilot1_CHB",pilot1_interval_list, b36_ref)
    runEval(annotate_421.out,gunzip_p1_chb.outFile,"OMNI_421","Pilot1_CHB",pilot1_interval_list, b36_ref)
    runEval(annotate_764.out,gunzip_p1_yri.outFile,"OMNI_764","Pilot1_YRI",pilot1_interval_list, b36_ref)
    runEval(annotate_421.out,gunzip_p1_yri.outFile,"OMNI_421","Pilot1_YRI",pilot1_interval_list, b36_ref)
    runEval(annotate_764.out,pilot3_release_vcf,"OMNI_764","Pilot3",pilot3_interval_list, b36_ref)
    runEval(annotate_421.out,pilot3_release_vcf,"OMNI_421","Pilot3",pilot3_interval_list, b36_ref)
    runEval(annotate_764_b37.out,gunzip_ag_eur.outFile,"OMNI_764","August_EUR",production_interval_list, b37_ref)
    runEval(annotate_764_b37.out,gunzip_ag_asn.outFile,"OMNI_764","August_ASN",production_interval_list, b37_ref)
    runEval(annotate_764_b37.out,gunzip_ag_afr.outFile,"OMNI_764","Ausust_AFR",production_interval_list, b37_ref)
    runEval(annotate_764.out,gunzip_hiseq.outFile,"OMNI_764","HiSeq",hiseq_interval_list, b36_ref)
    runEval(annotate_764.out,OMNI_hapmap_b36_vcf,"OMNI_764","OMNI_421",pilot1_interval_list,b36_ref)

    var eval1KG_exclude = new VariantEval with OmniArgs
    eval1KG_exclude.samples :+= "/broad/shptmp/chartl/omni/scratch/OMNI_764_vs_Pilot3.sample_overlap.exclude.mixups.txt"
    eval1KG_exclude.rodBind :+= new RodBind("evalOMNI_764","VCF",annotate_764.out)
    eval1KG_exclude.rodBind :+= new RodBind("compPilot3","VCF",pilot3_release_vcf)
    eval1KG_exclude.evalModule :+= "GenotypeConcordance"
    eval1KG_exclude.evalModule :+= "SimpleMetricsBySample"
    eval1KG_exclude.reference_sequence = b36_ref
    eval1KG_exclude.reportType = Some(VE2TemplateType.CSV)
    eval1KG_exclude.intervalsString :+= pilot3_interval_list
    eval1KG_exclude.out = new File(eval_dir+"%s_vs_%s.%s".format("OMNI_764","Pilot3","exclude.mixups.eval.csv"))

    add(eval1KG_exclude)

    runAFComparison(annotate_764.out, gunzip_p1_ceu.outFile, gunzip_p1_chb.outFile, gunzip_p1_yri.outFile)

  }

  def processAuxiliaryChipData(otherChips: File) : List[(String,File)] = {
    // todo ==== me
    return Nil
  }

  def runEval(eval: File, comp: File, eBase: String, cBase: String, intervals: String, reference: File) = {
    var base = "%s_vs_%s".format(eBase,cBase)
    var getOverlap = new GetSampleOverlap
    getOverlap.in_vcfs :+= eval
    getOverlap.in_vcfs :+= comp
    getOverlap.outFile = new File(scratch_dir+base+".sample_overlap.txt")
    add(getOverlap)

    var vEval = new VariantEval with OmniArgs
    vEval.samples :+= getOverlap.outFile.getAbsolutePath
    vEval.rodBind :+= new RodBind("eval"+eBase,"VCF",eval)
    vEval.rodBind :+= new RodBind("comp"+cBase,"VCF",comp)
    vEval.evalModule :+= "GenotypeConcordance"
    vEval.evalModule :+= "SimpleMetricsBySample"
    vEval.intervalsString :+= intervals
    vEval.reference_sequence = reference
    vEval.reportType = Some(VE2TemplateType.CSV)

    vEval.out = new File(eval_dir+base+".eval.csv")

    add(vEval)

  }

  def swapExt(s: String, d: String, f: String) : String = {
    return s.stripSuffix(d)+f
  }

  def runAFComparison(omni: File, p1ceu: File, p1asn: File, p1yri:File ) : Boolean = {
    // step one, set up some of the info
    var populations : List[String] = Nil // these are the pilot 1 populations
    populations :+= "CEU"
    populations :+= "CHBJPT"
    populations :+= "YRI"
    var panels : List[String] = Nil // these are the analysis panels
    panels :+= "EUR"
    panels :+= "ASN"
    panels :+= "ASW"
    panels :+= "AFR"
    panels :+= "ADM"
    // step two -- subset the OMNI chip to the actual sample names
    var nameToSubset: HashMap[String,SelectVariants] = new HashMap[String,SelectVariants]
    for ( p <- populations ) {
      nameToSubset += p -> sampleSubset(p,omni)
    }

    for ( p <- panels ) {
      nameToSubset += p -> sampleSubset(p,omni)
    }

    // step three -- compare the pilot 1 files against all populations and panels

    runComps("Pilot1CEU",p1ceu,"CEU",nameToSubset("CEU").out)
    runComps("Pilot1CEU",p1ceu,"EUR",nameToSubset("EUR").out)
    runComps("Pilot1CHBJPT",p1asn,"CHBJPT",nameToSubset("CHBJPT").out)
    runComps("Pilot1CHBJPT",p1asn,"ASN",nameToSubset("ASN").out)
    runComps("Pilot1YRI",p1yri,"YRI",nameToSubset("YRI").out)
    runComps("Pilot1YRI",p1yri,"AFR",nameToSubset("AFR").out)
    runComps("EUR",nameToSubset("EUR").out,"AFR",nameToSubset("AFR").out)
    runComps("EUR",nameToSubset("EUR").out,"ASN",nameToSubset("ASN").out)
    runComps("EUR",nameToSubset("EUR").out,"ASW",nameToSubset("ASW").out)
    runComps("EUR",nameToSubset("EUR").out,"AMR",nameToSubset("ADM").out)

    return true

  }

  def getOmniSampleListByPanel(panel: String) : String = {
    return scratch_dir+"OMNI_764_%s.txt".format(panel)
  }

  def sampleSubset(panel: String, omni: File) : SelectVariants = {
    var sv : SelectVariants = new SelectVariants with OmniArgs
    sv.reference_sequence = b36_ref
    sv.variantVCF = omni
    sv.sample :+= getOmniSampleListByPanel(panel)
    sv.out = new File(vcf_dir+swapExt(omni.getName,".vcf",".%s.vcf".format(panel)))
    add(sv)
    return sv
  }

  def runComps(eBase: String, evalVCF: File, cBase: String, compVCF: File) = {
    var eval: VariantEval = new VariantEval with OmniArgs
    eval.reference_sequence = b36_ref
    eval.rodBind :+= new RodBind("eval%s".format(eBase),"VCF",evalVCF)
    eval.rodBind :+= new RodBind("comp%s".format(cBase),"VCF",compVCF)
    eval.noStandard = true
    eval.E :+= "AlleleFrequencyComparison"

    add(eval)

    var combine: CombineVariants = new CombineVariants with OmniArgs
    combine.reference_sequence = b36_ref
    combine.rodBind :+= new RodBind(eBase,"VCF",evalVCF)
    combine.rodBind :+= new RodBind(cBase,"VCF",compVCF)
    combine.out = new File(vcf_dir+"%s_plus_%s.vcf".format(eBase,cBase))
    combine.genotypeMergeOptions = Some(VariantContextUtils.GenotypeMergeType.UNIQUIFY)
    combine.priority = "%s,%s".format(eBase,cBase)

    add(combine)

  }
}
