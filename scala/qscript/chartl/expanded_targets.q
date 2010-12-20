import org.broadinstitute.sting.commandline.ArgumentCollection
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.library.ipf.ExpandIntervals
import org.broadinstitute.sting.queue.pipeline.PipelineArgumentCollection
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._

class expanded_targets extends QScript {
  @ArgumentCollection var args : PipelineArgumentCollection = new PipelineArgumentCollection
  @Argument(shortName="bait",doc="The list of baits associated with the target list",required=false) var baitFile : File = _
  @Argument(shortName="thisTrigger",doc="The trigger track to use",required=false) var thisTrigger : File = new File("/humgen/gsa-hphome1/chartl/projects/exome/expanded/triggers/joined.omni.hiseq.vcf")

  def script = {

    val intervalExpands : List[ExpandIntervals] = (new Range(0,40,1)).toList.map( u => {
      new ExpandIntervals(args.projectIntervals,1+5*u,5,new File(System.getProperty("user.dir")+"/"+args.projectName+"_expanded_%d_%d.interval_list".format(1+5*u,6+5*u)),args.projectRef,"INTERVALS")
    })

    addAll(intervalExpands)

    val callFiles: List[File] = intervalExpands.map(_.outList).map(makeCalls(_,20))

  }

  def makeCalls(iList: File, scatterNo: Int): List[File] = {
    var scatters : List[(File,File)] = _

    var eval : VariantEval = new VariantEval
    eval.rodBind :+= new RodBind("evalInterval","vcf",filter.out)
    eval.rodBind :+= new RodBind("compHiSeq","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/NA12878/NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.vcf"))
    eval.rodBind :+= new RodBind("compHiSeq_atSites","vcf",callHiseq.out)
    eval.rodBind :+= new RodBind("compOMNI","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni_2.5_764_samples.b37.deduped.annot.vcf"))
    eval.out = swapExt(iList,".interval_list",".eval")
    eval.reportType = Option(org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType.CSV)

  }

  def makeCalls(iList: File) : File = {

    var bams : List[File] = asScalaIterable(new XReadLines(args.projectBams)).toList.map(u => new File(u))

    trait GATKArgs extends CommandLineGATK {
      this.reference_sequence = args.projectRef
      this.DBSNP = args.projectDBSNP
      this.intervals = List(iList)
      this.jarFile = args.gatkJar
      this.memoryLimit = Some(6)
    }

    var rtc : RealignerTargetCreator = new RealignerTargetCreator with GATKArgs
    rtc.out = swapExt(iList,".interval_list",".targets.txt")
    rtc.input_file = bams
    var clean : IndelRealigner = new IndelRealigner with GATKArgs
    clean.targetIntervals = rtc.out
    clean.out = swapExt(iList,".interval_list",".cleaned.bam")
    clean.input_file = bams
    clean.maxReads = Some(100000)
    var call : UnifiedGenotyper = new UnifiedGenotyper with GATKArgs
    call.scatterCount = 4
    call.input_file = List(clean.out)
    call.out = swapExt(iList,".interval_list",".raw.vcf")
    call.trig_emit_conf = Some(0.0)
    call.rodBind :+= new RodBind("trigger","vcf",thisTrigger)
    var filter : VariantFiltration = new VariantFiltration with GATKArgs
    filter.rodBind :+= new RodBind("variant","vcf",call.out)
    filter.filterExpression :+= "\"QD<5.0\""
    filter.filterName :+= "LowQualByDepth"
    filter.filterExpression :+= "\"SB>-0.10\""
    filter.filterName :+= "HighStrandBias"
    filter.out = swapExt(iList,".interval_list",".filtered.vcf")
    var callHiseq : UnifiedGenotyper = new UnifiedGenotyper with GATKArgs
    callHiseq.input_file = List(new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"))
    callHiseq.rodBind :+= new RodBind("trigger","vcf",filter.out)
    callHiseq.out = swapExt(iList,".interval_list",".hiSeq.genotypes.vcf")
    callHiseq.trig_emit_conf = Some(0.0)

    add(rtc,clean,call,filter,callHiseq)

    return (filter.out,callHiSeq.out)

  }
}
