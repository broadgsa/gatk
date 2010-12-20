import org.broadinstitute.sting.commandline.ArgumentCollection
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.library.ipf.ExpandIntervals
import org.broadinstitute.sting.queue.library.ipf.intervals.GroupIntervals
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
      new ExpandIntervals(args.projectIntervals,1+5*u,5,new File(System.getProperty("user.dir")+"/"+args.projectName+"_expanded_%d_%d.interval_list".format(1+5*u,6+5*u)),args.projectRef,"TSV","INTERVALS")
    })

     trait GATKArgs extends CommandLineGATK {
      this.reference_sequence = args.projectRef
      this.DBSNP = args.projectDBSNP
      this.jarFile = args.gatkJar
    }

    val userDir = System.getProperty("user.dir")

    addAll(intervalExpands)

    val cleanIntervals : ExpandIntervals = new ExpandIntervals(args.projectIntervals,1,210,new File(userDir+"/"+args.projectName+"_expanded_full.interval_list"),args.projectRef,"TSV","INTERVALS")

    add(cleanIntervals)

    val uncleanBams : List[File] = asScalaIterable(new XReadLines(args.projectBams)).toList.map(u => new File(u))
    val realign : List[RealignerTargetCreator] = uncleanBams.map(u => {
      var rtc : RealignerTargetCreator = new RealignerTargetCreator with GATKArgs
      rtc.out = swapExt(userDir,u,".bam",".expanded.targets.txt")
      rtc.input_file :+= u.getAbsoluteFile
      rtc.intervals :+= cleanIntervals.outList
      rtc
    })
    val clean : List[IndelRealigner] = realign.map( u => {
      var clean : IndelRealigner = new IndelRealigner with GATKArgs
      clean.targetIntervals = u.out
      clean.input_file = u.input_file
      clean.memoryLimit = Some(4)
      clean.out = new File(userDir+"/"+swapExt(u.out,".bam",".expanded.targets.bam").getName)
      clean
    })

    addAll(realign)
    addAll(clean)

    val callFiles: List[File] = intervalExpands.map(u => makeCalls(u.outList,20,clean.map(h => h.out)))

  }

  def makeCalls(iList: File, scatterNo: Int, bams: List[File]): File = {
    var scatters : GroupIntervals = new GroupIntervals(iList,20,true,Some(System.getProperty("user.dir")))
    var filtered : List[(File,File)] = scatters.outputList.zipWithIndex.map(v => reallyMakeCalls(v._1,bams,v._2))
    var gatherStandard : VcfGatherFunction = new VcfGatherFunction
    gatherStandard.gatherParts = filtered.map(u => u._1)
    gatherStandard.originalOutput = swapExt(iList,".interval_list",".filtered.calls.vcf")
    var gatherHiseq : VcfGatherFunction = new VcfGatherFunction
    gatherHiseq.gatherParts = filtered.map(u => u._2)
    gatherHiseq.originalOutput = swapExt(iList,".interval_list",".filtered.hiseq.vcf")

    add(scatters)
    add(gatherStandard)
    add(gatherHiseq)

    trait GATKArgs extends CommandLineGATK {
      this.reference_sequence = args.projectRef
      this.DBSNP = args.projectDBSNP
      this.jarFile = args.gatkJar
    }


    var eval : VariantEval = new VariantEval with GATKArgs
    eval.rodBind :+= new RodBind("evalInterval","vcf",gatherStandard.originalOutput)
    eval.rodBind :+= new RodBind("compHiSeq","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/NA12878/NA12878.hg19.HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.optimized.cut.vcf"))
    eval.rodBind :+= new RodBind("compHiSeq_atSites","vcf",gatherHiseq.originalOutput)
    eval.rodBind :+= new RodBind("compOMNI","vcf",new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni_2.5_764_samples.b37.deduped.annot.vcf"))
    eval.out = swapExt(iList,".interval_list",".eval")
    eval.reportType = Option(org.broadinstitute.sting.utils.report.VE2ReportFactory.VE2TemplateType.CSV)

    eval.out
  }

  def reallyMakeCalls(iList: File, bams : List[File], scatterNo: Int) : (File,File) = {

    trait GATKArgs extends CommandLineGATK {
      this.reference_sequence = args.projectRef
      this.DBSNP = args.projectDBSNP
      this.jarFile = args.gatkJar
      this.intervals :+= iList
    }

    var call : UnifiedGenotyper = new UnifiedGenotyper with GATKArgs
    call.input_file = bams
    call.out = swapExt(iList,".interval_list",".scatter%d.raw.vcf".format(scatterNo))
    call.trig_emit_conf = Some(0.0)
    call.rodBind :+= new RodBind("trigger","vcf",thisTrigger)
    var filter : VariantFiltration = new VariantFiltration with GATKArgs
    filter.rodBind :+= new RodBind("variant","vcf",call.out)
    filter.filterExpression :+= "\"QD<5.0\""
    filter.filterName :+= "LowQualByDepth"
    filter.filterExpression :+= "\"SB>-0.10\""
    filter.filterName :+= "HighStrandBias"
    filter.out = swapExt(iList,".interval_list",".scatter%d.filtered.vcf".format(scatterNo))
    var callHiseq : UnifiedGenotyper = new UnifiedGenotyper with GATKArgs
    callHiseq.input_file = List(new File("/humgen/1kg/analysis/bamsForDataProcessingPapers/NA12878.HiSeq.WGS.bwa.cleaned.recal.bam"))
    callHiseq.rodBind :+= new RodBind("trigger","vcf",filter.out)
    callHiseq.out = swapExt(iList,".interval_list",".scatter%d.hiSeq.genotypes.vcf".format(scatterNo))
    callHiseq.trig_emit_conf = Some(0.0)

    add(call,filter,callHiseq)

    return (filter.out,callHiseq.out)

  }
}
