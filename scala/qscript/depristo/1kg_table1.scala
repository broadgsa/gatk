import org.broadinstitute.sting.queue.QScript._
// Other imports can be added here

val UNIVERSAL_GATK_ARGS = " -l INFO " // -L 1 
val unusedArgs = setArgs(args)

class Target(project: String, snpVCF: String, indelVCF: String, calledGenome: Double, targetGenome: Double, pop: String, pilot : String, bam: String = null) {
    def reportFile: String = List(pop, pilot, "report").mkString(".")
    def extraArgs = { 
      val basic = "--project %s --snps %s --calledGenome %f --totalGenome %f --pop %s".format(project, snpVCF, calledGenome, targetGenome, pop)
      basic + (if ( indelVCF == null ) "" else " --indels " + indelVCF)
    }
    
    def getPilot = pilot
    def getProject = project
    def getPop = pop
    def getSNPVCF = snpVCF
    def getIndelVCF = indelVCF
    def hasIndelVCF = indelVCF != null 
    def getBAM = bam
    def hasBAM = bam != null
    def getDOC = List(getPilot, getPop, getProject, "doc").mkString(".")
    def getDOCSummaryFile = "doc/" + getDOC + ".sample_summary"
    def hasDOC = hasBAM
    private def getEval(t: String) = List(getPilot, getPop, getProject, t, "eval").mkString(".")
    def getSNPEval = getEval("snps")
    def getIndelEval = getEval("indels")
}

var targets: List[Target] = List()

val p1Targets = List(("CEU", 2.43e9), ("YRI", 2.39e9), ("CHBJPT", 2.41e9))

for ( (pop: String,called) <- p1Targets ) 
  targets ::= new Target("SRP000031", pop + ".pilot1.vcf", "v1/dindel-v2/"+pop+".low_coverage.2010_06.indel.genotypes.vcf", called, 2.85e9, pop, "pilot1")

// pilot 2
val p2Targets = List(("CEU", 2.264e9), ("YRI", 2.214e9))
for ( (pop: String, called) <- p2Targets )
  targets ::= new Target("SRP000032", "/humgen/gsa-hpprojects/1kg/releases/pilot_paper_calls/trio/snps/" + pop + ".trio.2010_03.genotypes.vcf.gz", "v1/dindel-v2/"+pop+".trio.2010_06.indel.genotypes.vcf", called, 2.85e9, pop, "pilot2")

// pilot 3
for (POP <- List("CEU", "CHB", "CHD", "JPT", "LWK", "TSI", "YRI")) {
  val indels = if ( POP != "LWK" ) "/humgen/gsa-hpprojects/1kg/releases/pilot_paper_calls/exon/indel/"+POP+".exon.2010_06.genotypes.vcf.gz" else null
  targets ::= new Target("SRP000033", "/humgen/gsa-hpprojects/1kg/releases/pilot_paper_calls/exon/snps/" + POP + ".exon.2010_03.genotypes.vcf.gz", indels, 1.43e6, 1.43e6, POP, "pilot3", "/humgen/gsa-hpprojects/1kg/1kg_pilot3/useTheseBamsForAnalysis/pilot3.%s.cleaned.bam".format(POP))
}

// merged files
targets ::= new Target("SRP000031", "pilot1.snps.merged.vcf", "pilot1.indels.merged.vcf", 2.42e9, 2.85e9, "all", "pilot1.merged")
targets ::= new Target("SRP000032", "pilot2.snps.merged.vcf", "pilot2.indels.merged.vcf", 2.565e9, 2.85e9, "all", "pilot2.merged")
targets ::= new Target("SRP000033", "pilot3.snps.merged.vcf", "pilot3.indels.merged.vcf", 1.43e7, 1.43e7, "all", "pilot3.merged")
targets ::= new Target("SRP00003.", "1kg.snps.merged.vcf", "1kg.indels.merged.vcf", 2.42e7, 2.85e9, "all", "1kg.merged")

val INTERVALS = Map(
    "pilot1" -> null,
    "pilot2" -> null,
    "pilot3" -> "/humgen/gsa-hpprojects/1kg/1kg_pilot3/documents/CenterSpecificTargetLists/results/p3overlap.targets.b36.interval_list"
    )

def setupStage(stage: String) = stage match { 
    case "ALL" =>
        // initial pilot1 merge -- autosomes + x
        for ( (pop: String,called) <- p1Targets ) {
            val auto = "/humgen/gsa-hpprojects/1kg/releases/pilot_paper_calls/low_coverage/snps/"+ pop +".low_coverage.2010_07.genotypes.vcf.gz"
            // todo -- remove fixed when Laura gives us the official calls
            val x = "/humgen/gsa-hpprojects/1kg/releases/pilot_paper_calls/low_coverage/snps/"+ pop +".low_coverage.2010_07.xchr.fixed.genotypes.vcf.gz"
            val combineSNPs = new Combine(List(auto, x), pop + ".pilot1.vcf")
            add(combineSNPs)
        }

        // create pilot wide merges
        val pilots = List("pilot2", "pilot1", "pilot3") // order of perference in merging
        for ( pilot <- pilots ) {
            val pilotTargets = targets filter (_.getPilot == pilot)
            val combineSNPs = new Combine(pilotTargets.map(_.getSNPVCF), pilot + ".snps.merged.vcf")
            add(combineSNPs)
    
            if ( pilotTargets(0).getIndelVCF != null ) {
                val combineIndels = new Combine(pilotTargets.map(_.getIndelVCF).filter((x: String) => x != null), pilot + ".indels.merged.vcf")
                add(combineIndels)
            }
        }

        // create project wide merges
        val snps = "1kg.snps.merged.vcf"
        val indels = "1kg.indels.merged.vcf"

        add(new Combine(pilots.map(_ + ".snps.merged.vcf"), snps))
        add(new Combine(pilots.map(_ + ".indels.merged.vcf"), indels))

        // VariantEval of the SNPs
        for (target <- targets) {
          add(new VariantEval(target.getSNPVCF, target.getSNPEval))
          add(new StatPop(target))
        }

    case "DOC" => 
        for (target <- targets) {
          if ( target.hasBAM ) 
            add(new DepthOfCoverage(target.getBAM, target.getDOC, INTERVALS(target.getPilot)))
        }
    case "MASK" => 
        for ( pop <- List("CEU", "YRI", "CHBJPT") ) 
            add(new MaskStats(pop))

    case _ => throw new Exception("Unknown stage" + stage)
}

setupStage(unusedArgs(0))

// Populate parameters passed in via -P
setParams

// Run the pipeline
run

// Using scala anonymous classes
class VariantEval(vcfIn: String, evalOut: String, vcfType: String = "VCF") extends GatkFunction {
    @Input(doc="foo") var vcfFile: File = new File(vcfIn)
    @Output(doc="foo") var evalFile: File = new File(evalOut)
    override def dotString = "VariantEval: " + vcfFile.getName
    def commandLine = gatkCommandLine("VariantEval") + UNIVERSAL_GATK_ARGS + "-D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_b36.rod -reportType Grep -B eval,%s,%s -o %s -E CompOverlap".format(vcfType, vcfFile, evalFile)
}

class StatPop(target: Target) extends CommandLineFunction {
    @Input(doc="foo") var snpVCF = new File(target.getSNPVCF)
    @Input(doc="foo") var snpEval = new File(target.getSNPEval)
    @Input(doc="foo") var indelVCF = if (target.hasIndelVCF) new File(target.getIndelVCF) else {}
    @Output(doc="foo") var reportFile: File = new File(target.reportFile)
    override def dotString = "1kgStats: " + reportFile
    def commandLine = "python ~/dev/GenomeAnalysisTK/trunk/python/1kgStatsForCalls.py -v -a pilot_data.alignment.index -s pilot_data.sequence.index -r /broad/1KG/DCC/ftp/ -o " + target.reportFile + " " + target.extraArgs + (if (target.hasDOC) " -c " + target.getDOCSummaryFile else "") + " --snpsEval " + target.getSNPEval + (if (target.hasIndelVCF) " --indels " + target.getIndelVCF else "")
}

class Combine(vcfsInArg: List[String], vcfOutPath: String) extends GatkFunction {
  @Input(doc="foo") var vcfs = vcfsInArg.map((x: String) => new File(x))
  @Output(doc="foo") var vcfFile: File = new File(vcfOutPath)
  override def dotString = "CombineVariants: " + vcfs.map(_.getName).mkString(",") + " => " + vcfFile.getName
  def commandLine = gatkCommandLine("CombineVariants") + UNIVERSAL_GATK_ARGS + "-variantMergeOptions UNION -genotypeMergeOptions PRIORITIZE -o %s %s -priority %s".format(vcfFile, vcfs.map( input => " -B %s,VCF,%s".format(input.getName,input)).mkString(""), vcfs.map( _.getName ).mkString(","))
}

class MaskStats(pop: String) extends CommandLineFunction {
    @Output(doc="foo") var outFile: File = new File(pop + ".stats")
    def commandLine = "python ~/dev/GenomeAnalysisTK/trunk/python/maskStats.py masks/" + pop + ".mask.fa.gz -x MT -x Y -o " + outFile
}

class DepthOfCoverage(bam: String, docOutPath: String, interval: String) extends GatkFunction {
  @Input(doc="foo") var bamFile: File = new File(bam)
  @Output(doc="foo") var docFile: File = new File(docOutPath)
  override def dotString = "DOC: " + bamFile.getName
  def commandLine = gatkCommandLine("DepthOfCoverage") + UNIVERSAL_GATK_ARGS + "-omitIntervals -omitBaseOutput -mbq 0 -mmq 0 -o %s -I %s".format(docFile, bamFile) + (if (interval != null) " -XL MT -XL Y -L " + interval else "") 
}
