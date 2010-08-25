import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * An example building on the intro ExampleCountReads.scala.
 * Runs an INCOMPLETE version of the UnifiedGenotyper with VariantEval and optional VariantFiltration.
 * If run on a compute cluster, splits the UnifiedGenotyper into 3 parts.
 */
class ExampleUnifiedGenotyper extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="The path to the GenomeAnalysisTK.jar file.", shortName="gatk")
  var gatkJar: File = null // The command line must pass the gatk jar to this script via -gatk.

  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _ // _ is scala shorthand for null

  @Input(doc="Bam file to genotype.", shortName="I")
  var bamFile: File = _

  // The following arguments are all optional.

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Argument(doc="A optional list of filter names.", shortName="filter", required=false)
  var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(doc="An optional list of filter expressions.", shortName="filterExpression", required=false)
  var filterExpressions: List[String] = Nil


  // This trait allows us set the variables below in one place,
  // and then reuse this trait on each CommandLineGATK function below.
  trait UnifiedGenotyperArguments extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.referenceFile
    this.intervals = qscript.intervals
    // Some() is how you set the value for an scala Option.
    // Set the memory limit to 2 gigabytes on each command.
    this.memoryLimit = Some(2)
  }


  def script = {
    // Create the four function that we can run.
    val genotyper = new UnifiedGenotyper with UnifiedGenotyperArguments
    val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    val evalFiltered = new VariantEval with UnifiedGenotyperArguments

    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    genotyper.scatterCount = 3
    genotyper.input_file :+= qscript.bamFile.toNamedFile
    genotyper.out = swapExt(qscript.bamFile, "bam", "unfiltered.vcf")

    evalUnfiltered.rodBind :+= RodBind("vcf", "VCF", genotyper.out)
    evalUnfiltered.out = swapExt(genotyper.out, "vcf", "eval")

    variantFilter.rodBind :+= RodBind("vcf", "VCF", genotyper.out)
    variantFilter.out = swapExt(qscript.bamFile, "bam", "filtered.vcf")
    variantFilter.filterName = filterNames
    variantFilter.filterExpression = filterExpressions

    evalFiltered.rodBind :+= RodBind("vcf", "VCF", variantFilter.out)
    evalFiltered.out = swapExt(variantFilter.out, "vcf", "eval")

    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    if (filterNames.size > 0)
      add(variantFilter, evalFiltered)
  }
}
