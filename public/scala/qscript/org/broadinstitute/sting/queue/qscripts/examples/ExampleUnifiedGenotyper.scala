package org.broadinstitute.sting.queue.qscripts.examples

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

/**
 * An example building on the intro ExampleCountReads.scala.
 * Runs an INCOMPLETE version of the UnifiedGenotyper with VariantEval and optional VariantFiltration.
 */
class ExampleUnifiedGenotyper extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

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
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    // Set the memory limit to 2 gigabytes on each command.
    this.memoryLimit = 2
  }

  def script() {
    // Create the four function that we can run.
    val genotyper = new UnifiedGenotyper with UnifiedGenotyperArguments
    val variantFilter = new VariantFiltration with UnifiedGenotyperArguments
    val evalUnfiltered = new VariantEval with UnifiedGenotyperArguments
    val evalFiltered = new VariantEval with UnifiedGenotyperArguments

    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    genotyper.scatterCount = 3
    genotyper.input_file :+= qscript.bamFile
    genotyper.out = swapExt(qscript.bamFile, "bam", "unfiltered.vcf")

    evalUnfiltered.eval :+= genotyper.out
    evalUnfiltered.out = swapExt(genotyper.out, "vcf", "eval")

    variantFilter.variant = genotyper.out
    variantFilter.out = swapExt(qscript.bamFile, "bam", "filtered.vcf")
    variantFilter.filterName = filterNames
    variantFilter.filterExpression = filterExpressions.map("\"" + _ + "\"")

    evalFiltered.eval :+= variantFilter.out
    evalFiltered.out = swapExt(variantFilter.out, "vcf", "eval")

    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    if (filterNames.size > 0)
      add(variantFilter, evalFiltered)
  }
}
