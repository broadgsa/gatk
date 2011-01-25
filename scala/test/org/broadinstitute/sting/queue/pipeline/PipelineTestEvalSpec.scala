package org.broadinstitute.sting.queue.pipeline

import java.io.File

/**
 * Data validations to evaluate on a VCF using VariantEval.
 */
class PipelineTestEvalSpec {
  // TODO: Reuse the Project "YAML" object for reference, intervals, etc.

  /** VCF to eval */
  var vcf: File = _

  /** Reference for the VCF */
  var reference: File = _

  /** Intervals for the VCF */
  var intervals: File = _

  /** DBSNP to use for comparisons, via -B:dbsnp,VCF or -D */
  var dbsnp: File = _

  /** List of eval modules to output. */
  var evalModules = List("CompOverlap", "CountFunctionalClasses", "CountVariants", "SimpleMetricsBySample", "TiTvVariantEvaluator")

  /** Validations to assert. */
  var validations: List[PipelineValidation] = Nil
}

/** A VariantEval JEXL and range of values to validate. */
class PipelineValidation(val metric: String, val min: String, val max: String) {
}

/** A VariantEval JEXL and target to validate within a 1% tolerance. */
class IntegerValidation(metric: String, target: Int)
        extends PipelineValidation(metric,
          (target * .99).floor.toInt.toString, (target * 1.01).ceil.toInt.toString) {
}

/** A VariantEval JEXL and target to validate within a 1% tolerance. */
class DoubleValidation(metric: String, target: Double)
        extends PipelineValidation(metric,
          "%.2f".format((target * 99).floor / 100), "%.2f".format((target * 101).ceil / 100)) {
}
