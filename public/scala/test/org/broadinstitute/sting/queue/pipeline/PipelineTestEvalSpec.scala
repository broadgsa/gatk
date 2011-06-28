package org.broadinstitute.sting.queue.pipeline

/**
 * Data validations to evaluate on a GATKReport.
 */
class PipelineTestEvalSpec {
  /** List of eval modules to output. */
  var evalReport: String = _

  /** Validations to assert. */
  var validations: List[PipelineValidation[_]] = Nil
}

/** A VariantEval JEXL and range of values to validate. */
abstract class PipelineValidation[T <: AnyVal](val table: String, val key: String, val metric: String, val target: T, val min: T, val max: T) {
  def parse(x: String): T
  def inRange(x: String): Boolean
}

/** A VariantEval JEXL and target to validate within a 1% tolerance. */
class IntegerValidation(table: String, key: String, metric: String, target: Int)
        extends PipelineValidation[Int](table, key, metric, target,
          (target * .99).floor.toInt, (target * 1.01).ceil.toInt) {
  def parse(x: String) = x.toInt
  def inRange(x: String) = parse(x) >= min && parse(x) <= max
}

/** A VariantEval JEXL and target to validate within a 1% tolerance. */
class DoubleValidation(table: String, key: String, metric: String, target: Double)
        extends PipelineValidation(table, key, metric, target,
          (target * 99).floor / 100, (target * 101).ceil / 100) {
  def parse(x: String) = x.toDouble
  def inRange(x: String) = parse(x) >= min && parse(x) <= max
}
