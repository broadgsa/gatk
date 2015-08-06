/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.pipeline

/**
 * Data validations to evaluate on a GATKReport.
 */
class QueueTestEvalSpec {
  /** Eval modules to output. */
  var evalReport: String = _

  /** Validations to assert. */
  var validations: Seq[PipelineValidation[_]] = Nil
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
