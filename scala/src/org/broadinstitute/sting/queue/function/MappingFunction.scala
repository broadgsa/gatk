package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.queue.engine.QGraph
import scala.collection.immutable.ListMap

/**
 * Utility class to map a set of inputs to set of outputs.
 * The QGraph uses this function internally to return
 */
class MappingFunction(private val in: ListMap[String, Any], private val out: ListMap[String, Any]) extends QFunction {
  def inputs = in
  def outputs = out
  def run(qGraph: QGraph) = null
}
