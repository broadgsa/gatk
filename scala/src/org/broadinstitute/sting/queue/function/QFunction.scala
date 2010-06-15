package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.queue.engine.QGraph
import scala.collection.immutable.ListMap

trait QFunction {
  def inputs: ListMap[String, Any]
  def outputs: ListMap[String, Any]
  def missingValues = Set.empty[String]
}
