package org.broadinstitute.sting.queue.engine

import graphing.{ExplicitJobGrapher, TreeJobGrapher}
import scheduling.ResourceNode

/**
 * Syntactic sugar for filling in a pipeline using a Scala script.
 */
object Pipeline {
  def addRule(rule: (Any, Any), commandString: String): Unit = TreeJobGrapher.addRule(rule, commandString)
  def run(args: Array[String], sources: Any, targets: Any): Unit = TreeJobGrapher.run(args, sources, targets)

  def node() = ExplicitJobGrapher.node()
  def addEdge(rule: (ResourceNode, ResourceNode), commandString: String): Unit = ExplicitJobGrapher.addEdge(rule, commandString)
  def run(): Unit = ExplicitJobGrapher.run()
}
