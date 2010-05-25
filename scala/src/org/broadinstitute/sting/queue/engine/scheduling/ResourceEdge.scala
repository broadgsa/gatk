package org.broadinstitute.sting.queue.engine.scheduling

abstract class ResourceEdge {
  def traverse(graph: JobScheduler): Unit
  def lookup(key: String) : Option[String] = None
}
