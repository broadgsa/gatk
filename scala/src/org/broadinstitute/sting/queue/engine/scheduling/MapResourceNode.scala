package org.broadinstitute.sting.queue.engine.scheduling

class MapResourceNode (private val resources: Map[String,String]) extends ResourceNode {
  override def lookup(key: String) : Option[String] = this.resources.get(key)
}
