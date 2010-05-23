package org.broadinstitute.sting.queue.engine.scheduling

import org.apache.commons.lang.builder.{HashCodeBuilder, EqualsBuilder}

class ResourceNode(val resources: Map[String,String]) {
  override def equals(p1: Any) = EqualsBuilder.reflectionEquals(this, p1)
  override def hashCode = HashCodeBuilder.reflectionHashCode(this)
}
