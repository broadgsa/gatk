package org.broadinstitute.sting.queue.engine.scheduling

import org.apache.commons.lang.builder.{HashCodeBuilder, EqualsBuilder}

class ResourceNode() {
  override def equals(p1: Any) = EqualsBuilder.reflectionEquals(this, p1)
  override def hashCode = HashCodeBuilder.reflectionHashCode(this)
  def lookup(key: String) : Option[String] = None
}
