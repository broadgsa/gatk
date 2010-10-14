package org.broadinstitute.sting.queue.util

import java.lang.management.ManagementFactory
import java.net.InetAddress

/**
 * A collection of various system utilites.
 */
object SystemUtils {
  val pidAtHost = ManagementFactory.getRuntimeMXBean.getName
  val hostName = InetAddress.getLocalHost.getCanonicalHostName
  val domainName = hostName.split('.').takeRight(2).mkString(".")
}
