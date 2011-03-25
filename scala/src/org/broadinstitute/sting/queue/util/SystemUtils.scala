package org.broadinstitute.sting.queue.util

import java.lang.management.ManagementFactory
import java.net.InetAddress

/**
 * A collection of various system utilities.
 */
object SystemUtils {
  val inetAddress = InetAddress.getLocalHost.getHostAddress
  val hostName = InetAddress.getLocalHost.getCanonicalHostName

  val domainName = {
    if (hostName == inetAddress)
      inetAddress
    else
      hostName.split('.').takeRight(2).mkString(".")
  }

  val pidAtHost = {
    val mxBeanName = ManagementFactory.getRuntimeMXBean.getName
    if (hostName == inetAddress)
      mxBeanName
    else
      mxBeanName.split('.').head
  }
}
