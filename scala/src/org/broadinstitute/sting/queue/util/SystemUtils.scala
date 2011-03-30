/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.util

import java.lang.management.ManagementFactory
import java.net.InetAddress
import java.io.File
import io.Source

/**
 * A collection of various system utilities.
 */
object SystemUtils {
  val inetAddress = InetAddress.getLocalHost.getHostAddress
  val hostName = InetAddress.getLocalHost.getCanonicalHostName

  val mailName = {
    val mailnameFile = new File("/etc/mailname")
    if (mailnameFile.exists)
      Source.fromFile(mailnameFile).mkString.trim
    else if (hostName == inetAddress)
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
