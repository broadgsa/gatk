/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.util

import java.net.InetAddress
import java.io.File
import io.Source

/**
 * A collection of various system utilities.
 */
object SystemUtils extends Logging {
  private val localAddress = {
    try {
      InetAddress.getLocalHost
    } catch {
      case e: Exception =>
        InetAddress.getLoopbackAddress
    }
  }
  val inetAddress = localAddress.getHostAddress
  val canonicalHostName = localAddress.getCanonicalHostName

  val hostName = {
    if (canonicalHostName != inetAddress)
      canonicalHostName
    else
      localAddress.getHostName
  }

  val mailName = {
    val mailnameFile = new File("/etc/mailname")
    if (mailnameFile.exists)
      try {
        Source.fromFile(mailnameFile).mkString.trim
      } catch {
        case e: Exception =>
          logger.error("Unabled to read mail domain. Using hostname.", e)
          hostName.split('.').takeRight(2).mkString(".")
      }
    else
      hostName.split('.').takeRight(2).mkString(".")
  }
}
