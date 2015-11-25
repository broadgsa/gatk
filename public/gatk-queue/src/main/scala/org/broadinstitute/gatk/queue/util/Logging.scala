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

import org.apache.log4j._
import java.lang.{Throwable, String}
import org.apache.log4j.spi.{ThrowableInformation, LocationInfo, LoggingEvent, LoggerFactory}
import java.util.Hashtable

/**
 * A mixin to add logging to a class
 */
trait Logging {
  private val className = this.getClass.getName
  protected lazy val logger = Logger.getLogger(className, Logging.factory)
}

object Logging {
  private object factory extends LoggerFactory {
    /**
     * Do not log the scala generated class name after "$".
     * All these shenanigans can be avoided if we use the format
     * string c{1} (the name passed to Logger.getLogger)
     * versus C{1} (the class name of the calling method)
     */
    def makeNewLoggerInstance(name: String) = new Logger(name) {
      override def forcedLog(fqcn: String, level: Priority, message: AnyRef, t: Throwable) = {
        // Strip off the "$" from the class name.
        val info = new LocationInfo(new Throwable, fqcn) { override def getClassName = super.getClassName.takeWhile(_ != '$') }
        // If the following aren't retrieved here, LoggingEvent will NOT lazy load them.
        val throwable = if (t == null) null else new ThrowableInformation(t)
        val ndc = NDC.get
        val mdc = MDC.getContext().asInstanceOf[Hashtable[_,_]];
        callAppenders(new LoggingEvent(fqcn, this, System.currentTimeMillis,
          level.asInstanceOf[Level], message, null, throwable, ndc, info, mdc))
      }
    }
  }
}
