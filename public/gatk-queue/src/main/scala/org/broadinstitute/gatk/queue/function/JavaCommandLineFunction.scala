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

package org.broadinstitute.gatk.queue.function

import org.broadinstitute.gatk.utils.commandline.Argument
import org.broadinstitute.gatk.utils.io.IOUtils
import java.io.File

/**
 * Defines a command line function that runs java code.
 */
trait JavaCommandLineFunction extends CommandLineFunction {
  @Argument(doc="jar", exclusiveOf="javaMainClass")
  var jarFile: File = _

  @Argument(doc="Main class to run from javaClasspath", exclusiveOf="jarFile")
  var javaMainClass: String = _

  /**
   * Class path for the main class.
   * Defaults to the current classpath.
   */
  var javaClasspath: Seq[String] = Nil

  /**
   * Memory limit for the java executable, or if None will use the default memoryLimit.
   */
  @Argument(doc="Java memory limit", required=false)
  var javaMemoryLimit: Option[Double] = None

  /**
   * Max number of GC threads
   */
  var javaGCThreads: Option[Int] = None

  /**
   * Max percent of time spent in garbage collection
   */
  var javaGCTimeLimit: Option[Int] = None

  /**
   * Min percent of max heap freed during a garbage collection
   */
  var javaGCHeapFreeLimit: Option[Int] = None

  override def freezeFieldValues() {
    super.freezeFieldValues()

    if (javaMemoryLimit.isEmpty && memoryLimit.isDefined)
      javaMemoryLimit = memoryLimit

    if (javaMainClass != null && javaClasspath.isEmpty)
      javaClasspath = JavaCommandLineFunction.currentClasspath

    if (!this.qSettings.disableDefaultJavaGCOptimizations) {
      // By default set the GC threads to 4
      if (javaGCThreads.isEmpty)
        javaGCThreads = Some(4)

      // By default exit if more than 50% of time in GC
      if (javaGCTimeLimit.isEmpty)
        javaGCTimeLimit = Some(50)

      // By default exit if GC does not free up 10% of the heap
      if (javaGCHeapFreeLimit.isEmpty)
        javaGCHeapFreeLimit = Some(10)
    }
  }


  override def copySettingsTo(function: QFunction) {
    super.copySettingsTo(function)
    function match {
      case java: JavaCommandLineFunction =>
        if (java.javaMemoryLimit.isEmpty)
          java.javaMemoryLimit = this.javaMemoryLimit
        if (java.javaGCThreads.isEmpty)
          java.javaGCThreads = this.javaGCThreads
        if (java.javaGCTimeLimit.isEmpty)
          java.javaGCTimeLimit = this.javaGCTimeLimit
        if (java.javaGCHeapFreeLimit.isEmpty)
          java.javaGCHeapFreeLimit = this.javaGCHeapFreeLimit
      case _ => /* ignore */
    }
  }

  /**
   * Returns the java executable to run.
   */
  def javaExecutable: String = {
    if (jarFile != null)
      required("-jar", jarFile)
    else if (javaMainClass != null)
      required("-cp", javaClasspath.mkString(File.pathSeparator)) +
      required(javaMainClass)
    else
      null
  }

  def javaOpts = Array(
    optional("-Xmx", javaMemoryLimit.map(gb => (gb * 1024).ceil.toInt), "m", spaceSeparated=false),
    conditional(javaGCThreads.isDefined || javaGCTimeLimit.isDefined || javaGCHeapFreeLimit.isDefined, "-XX:+UseParallelOldGC"),
    optional("-XX:ParallelGCThreads=", javaGCThreads, spaceSeparated=false),
    optional("-XX:GCTimeLimit=", javaGCTimeLimit, spaceSeparated=false),
    optional("-XX:GCHeapFreeLimit=", javaGCHeapFreeLimit, spaceSeparated=false),
    required("-Djava.io.tmpdir=", jobTempDir, spaceSeparated=false)).mkString("")

  def commandLine = required("java") +
                    javaOpts +
                    javaExecutable
}

object JavaCommandLineFunction {
  val currentClasspath = System.getProperty("java.class.path")
    .split(File.pathSeparatorChar).map(path => IOUtils.absolute(new File(path)).getPath).toSeq
}
