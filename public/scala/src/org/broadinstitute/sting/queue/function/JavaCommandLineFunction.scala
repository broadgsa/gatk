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

package org.broadinstitute.sting.queue.function

import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.utils.io.IOUtils
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
  var javaClasspath: List[String] = Nil

  /**
   * Memory limit for the java executable, or if None will use the default memoryLimit.
   */
  var javaMemoryLimit: Option[Double] = None

  /**
   * Returns the java executable to run.
   */
  def javaExecutable: String = {
    if (jarFile != null)
      "-jar " + jarFile
    else if (javaMainClass != null)
      "-cp \"%s\" %s".format(javaClasspath.mkString(File.pathSeparator), javaMainClass)
    else
      null
  }

  override def freezeFieldValues() {
    super.freezeFieldValues()

    if (javaMemoryLimit.isEmpty && memoryLimit.isDefined)
      javaMemoryLimit = memoryLimit

    if (javaMainClass != null && javaClasspath.isEmpty)
      javaClasspath = JavaCommandLineFunction.currentClasspath
  }

  def javaOpts = "%s -Djava.io.tmpdir=%s"
    .format(optional(" -Xmx", javaMemoryLimit.map(gb => (gb * 1024).ceil.toInt), "m"), jobTempDir)

  def commandLine = "java%s %s"
    .format(javaOpts, javaExecutable)
}

object JavaCommandLineFunction {
  val currentClasspath = System.getProperty("java.class.path")
    .split(File.pathSeparatorChar).map(path => IOUtils.absolute(new File(path)).getPath).toList
}
