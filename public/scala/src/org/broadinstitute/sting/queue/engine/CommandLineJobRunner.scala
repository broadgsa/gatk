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

package org.broadinstitute.sting.queue.engine

import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.io.File
import org.broadinstitute.sting.queue.util.{Logging, IOUtils}

/**
 * Runs a command line function.
 */
trait CommandLineJobRunner extends JobRunner[CommandLineFunction] with Logging {

  /** A generated exec shell script. */
  protected var jobScript: File = _

  /** Which directory to use for the job status files. */
  protected def jobStatusDir = function.jobTempDir

  override def init() {
    super.init()
    var exec = new StringBuilder
    
    var dirs = Set.empty[File]
    for (dir <- function.jobDirectories)
      dirs += IOUtils.dirLevel(dir, 2)
    if (dirs.size > 0) {
      // prepend "cd '<dir_1>' [&& cd '<dir_n>']" to automount the directories.
      exec.append(dirs.mkString("cd '", "' && cd '", "'"))
      exec.append(" && cd '%s' && \\%n".format(function.commandDirectory))
    }
    exec.append(function.commandLine)

    this.jobScript = IOUtils.writeTempFile(exec.toString, ".exec", "", jobStatusDir)
  }

  override def cleanup() {
    super.cleanup()
    IOUtils.tryDelete(jobScript)
  }
}
