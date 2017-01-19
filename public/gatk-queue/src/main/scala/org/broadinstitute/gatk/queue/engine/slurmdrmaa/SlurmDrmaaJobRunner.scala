/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.queue.engine.slurmdrmaa

import org.broadinstitute.gatk.queue.engine.drmaa.DrmaaJobRunner
import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.util.Logging
import org.ggf.drmaa.Session

/**
 * Runs jobs on a SLURM compute cluster.
 * This requires to install this plugin for SLURM: http://apps.man.poznan.pl/trac/slurm-drmaa
 *
 * @author Peter van 't Hof
 */
class SlurmDrmaaJobRunner(session: Session, function: CommandLineFunction) extends DrmaaJobRunner(session, function) with Logging {
  // Pbs Engine disallows certain characters from being in job names.
  // This replaces all illegal characters with underscores
  protected override val jobNameFilter = """[\n\t\r/:,@\\*?]"""
  protected override val minRunnerPriority = -1023
  protected override val maxRunnerPriority = 0

  override protected def functionNativeSpec = {
  
  	// create nativeSpec variable
  		var nativeSpec: String = ""

    // If the qualityOfSerice is set, specify the qualityOfSerice
    function.qualityOfSerice.foreach(nativeSpec += " --qos=" + _)

    // If the resident set size limit is defined specify the memory limit
    function.residentLimit.map(l => (l * 1024).ceil.toInt).foreach(l => nativeSpec += s" --mem=${l}")

    // If walltime is set specify the walltime
    function.wallTime.foreach(t => nativeSpec += s" --time=${t}:00:00")

    // If more than 1 core is requested, set the proper request
    // the cores will be requested as part of a single node
    if ( function.nCoresRequest.getOrElse(1) > 1 ) {
      if ( function.qSettings.dontRequestMultipleCores )
        logger.warn("Sending multicore job %s to farm without requesting appropriate number of cores (%d)".format(
          function.shortDescription, function.nCoresRequest.get))
      else
        nativeSpec += " --ntasks=%d".format(function.nCoresRequest.getOrElse(1))
    }

    val n = nativeSpec + " " + super.functionNativeSpec
    logger.info(s"Native spec is: $n")
    n.trim()
  }
}
