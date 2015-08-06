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

package org.broadinstitute.gatk.queue.engine.pbsengine

import org.broadinstitute.gatk.queue.util.Logging
import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.engine.drmaa.DrmaaJobRunner
import org.ggf.drmaa.Session

/**
 * Runs jobs on a PBS/Torque Engine compute cluster.
 * NB - THIS FILE HAS BEEN MODIFIED from the original code
 * of the GridEngine package
 */
class PbsEngineJobRunner(session: Session, function: CommandLineFunction) extends DrmaaJobRunner(session, function) with Logging {
  // Pbs Engine disallows certain characters from being in job names.
  // This replaces all illegal characters with underscores
  protected override val jobNameFilter = """[\n\t\r/:,@\\*?]"""
  protected override val minRunnerPriority = -1023
  protected override val maxRunnerPriority = 0

  override protected def functionNativeSpec = {
  
  	// create nativeSpec variable
  		var nativeSpec: String = ""

    // If a project name is set specify the project name
    if (function.jobProject != null)
      nativeSpec += " -P " + function.jobProject

    // If the job queue is set specify the job queue
    if (function.jobQueue != null)
      nativeSpec += " -q " + function.jobQueue

    // If the resident set size is requested pass on the memory request
    // mem_free is the standard, but may also be virtual_free or even not available
    
    if (function.qSettings.residentRequestParameter != null && function.residentRequest.isDefined)
      nativeSpec += " -l %s=%dM".format(function.qSettings.residentRequestParameter, function.residentRequest.map(_ * 1024).get.ceil.toInt)

    // If the resident set size limit is defined specify the memory limit
    if (function.residentLimit.isDefined)
      nativeSpec += " -l mem=%dM".format(function.residentLimit.map(_ * 1024).get.ceil.toInt)

    // If more than 1 core is requested, set the proper request
    // the cores will be requested as part of a single node
    
    if ( function.nCoresRequest.getOrElse(1) > 1 ) {
      if ( function.qSettings.dontRequestMultipleCores )
        logger.warn("Sending multicore job %s to farm without requesting appropriate number of cores (%d)".format(
          function.shortDescription, function.nCoresRequest.get))
      else
        nativeSpec += " -l nodes=1:ppn=%d".format(function.nCoresRequest.get)
    }

    // Pass on any job resource requests
    // NB: blank because resource requests in PBS can be preceded by different
    // arguments, i.e. -l but also -o or -j if they are not exactly "resources" strictly speaking
    // therefore the user will add them in the request, i.e. -jobResReq "-j oe"
    // but this will allow more flexibility in setting the options for PBS jobs on different Clusters
    
    nativeSpec += function.jobResourceRequests.map(" " + _).mkString

    // Pass on any job environment names
    nativeSpec += function.jobEnvironmentNames.map(" " + _).mkString

    // If the priority is set specify the priority
    val priority = functionPriority
    if (priority.isDefined)
      nativeSpec += " -p " + priority.get

    logger.debug("Native spec is: %s".format(nativeSpec))
    (nativeSpec + " " + super.functionNativeSpec).trim()
  }
}
