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

package org.broadinstitute.gatk.queue.engine.drmaa

import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.engine.CommandLineJobManager
import org.broadinstitute.gatk.utils.jna.drmaa.v1_0.JnaSessionFactory
import org.ggf.drmaa.Session

/**
 * Runs jobs using DRMAA
 */
class DrmaaJobManager extends CommandLineJobManager[DrmaaJobRunner] {
  protected var session: Session = _

  protected def newSession() = new JnaSessionFactory().getSession
  protected def contact = null

  override def init() {
    session = newSession()
    session.init(contact)
  }

  override def exit() {
    session.exit()
  }

  def runnerType = classOf[DrmaaJobRunner]
  def create(function: CommandLineFunction) = new DrmaaJobRunner(session, function)

  override def updateStatus(runners: Set[DrmaaJobRunner]) = {
    var updatedRunners = Set.empty[DrmaaJobRunner]
    runners.foreach(runner => if (runner.updateJobStatus()) {updatedRunners += runner})
    updatedRunners
  }
  override def tryStop(runners: Set[DrmaaJobRunner]) {
    runners.foreach(_.tryStop())
  }
}
