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

package org.broadinstitute.gatk.queue.function.scattergather

import java.io.File
import org.broadinstitute.gatk.utils.commandline.Input
import org.broadinstitute.gatk.queue.function.QFunction

/**
 * Base class for Scatter functions.
 */
trait ScatterFunction extends QFunction {
  analysisName = "Scatter"

  var originalFunction: ScatterGatherableFunction = _

  @Input(doc="Original inputs to scatter")
  var originalInputs: Set[File] = _

  override def shortDescription = analysisName + ": %s ...".format(firstOutput.getName)

  /**
   * Called to initialize scatter function values after all other values have been setup for this function.
   */
  def init() {}

  /**
   * Returns true if the scatter function can scatter this original function.
   * @return true if the scatter function can scatter this original function.
   */
  def isScatterGatherable = true

  /**
   * Returns the number of clones that should be created.
   */
  def scatterCount: Int = originalFunction.scatterCount

  /**
   * Initializes the input fields for the clone function.
   * The input values should be set to their defaults
   * and may be changed by the user.
   * @param cloneFunction CloneFunction to initialize.
   * @param index The one based scatter index.
   */
  def initCloneInputs(cloneFunction: CloneFunction, index: Int) {}

  /**
   * Binds the input fields for the clone function to this scatter function.
   * The input values should be set to their absolute values and added
   * to scatter parts.
   * @param cloneFunction CloneFunction to bind.
   * @param index The one based scatter index.
   */
  def bindCloneInputs(cloneFunction: CloneFunction, index: Int) {}
}
