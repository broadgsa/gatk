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

package org.broadinstitute.gatk.queue.extensions.gatk

import org.broadinstitute.gatk.queue.function.scattergather.GatherFunction
import org.broadinstitute.gatk.queue.function.RetryMemoryLimit

/**
 * Merges a vcf text file.
 */
class VcfGatherFunction extends CombineVariants with GatherFunction with RetryMemoryLimit {
  this.assumeIdenticalSamples = true
  this.suppressCommandLineHeader = true

  private lazy val originalGATK = this.originalFunction.asInstanceOf[CommandLineGATK]

  override def freezeFieldValues() {
    this.variant = this.gatherParts.zipWithIndex map { case (input, index) => new TaggedFile(input, "input"+index) }
    this.out = this.originalOutput
    GATKIntervals.copyIntervalArguments(this.originalGATK, this)

    this.no_cmdline_in_header = originalGATK.no_cmdline_in_header
    this.sites_only = originalGATK.sites_only

    // ensure that the gather function receives the same unsafe parameter as the scattered function
    this.unsafe = this.originalGATK.unsafe

    super.freezeFieldValues()
  }
}
