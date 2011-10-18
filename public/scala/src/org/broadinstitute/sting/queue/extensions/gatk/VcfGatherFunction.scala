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

package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.scattergather.GatherFunction
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterArgumentTypeDescriptor

/**
 * Merges a vcf text file.
 */
class VcfGatherFunction extends CombineVariants with GatherFunction {

  private lazy val originalGATK = this.originalFunction.asInstanceOf[CommandLineGATK]

  override def freezeFieldValues {
    this.jarFile = this.originalGATK.jarFile
    this.reference_sequence = this.originalGATK.reference_sequence
    this.intervals = this.originalGATK.intervals
    this.intervalsString = this.originalGATK.intervalsString

    this.variant = this.gatherParts.zipWithIndex map { case (input, index) => new TaggedFile(input, "input"+index) }
    this.out = this.originalOutput
    this.assumeIdenticalSamples = true

    // NO_HEADER and sites_only from VCFWriterArgumentTypeDescriptor
    // are added by the GATKExtensionsGenerator to the subclass of CommandLineGATK

    val noHeader = QFunction.findField(originalFunction.getClass, VCFWriterArgumentTypeDescriptor.NO_HEADER_ARG_NAME)
    this.NO_HEADER = originalGATK.getFieldValue(noHeader).asInstanceOf[Boolean]

    val sitesOnly = QFunction.findField(originalFunction.getClass, VCFWriterArgumentTypeDescriptor.SITES_ONLY_ARG_NAME)
    this.sites_only = originalGATK.getFieldValue(sitesOnly).asInstanceOf[Boolean]

    super.freezeFieldValues
  }
}
