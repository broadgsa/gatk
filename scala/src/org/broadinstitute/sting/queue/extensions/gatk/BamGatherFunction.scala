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
import org.broadinstitute.sting.queue.extensions.picard.PicardBamFunction
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor

/**
 * Merges BAM files using net.sf.picard.sam.MergeSamFiles.
 */
class BamGatherFunction extends GatherFunction with PicardBamFunction {
  this.javaMainClass = "net.sf.picard.sam.MergeSamFiles"
  this.assumeSorted = Some(true)
  protected def inputBams = gatherParts
  protected def outputBam = originalOutput

  override def freezeFieldValues {
    val originalGATK = originalFunction.asInstanceOf[CommandLineGATK]

    // Whatever the original function can handle, merging *should* do less.
    this.memoryLimit = originalFunction.memoryLimit

    // bam_compression and index_output_bam_on_the_fly from SAMFileWriterArgumentTypeDescriptor
    // are added by the GATKExtensionsGenerator to the subclass of CommandLineGATK

    val compression = QFunction.findField(originalFunction.getClass, SAMFileWriterArgumentTypeDescriptor.COMPRESSION_FULLNAME)
    this.compressionLevel = originalGATK.getFieldValue(compression).asInstanceOf[Option[Int]]

    val disableIndex = QFunction.findField(originalFunction.getClass, SAMFileWriterArgumentTypeDescriptor.DISABLE_INDEXING_FULLNAME)
    this.createIndex = Some(!originalGATK.getFieldValue(disableIndex).asInstanceOf[Boolean])

    super.freezeFieldValues
  }
}
