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

package org.broadinstitute.gatk.queue.extensions.picard

import org.broadinstitute.gatk.utils.commandline._

import picard.sam.ValidateSamFile.Mode

import java.io.File

/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 6/22/11
 * Time: 10:35 AM
 */
class ValidateSamFile extends org.broadinstitute.gatk.queue.function.JavaCommandLineFunction with PicardBamFunction {
  analysisName = "ValidateSamFile"
  javaMainClass = "picard.sam.ValidateSamFile"

  @Input(doc="The input SAM or BAM files to analyze.  Must be coordinate sorted.", shortName = "input", fullName = "input_bam_files", required = true)
  var input: Seq[File] = Nil

  @Output(doc="Send output to a file instead of stdout", shortName = "output", fullName = "output_file", required = false)
  var output: File = _

  @Argument(doc="Mode of output", shortName = "mode", fullName = "mode_of_output", required = false)
  var MODE: Mode = Mode.VERBOSE

  @Argument(doc="List of validation error types to ignore.", shortName = "ignore", fullName = "ignore_error_types", required = false)
  var IGNORE: Seq[String] = Nil

  @Argument(doc = "The maximum number of lines output in verbose mode.", shortName = "max", fullName = "max_output", required = false)
  var MAX_OUTPUT: Int = 100

  @Argument(doc = "Reference sequence file, the NM tag check will be skipped if this is missing.", shortName = "ref", fullName ="reference_sequence", required=false)
  var REFERENCE_SEQUENCE: File = _

  @Argument(doc = "If true, only report errors, and ignore warnings.", shortName = "iw", fullName = "ignore_warnings", required = false)
  var IGNORE_WARNINGS: Boolean = false

  @Argument(doc = "If true and input is a BAM file with an index file, also validates the index.", shortName = "vi", fullName = "validate_index", required = false)
  var VALIDATE_INDEX: Boolean = true

  @Argument(doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads. If so, C->T is not counted as an error in computing the value of the NM tag.", shortName = "bs", fullName = "is_bisulfite_sequenced", required = false)
  var IS_BISULFITE_SEQUENCED: Boolean = false

  @Argument(doc = "Relevant for a coordinate-sorted file containing read pairs only. Maximum number of file handles to keep open when spilling mate info to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system.", shortName = "max_files", fullName = "max_open_temp_files", required = false)
  var MAX_OPEN_TEMP_FILES: Int = 8000

  this.sortOrder = null
  override def inputBams = input
  override def outputBam = output
  override def commandLine = super.commandLine +
                             required("MODE=" + MODE) +
                             required("MAX_OUTPUT=" + MAX_OUTPUT) +
                             required("MAX_OPEN_TEMP_FILES=" + MAX_OPEN_TEMP_FILES) +
                             conditional(!VALIDATE_INDEX, "VALIDATE_INDEX=false") +
                             conditional(IGNORE_WARNINGS, "IGNORE_WARNINGS=true") +
                             conditional(IS_BISULFITE_SEQUENCED, "IS_BISULFITE_SEQUENCED=true") +
                             repeat("IGNORE=", IGNORE, spaceSeparated=false)  // repeat() returns "" for null/empty list
}