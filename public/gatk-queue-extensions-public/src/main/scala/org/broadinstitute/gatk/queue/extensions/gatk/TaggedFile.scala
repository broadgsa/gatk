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

import java.io.File
import org.broadinstitute.gatk.utils.io.FileExtension
import org.broadinstitute.gatk.queue.util.ShellUtils

/**
 * Used to provide tagged -I input_file arguments to the GATK.
 */
class TaggedFile(path: String, val tag: String) extends File(path) with FileExtension {
  def this(file: File, tag: String) =
    this(file.getPath, tag)
  def withPath(path: String) = new TaggedFile(path, tag)
}

/**
 * Used to provide -I input_file arguments to the GATK.
 */
object TaggedFile {
  def apply(path: String, tag: String) = new TaggedFile(path, tag)
  def apply(file: File, tag: String) = new TaggedFile(file, tag)

  def formatCommandLineParameter(cmdLineParam: String, value: Any) = {
    value match {
      case taggedFile: TaggedFile if (taggedFile.tag != null) =>
        "%s:%s".format(cmdLineParam, taggedFile.tag)
      case file: File =>
        cmdLineParam
      case x =>
        ""
    }
  }
}
