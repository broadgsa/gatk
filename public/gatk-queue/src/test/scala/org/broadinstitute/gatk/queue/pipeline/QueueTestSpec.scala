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

package org.broadinstitute.gatk.queue.pipeline

class QueueTestSpec(var name: String = null) {

  /** The arguments to pass to the Queue test, ex: "-S scala/qscript/examples/HelloWorld.scala" */
  var args: String = _

  /** Job Queue to run the test.  Default is null which means use hour. */
  var jobQueue: String = _

  /** Job runners to run the test.  Default is null which means use the default. */
  var jobRunners: Seq[String] = _

  /** Expected MD5 results for each file path. */
  var fileMD5s = Map.empty[String, String]

  /** VariantEval validations to run on a VCF after the pipeline has completed. */
  var evalSpec: QueueTestEvalSpec = _

  /** Expected exception from the test. */
  var expectedException: Class[_ <: Exception] = null

  /** Expected files. The file name may contain wildcards acceptable by the WildcardFileFilter. */
  var expectedFilePaths: Seq[String] = Seq.empty

  /** Unexpected files. The file name may contain wildcards acceptable by the WildcardFileFilter. */
  var unexpectedFilePaths: Seq[String] = Seq.empty

  /** If true will check the MD5s without failing. */
  var parameterize = false

  def this(args: String, fileMD5s: Traversable[(String, String)]) = {
    this()
    this.args = args
    this.fileMD5s = fileMD5s.toMap
  }

  def this(args: String, expectedException: Class[_ <: Exception]) = {
    this()
    this.args = args
    this.expectedException = expectedException
  }
}
