package org.broadinstitute.sting.queue.pipeline

import java.io.File

class PipelineTestSpec {

  /** The arguments to pass to the Queue test, ex: "-S scala/qscript/examples/HelloWorld.scala" */
  var args: String = _

  /** Expected MD5 results for each file path. */
  var fileMD5s = Map.empty[File, String]

  /** VariantEval validations to run on a VCF after the pipeline has completed. */
  var evalSpec: PipelineTestEvalSpec = _

  /** Expected exception from the test. */
  var expectedException: Class[_ <: Exception] = null

  def this(args: String, fileMD5s: Traversable[(File, String)]) = {
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
