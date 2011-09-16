package org.broadinstitute.sting.queue.pipeline

class PipelineTestSpec(var name: String = null) {

  /** The arguments to pass to the Queue test, ex: "-S scala/qscript/examples/HelloWorld.scala" */
  var args: String = _

  /** Job Queue to run the test.  Default is null which means use hour. */
  var jobQueue: String = _

  /** Job runners to run the test.  Default is null which means use the default. */
  var jobRunners: List[String] = _

  /** Expected MD5 results for each file path. */
  var fileMD5s = Map.empty[String, String]

  /** VariantEval validations to run on a VCF after the pipeline has completed. */
  var evalSpec: PipelineTestEvalSpec = _

  /** Expected exception from the test. */
  var expectedException: Class[_ <: Exception] = null

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
