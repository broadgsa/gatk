package org.broadinstitute.sting.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}

class HelloWorldPipelineTest {
  @Test
  def testHelloWorld {
    var testName = "helloworld"
    val spec = new PipelineTestSpec
    spec.args = "-S scala/qscript/examples/HelloWorld.scala -jobPrefix HelloWorld -jobQueue hour"
    PipelineTest.executeTest(testName, spec)
  }
}
