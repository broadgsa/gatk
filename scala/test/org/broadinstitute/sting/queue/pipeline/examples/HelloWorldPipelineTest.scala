package org.broadinstitute.sting.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}

class HelloWorldPipelineTest {
  @Test
  def testHelloWorld {
    val spec = new PipelineTestSpec
    spec.name = "helloworld"
    spec.args = "-S scala/qscript/examples/HelloWorld.scala -jobPrefix HelloWorld"
    PipelineTest.executeTest(spec)
  }
}
