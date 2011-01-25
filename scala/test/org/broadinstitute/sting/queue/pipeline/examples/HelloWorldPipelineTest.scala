package org.broadinstitute.sting.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}

class HelloWorldPipelineTest {
  @Test
  def testHelloWorld {
    var testName = "helloworld"
    val spec = new PipelineTestSpec
    spec.args = "-S scala/qscript/examples/HelloWorld.scala -jobPrefix HelloWorld -jobQueue hour"
    // TODO: working example of MD5 usage.
    // spec.fileMD5s += new File(PipelineTest.runDir(testName) + "hello.out") -> "0123456789abcdef0123456789abcdef"
    PipelineTest.executeTest(testName, spec)
  }
}
