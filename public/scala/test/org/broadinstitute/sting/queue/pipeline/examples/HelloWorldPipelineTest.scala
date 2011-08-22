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

package org.broadinstitute.sting.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}

class HelloWorldPipelineTest {
  @Test
  def testHelloWorld() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorld"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala"
    spec.jobRunners = PipelineTest.allJobRunners
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithPrefix() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithPrefix"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -jobPrefix HelloWorld"
    spec.jobRunners = PipelineTest.allJobRunners
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithMemoryLimit() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldMemoryLimit"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -memLimit 1.25"
    spec.jobRunners = PipelineTest.allJobRunners
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithPriority() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithPriority"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -jobPriority 100"
    spec.jobRunners = PipelineTest.allJobRunners
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithLsfResource() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithLsfResource"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -jobResReq rusage[iodine_io=1] -jobResReq select[swp>0] -jobResReq order[swp]"
    spec.jobRunners = List("Lsf706")
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithLsfResourceAndMemoryLimit() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithLsfResourceAndMemoryLimit"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -memLimit 1.25 -jobResReq rusage[iodine_io=1] -jobResReq select[swp>0] -jobResReq order[swp]"
    spec.jobRunners = List("Lsf706")
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithLsfEnvironment() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithLsfEnvironment"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -jobEnv tv"
    spec.jobRunners = List("Lsf706")
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithGridEngineResource() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithGridEngineResource"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -jobResReq s_core=1000M"
    spec.jobRunners = List("GridEngine")
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithGridEngineResourceAndMemoryLimit() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithGridEngineResourceAndMemoryLimit"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -memLimit 1.25 -jobResReq s_core=1000M"
    spec.jobRunners = List("GridEngine")
    PipelineTest.executeTest(spec)
  }

  @Test
  def testHelloWorldWithGridEngineEnvironment() {
    val spec = new PipelineTestSpec
    spec.name = "HelloWorldWithGridEngineEnvironment"
    spec.args = "-S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/HelloWorld.scala" +
      " -jobEnv \"make 1\""
    spec.jobRunners = List("GridEngine")
    PipelineTest.executeTest(spec)
  }
}
