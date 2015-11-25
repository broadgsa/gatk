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

package org.broadinstitute.gatk.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.gatk.queue.pipeline.{QueueTest, QueueTestSpec}

class HelloWorldQueueTest {
  @Test(timeOut=36000000)
  def testHelloWorld() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorld"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala"
    spec.jobRunners = QueueTest.allJobRunners
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldWithRunName() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithRunName"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -runName HelloWorld"
    spec.jobRunners = QueueTest.allJobRunners
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldWithMemoryLimit() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldMemoryLimit"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -memLimit 1.25"
    spec.jobRunners = QueueTest.allJobRunners
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldWithPriority() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithPriority"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -jobPriority 100"
    spec.jobRunners = QueueTest.allJobRunners
    QueueTest.executeTest(spec)
  }

  @Test(enabled=false, timeOut=36000000)
  def testHelloWorldWithLsfResource() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithLsfResource"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -jobResReq rusage[iodine_io=1] -jobResReq select[swp>0] -jobResReq order[swp]"
    spec.jobRunners = Seq("Lsf706")
    QueueTest.executeTest(spec)
  }

  @Test(enabled=false, timeOut=36000000)
  def testHelloWorldWithLsfResourceAndMemoryLimit() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithLsfResourceAndMemoryLimit"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -memLimit 1.25 -jobResReq rusage[iodine_io=1] -jobResReq select[swp>0] -jobResReq order[swp]"
    spec.jobRunners = Seq("Lsf706")
    QueueTest.executeTest(spec)
  }

  @Test(enabled=false, timeOut=36000000)
  def testHelloWorldWithLsfEnvironment() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithLsfEnvironment"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -jobEnv tv"
    spec.jobRunners = Seq("Lsf706")
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldWithGridEngineResource() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithGridEngineResource"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -jobResReq s_core=1000M"
    spec.jobRunners = Seq("GridEngine")
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldWithGridEngineResourceAndMemoryLimit() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithGridEngineResourceAndMemoryLimit"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -memLimit 1.25 -jobResReq s_core=1000M"
    spec.jobRunners = Seq("GridEngine")
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldWithGridEngineEnvironment() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithGridEngineEnvironment"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -jobEnv \"make 1\""
    spec.jobRunners = Seq("GridEngine")
    QueueTest.executeTest(spec)
  }

  // disabled because our DRMAA implementation doesn't support wallTime
  @Test(enabled=false, timeOut=36000000)
  def testHelloWorldWithWalltime() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithWalltime"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -wallTime 100"
    spec.jobRunners = QueueTest.allJobRunners
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldWithLogDirectory() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithLogDirectory"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -logDir pipelineLogDir"
    spec.jobRunners = QueueTest.allJobRunners
    spec.expectedFilePaths = Seq("pipelineLogDir/HelloWorld-1.out")
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldParallelShell() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithLogDirectory"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala"
    spec.jobRunners = Seq("ParallelShell")
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testHelloWorldParallelShellMaxConcurrentRun() {
    val spec = new QueueTestSpec
    spec.name = "HelloWorldWithLogDirectory"
    spec.args = "-S " + QueueTest.publicQScriptsPackageDir + "examples/HelloWorld.scala" +
      " -maxConcurrentRun 10"
    spec.jobRunners = Seq("ParallelShell")
    QueueTest.executeTest(spec)
  }
}
