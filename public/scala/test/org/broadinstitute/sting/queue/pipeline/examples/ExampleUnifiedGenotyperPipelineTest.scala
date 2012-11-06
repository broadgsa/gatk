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

import org.testng.annotations.{DataProvider, Test}
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}
import org.broadinstitute.sting.BaseTest

class ExampleUnifiedGenotyperPipelineTest {
  @Test(timeOut=36000000)
  def testUnifiedGenotyper() {
    val spec = new PipelineTestSpec
    spec.name = "unifiedgenotyper"
    spec.args = Array(
      " -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM.bam",
      " -filter QD",
      " -filterExpression 'QD < 2.0'").mkString
    spec.jobRunners = PipelineTest.allJobRunners
    PipelineTest.executeTest(spec)
  }

  @DataProvider(name = "ugIntervals")
  def getUnifiedGenotyperIntervals =
    Array(
      Array("gatk_intervals", BaseTest.validationDataLocation + "intervalTest.intervals"),
      Array("bed_intervals", BaseTest.validationDataLocation + "intervalTest.bed"),
      Array("vcf_intervals", BaseTest.validationDataLocation + "intervalTest.1.vcf")
    ).asInstanceOf[Array[Array[Object]]]

  @Test(dataProvider = "ugIntervals", timeOut=36000000)
  def testUnifiedGenotyperWithIntervals(intervalsName: String, intervalsPath: String) {
    val spec = new PipelineTestSpec
    spec.name = "unifiedgenotyper_with_" + intervalsName
    spec.args = Array(
      " -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala",
      " -I " + BaseTest.validationDataLocation + "OV-0930.normal.chunk.bam",
      " -R " + BaseTest.hg18Reference,
      " -L " + intervalsPath).mkString
    spec.jobRunners = Seq("Lsf706")
    PipelineTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testUnifiedGenotyperNoGCOpt() {
    val spec = new PipelineTestSpec
    spec.name = "unifiedgenotyper_no_gc_opt"
    spec.args = Array(
      " -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM.bam",
      " -noGCOpt").mkString
    spec.jobRunners = PipelineTest.allJobRunners
    PipelineTest.executeTest(spec)
  }

  @DataProvider(name="resMemReqParams")
  def getResMemReqParam = Array(Array("mem_free"), Array("virtual_free")).asInstanceOf[Array[Array[Object]]]

  @Test(dataProvider = "resMemReqParams", timeOut=36000000)
  def testUnifiedGenotyperResMemReqParam(reqParam: String) {
    val spec = new PipelineTestSpec
    spec.name = "unifiedgenotyper_" + reqParam
    spec.args = Array(
      " -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM.bam",
      " -resMemReqParam " + reqParam).mkString
    spec.jobRunners = Seq("GridEngine")
    PipelineTest.executeTest(spec)
  }
}
