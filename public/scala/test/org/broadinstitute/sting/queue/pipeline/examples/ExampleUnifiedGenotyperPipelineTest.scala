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
import org.broadinstitute.sting.BaseTest

class ExampleUnifiedGenotyperPipelineTest {
  @Test
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

  @Test
  def testUnifiedGenotyperWithGatkIntervals() {
    val spec = new PipelineTestSpec
    spec.name = "unifiedgenotyper_with_gatk_intervals"
    spec.args = Array(
      " -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala",
      " -I " + BaseTest.validationDataLocation + "OV-0930.normal.chunk.bam",
      " -R " + BaseTest.hg18Reference,
      " -L " + BaseTest.validationDataLocation + "intervalTest.intervals").mkString
    spec.jobRunners = Seq("Lsf706")
    PipelineTest.executeTest(spec)
  }

  @Test
  def testUnifiedGenotyperWithBedIntervals() {
    val spec = new PipelineTestSpec
    spec.name = "unifiedgenotyper_with_bed_intervals"
    spec.args = Array(
      " -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala",
      " -I " + BaseTest.validationDataLocation + "OV-0930.normal.chunk.bam",
      " -R " + BaseTest.hg18Reference,
      " -L " + BaseTest.validationDataLocation + "intervalTest.bed").mkString
    spec.jobRunners = Seq("Lsf706")
    PipelineTest.executeTest(spec)
  }

  @Test
  def testUnifiedGenotyperWithVcfIntervals() {
    val spec = new PipelineTestSpec
    spec.name = "unifiedgenotyper_with_vcf_intervals"
    spec.args = Array(
      " -S public/scala/qscript/org/broadinstitute/sting/queue/qscripts/examples/ExampleUnifiedGenotyper.scala",
      " -I " + BaseTest.validationDataLocation + "OV-0930.normal.chunk.bam",
      " -R " + BaseTest.hg18Reference,
      " -L " + BaseTest.validationDataLocation + "intervalTest.1.vcf").mkString
    spec.jobRunners = Seq("Lsf706")
    PipelineTest.executeTest(spec)
  }
}
