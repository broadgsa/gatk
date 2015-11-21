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

package org.broadinstitute.gatk.queue.extensions.gatk

import org.broadinstitute.gatk.queue.pipeline.{QueueTest, QueueTestSpec}
import org.broadinstitute.gatk.utils.BaseTest
import org.testng.annotations.Test

class QueueFeaturesQueueTest {

  @Test(timeOut=36000000)
  def testIncludeUnmapped(): Unit = {

    //First case: When no intervals are specified, unmapped reads should be included
    val testOut = "withunmapped.bam"
    val spec = new QueueTestSpec
    spec.name = "includeUnmapped"
    spec.args = Array(
      " -S " + QueueTest.publicQScriptsPackageDir + "examples/ExamplePrintReads.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM_with_unmapped.bam",
      " -out " + testOut).mkString
    spec.fileMD5s += testOut -> "3134a6c732d7f235373095586bc7d470"
    QueueTest.executeTest(spec)

    //Second case: When intervals are explicitly provided, unmapped reads should not be included
    val testOut2 = "withoutunmapped.bam"
    val spec2 = new QueueTestSpec
    spec2.name = "excludeUnmapped"
    spec2.args = Array(
      " -S " + QueueTest.publicQScriptsPackageDir + "examples/ExamplePrintReads.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM_with_unmapped.bam",
      " -L chr1",
      " -out " + testOut2).mkString
    spec2.fileMD5s += testOut2 -> "aa33e589879c4baf6a470d22da76d885"
    QueueTest.executeTest(spec2)
  }

}
