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

import org.testng.annotations.Test
import org.broadinstitute.gatk.queue.pipeline.{QueueTest, QueueTestSpec}
import org.broadinstitute.gatk.utils.BaseTest

class ExamplePrintReadsQueueTest {
  @Test(timeOut=36000000)
  def testDevNullOutput() {
    val spec = new QueueTestSpec
    spec.name = "devnulloutput"
    spec.args = Array(
      " -S " + QueueTest.publicQScriptsPackageDir + "examples/ExamplePrintReads.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM.bam",
      " -out /dev/null").mkString
    spec.jobRunners = QueueTest.allJobRunners
    QueueTest.executeTest(spec)
  }

  @Test(timeOut=36000000)
  def testCleanupBai() {
    val spec = new QueueTestSpec
    spec.name = "cleanupbai"
    spec.args = Array(
      " -S " + QueueTest.publicQScriptsPackageDir + "examples/ExamplePrintReads.scala",
      " -R " + BaseTest.publicTestDir + "exampleFASTA.fasta",
      " -I " + BaseTest.publicTestDir + "exampleBAM.bam",
      " -out exampleOut.bam").mkString
    spec.jobRunners = QueueTest.allJobRunners
    spec.unexpectedFilePaths :+= ".queue/scatterGather/ExamplePrintReads-1-sg/temp_1_of_1/exampleOut.bai"
    QueueTest.executeTest(spec)
  }
}
