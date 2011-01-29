package org.broadinstitute.sting.queue.pipeline.examples

import org.testng.annotations.Test
import org.broadinstitute.sting.queue.pipeline.{PipelineTest, PipelineTestSpec}
import org.broadinstitute.sting.BaseTest

class ExampleCountLociPipelineTest {
  @Test
  def testCountLoci {
    var testOut = "count.out"
    val spec = new PipelineTestSpec
    spec.name = "countloci"
    spec.args = "-S scala/qscript/examples/ExampleCountLoci.scala -gatk %s -R %s -I %s -o %s".format(
      PipelineTest.currentGATK, BaseTest.hg18Reference, BaseTest.validationDataLocation + "small_bam_for_countloci.bam", testOut
    )
    spec.fileMD5s += testOut -> "67823e4722495eb10a5e4c42c267b3a6"
    PipelineTest.executeTest(spec)
  }
}
