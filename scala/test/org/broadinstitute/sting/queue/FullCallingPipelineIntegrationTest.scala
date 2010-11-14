package org.broadinstitute.sting.queue

import org.testng.annotations.Test
import org.broadinstitute.sting.BaseTest

/**
 * A temporary integration test to ensure that the full calling pipeline compiles.  To be enhanced...
 */
class FullCallingPipelineIntegrationTest extends QScriptTest {
  @Test
  def testCompileFullCallingPipeline = {
    val command = ("-jobProject QueuePipelineTest -S %1$sscala/qscript/fullCallingPipeline.q -Y %2$s " +
            "-refseqTable /humgen/gsa-hpprojects/GATK/data/Annotations/refseq/refGene-big-table-hg18.txt " +
            "--gatkjar %1$sdist/GenomeAnalysisTK.jar -titv 3.0 -skipCleaning").format(
              stingDir, BaseTest.validationDataLocation + "QueuePipelineTestData/QueuePipelineTestData.yaml"
              )
    executeTest("fullCallingPipeline", command, null)
  }
}
