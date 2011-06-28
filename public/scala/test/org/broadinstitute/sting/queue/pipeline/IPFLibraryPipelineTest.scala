package org.broadinstitute.sting.queue.pipeline

import org.testng.annotations.Test
import org.broadinstitute.sting.BaseTest

class IPFLibraryPipelineTest {

  @Test
  def testVCFExtractSites {
    var testOut = "vcfes.vcf"
    val spec = new PipelineTestSpec
    spec.name = "vcfExtractSites"
    spec.args = "-S scala/qscript/oneoffs/QTools.q -T VCFExtractSites -ivcf %s -out %s".format(
      BaseTest.validationDataLocation + "omni_1212.subset.b37.vcf", testOut
    )
    spec.fileMD5s += testOut -> "4f496b8cf90302428a9edda486a337f4"
    PipelineTest.executeTest(spec)
  }

  @Test
  def testVCFExtractSamples {
    var testOut = "vcf.extract.samples.vcf"
    val spec = new PipelineTestSpec
    spec.name = "vcfExtractSamples"
    spec.args = "-S scala/qscript/oneoffs/QTools.q -T VCFExtractSamples -ivcf %s -out %s -sm HG00107,HG00500,NA18501,NA18942".format(
      BaseTest.validationDataLocation + "omni_1212.subset.b37.vcf", testOut
    )

    spec.fileMD5s += testOut -> "180d5a2e7a1fbc5117de6705bde3a7c8"
  }

  @Test
  def testVCFExtractIntervals {
    var testOut = "vcf.extract.intervals.list"
    val spec = new PipelineTestSpec
    spec.name = "vcfExtractIntervals"
    spec.args = "-S scala/qscript/oneoffs/QTools.q -T VCFExtractIntervals -ivcf %s -out %s".format(
      BaseTest.validationDataLocation + "omni_1212.subset.b37.vcf", testOut
    )

    spec.fileMD5s += testOut ->"28589eeb28fac753577c027abf939345"

  }

  @Test
  def testVCFSimpleMerge {
    var testOut = "vcf.simplemerge.vcf"
    val spec = new PipelineTestSpec
    val int1 = BaseTest.validationDataLocation + "omni.subset.interleaved.1.vcf"
    val int2 = BaseTest.validationDataLocation + "omni.subset.interleaved.2.vcf"
    spec.name = "vcfSimpleMerge"
    spec.args = "-S scala/qscript/oneoffs/QTools.q -T VCFSimpleMerge -vcfs %s,%s -out %s -ref %s".format(
      int1,int2,testOut,BaseTest.b37KGReference
    )

    spec.fileMD5s += testOut -> "c59b8de64ba787a2ca3a1734bdebf277"
  }

  @Test
  def testSortByRef {
    var testOut = "vcf.refsorted.vcf"
    val spec = new PipelineTestSpec
    val unsorted = BaseTest.validationDataLocation + "omni.pos_sorted.vcf"
    spec.name = "sortByRef"
    spec.args = "-S scala/qscript/oneoffs/QTools.q -T SortByRef -ivcf %s -out %s -ref %s".format(
      unsorted, testOut, BaseTest.b37KGReference
    )

    spec.fileMD5s += testOut -> "ee09af803bc94987d55d044c2ebbc0b8"
  }
}
