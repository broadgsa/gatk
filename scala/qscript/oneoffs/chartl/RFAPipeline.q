import org.broadinstitute.sting.commandline.ArgumentCollection
import org.broadinstitute.sting.gatk.CommandLineGATK
import org.broadinstitute.sting.oneoffprojects.walkers.newassociation.RFAArgumentCollection
import org.broadinstitute.sting.queue.extensions.gatk.{RodBind, RFCombine, RFExtractor}
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.IOUtils

class RFAPipeline extends QScript {
  rfapipeline =>

  @Argument(doc="bam list",shortName="bams") var bamList: File = _
  @Argument(doc="sting dir",shortName="sting") var stingDir: String = _
  @Argument(doc="reference file",shortName="ref") var ref: File = _
  @Argument(doc="output base",shortName="o") var baseOut: String = _
  @Argument(doc="interval list",shortName="L") var intervals: File = _
  @ArgumentCollection var rfaArgs : RFAArgumentCollection = new RFAArgumentCollection()
  @Argument(doc="Number of bams per text file output",shortName="br",required = false) var br = 20

  def script = {
    // step one, break the bam files up
    var subBams : Iterator[List[File]] = extractFileEntries(bamList).toList.grouped(br)

    trait ExtractorArgs extends RFExtractor {
      this.reference_sequence = rfapipeline.ref
      this.jarFile = new File(rfapipeline.stingDir+"/dist/GenomeAnalysisTK.jar")
      this.intervals :+= rfapipeline.intervals
      // copy the args into the extractor
      this.windowJump = Some(rfaArgs.windowJump)
      this.windowSize = Some(rfaArgs.windowSize)
      this.fixedZ = Some(rfaArgs.fixedZ)
      this.perSampleZ = Some(rfaArgs.sampleZThresh)
      this.HighInsertSize= Some(rfaArgs.highInsertSize)
      this.LowInsertSize = Some(rfaArgs.lowInsertSize)
      this.clippedBases = Some(rfaArgs.clippedBases)
      this.sampleEpsilon = Some(rfaArgs.EPSILON)
      this.memoryLimit = Some(2)
    }

    val extract : List[RFExtractor] = subBams.zipWithIndex.map( u => {
      var g = new RFExtractor with ExtractorArgs
      g.input_file ++= u._1
      g.out = new File("%s.%d.txt".format(rfapipeline.baseOut,u._2))
      g
    }).toList

    addAll(extract)

    trait CombineArgs extends RFCombine {
      this.reference_sequence = rfapipeline.ref
      this.jarFile = new File(rfapipeline.stingDir+"/dist/GenomeAnalysisTK.jar")
    }

    var combine : RFCombine = new RFCombine with CombineArgs
    var idx : Int = 0
    extract.foreach( ex => {
      val name = "s%d".format(idx)
      val exRB = new RodBind(name,"table",ex.out)
      combine.rodBind :+= exRB
      combine.memoryLimit = Some(6)
      idx+=1;
    })

    combine.out = new File(baseOut)

    add(combine)
  }
}