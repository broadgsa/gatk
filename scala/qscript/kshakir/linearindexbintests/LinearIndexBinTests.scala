import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.text._
import org.broadinstitute.sting.queue.extensions.gatk._
import collection.JavaConversions._

class LinearIndexBinTests extends QScript {
  qscript =>

  @Input(doc="The path to the GenomeAnalysisTK.jar file.", shortName="gatk")
  var gatkJar: File = null

  @Argument(doc="Rod list to test.  The first line in the list is the reference.", shortName="BL")
  var rodLists: List[File] = Nil

  @Argument(doc="Number of times to run the test.", shortName="numRuns", required=false)
  var numRuns = 1

  @Input(doc="memory limits", shortName="mem", required=true)
  var memoryLimits: List[Int] = Nil

  @Input(doc="max bin size", shortName="maxBin", required=false)
  var maxBinSize = 512

  def script = {
    var maxFeaturesPerBin = List.empty[String]
    var currMaxFeatures = maxBinSize
    while (currMaxFeatures > 1) {
      maxFeaturesPerBin +:= currMaxFeatures.toString
      currMaxFeatures /= 2
    }
    maxFeaturesPerBin :::= List("0.001", "0.01", "0.1", "1")

    for (run <- 1 to numRuns) {
      for (rodList <- rodLists) {
        val rodListName = rodList.getName
        val lines = new XReadLines(rodList).iterator
        val reference = new File(lines.next)
        val rodFiles = lines.map(path => new File(path)).toList

        for (memoryLimit <- memoryLimits) {
          for (maxFeatures <- maxFeaturesPerBin) {
            val dir = "%s_%smfpb_%02dg_run%02d".format(rodListName, "00000".take(5-maxFeatures.length) + maxFeatures, memoryLimit, run)

            val countRod = new CountRod {
              override def javaOpts = super.javaOpts + " -DMAX_FEATURES_PER_BIN=" + maxFeatures
            }
            countRod.jobOutputFile = new File(dir, "CountRod.out")
            countRod.out = new File(dir, "CountRod.txt")

            countRod.jarFile = qscript.gatkJar
            countRod.reference_sequence = reference
            countRod.memoryLimit = Some(memoryLimit)

            // Some of the BED files don't have a chrM, which makes the GATK angry.  Run unsafe.
            countRod.U = Some(org.broadinstitute.sting.gatk.arguments.ValidationExclusion.TYPE.ALL)

            for ((rodFile, index) <- rodFiles.zipWithIndex) {
               val rodType = rodFile.getName.split("\\.").last
               countRod.rodBind :+= new RodBind(rodType + (index+1), rodType.toUpperCase, rodFile.getAbsolutePath)
            }

            add(countRod)
          }
        }
      }
    }
  }
}
