import collection.JavaConversions
import org.broadinstitute.sting.queue.function.JarCommandLineFunction
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.IOUtils
import org.broadinstitute.sting.utils.text.XReadLines

class MultiFullCallingPipeline extends QScript {
  qscript =>

  @Input(doc="Sting home", shortName="stingHome")
  var stingHome: File = _

  @Input(doc="yaml lists to run", shortName="YL")
  var yamlList: File = _

  @Argument(doc="number of jobs per batch", shortName="BS")
  var batchSize: Int = _

  @Argument(doc="pipeline status to", shortName="PS", required = false)
  var pipelineStatusTo: String = _

  @Argument(doc="pipeline job queue", shortName="PJQ", required = false)
  var pipelineJobQueue: String = _

  @Argument(doc="pipeline short queue", shortName="PSQ", required = false)
  var pipelineShortQueue: String = _

  @Argument(doc="pipeline priority", shortName="PP", required = false)
  var pipelinePriority: Option[Int] = None

  @Argument(doc="pipeline retry", shortName="PR", required = false)
  var pipelineRetry: Option[Int] = None

  @Argument(doc="run with -tearScript", shortName="TS")
  var runWithTearScript = false

  def script {
    // Global arguments for all pipeline runs
    stingHome = IOUtils.absolute(stingHome)
    val queueJar = new File(stingHome, "dist/Queue.jar")
    val pipelineScript = new File(stingHome, "scala/qscript/playground/FullCallingPipeline.q")
    val gatkJar = new File(stingHome, "dist/GenomeAnalysisTK.jar")
    val tearScript = if (runWithTearScript) new File(stingHome, "R/DataProcessingReport/GetTearsheetStats.R") else null

    // Parse the yaml list
    var yamls = List.empty[File]
    for (yaml <- JavaConversions.asScalaIterator(new XReadLines(yamlList)))
      yamls :+= new File(yaml)

    // The list of previous outputs
    val lastOuts = new Array[File](batchSize)
    for (yamlGroup <- yamls.grouped(batchSize)) {
      for ((yaml, i) <- yamlGroup.zipWithIndex) {
        // Get the last output for index(i), which is null for the first job.
        val lastOut = lastOuts(i)

        // Run the pipeline on the yaml waiting for the last output.
        val runPipeline = new RunPipeline(yaml, lastOut)

        // Add this run to the graph.
        add(runPipeline)

        // Have the next job at index(i) wait for this output file.
        lastOuts(i) = runPipeline.jobOutputFile
      }
    }

    /**
     * Runs a yaml in a pipeline only after a previous pipeline
     * run has produced the passed in output file.
     */
    class RunPipeline(yamlFile: File, lastOutput: File) extends JarCommandLineFunction {
      private var yamlName = yamlFile.getName.stripSuffix(".yaml")

      @Input(doc="output file to wait for", required=false)
      var waitJobOutputFile = lastOutput

      @Output(doc="virtual output file tagging this pipeline as complete")
      var pipelineComplete = new File(yamlFile.getParentFile, yamlName + ".mfcp")

      commandDirectory = yamlFile.getParentFile
      jobOutputFile = IOUtils.absolute(commandDirectory, yamlName + ".queue.txt")
      jarFile = queueJar
      memoryLimit = Some(1)

      override def commandLine = super.commandLine +
        optional(" -statusTo ", qscript.pipelineStatusTo) +
        optional(" -jobQueue ", qscript.pipelineJobQueue) +
        optional(" -shortJobQueue ", qscript.pipelineShortQueue) +
        optional(" -jobPriority ", qscript.pipelinePriority) +
        optional(" -retry ", qscript.pipelineRetry) +
        optional(" -tearScript ", tearScript) +
        " -S %s --gatkjar %s -jobProject %s -jobPrefix %s -Y %s -bsub -run"
          .format(pipelineScript, gatkJar, yamlName, yamlName, yamlFile)

      override def dotString = "Queue: " + yamlName
    }
  }
}
