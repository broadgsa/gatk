import org.broadinstitute.sting.queue.QScript._
import org.apache.commons.io.FilenameUtils;
// Other imports can be added here

val unusedArgs = setArgs(args)

def runPipeline() = { 
    for (bamIn <- inputs(".bam")) {
      val root = bamIn.getPath()
      val bamRoot = FilenameUtils.removeExtension(root);
      val recalData = bamRoot + ".recal_data.csv"
      val recalBam = bamRoot + ".recal.bam"
      val recalRecalData = bamRoot + ".recal.recal_data.csv"
      //add(new CountCovariates(root, recalData, "-OQ"))
      val tableRecal = new TableRecalibrate(root, recalData, recalBam)
      tableRecal.intervals = new File("/humgen/gsa-hpprojects/GATK/data/chromosomes.hg18.interval_list")
      tableRecal.scatterCount = 25
      add(tableRecal)
      add(new Index(recalBam))
      add(new CountCovariates(recalBam, recalRecalData))
      add(new AnalyzeCovariates(recalData, recalData + ".analyzeCovariates"))
      add(new AnalyzeCovariates(recalRecalData, recalRecalData + ".analyzeCovariates"))
    }
}

runPipeline()

// Populate parameters passed in via -P
setParams

// Run the pipeline
run

class Index(bamIn: String) extends GatkFunction {
    @Input(doc="foo") var bam = bamIn
    memoryLimit = Some(1)
    def commandLine = "samtools index %s".format(bam)
}

class CountCovariates(bamIn: String, recalDataIn: String, args: String = "") extends GatkFunction {
    @Input(doc="foo") var bam = bamIn
    @Output(doc="foo") var recalData = recalDataIn
    memoryLimit = Some(4)
    def commandLine = gatkCommandLine("CountCovariates") + args + " -l INFO -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -I %s --max_reads_at_locus 20000 -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile %s".format(bam, recalData)
}

class TableRecalibrate(bamInArg: String, recalDataIn: String, bamOutArg: String) extends GatkFunction {
    @Input(doc="foo") var bamIn = bamInArg
    @Input(doc="foo") var recalData = recalDataIn
    @Gather(classOf[BamGatherFunction])
    @Output(doc="foo") var bamOut = new File(bamOutArg)
    memoryLimit = Some(2)
    def commandLine = gatkCommandLine("TableRecalibration") + "-l INFO -I %s -recalFile %s -outputBam %s".format(bamIn, recalData, bamOut.getPath())
}

class AnalyzeCovariates(recalDataIn: String, outputDir: String) extends GatkFunction {
    @Input(doc="foo") var recalData = recalDataIn
    memoryLimit = Some(4)
    def commandLine = "java -Xmx4g -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/AnalyzeCovariates.jar -recalFile %s -outputDir %s -resources /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/R/ -ignoreQ 5".format(recalData, outputDir)
}