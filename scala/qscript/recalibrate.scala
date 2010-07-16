import java.io.File
import org.broadinstitute.sting.queue.QScript._
import org.apache.commons.io.FilenameUtils;
// Other imports can be added here

val unusedArgs = setArgs(args)

def runPipeline(arg: String) = { 
    val scatter = arg == "scatter"

    for (bamIn <- inputs(".bam")) {
      val root = bamIn.getPath()
      val bamRoot = FilenameUtils.removeExtension(root);
      val recalData = new File(bamRoot + ".recal_data.csv")
      val recalBam = new File(bamRoot + ".recal.bam")
      val recalRecalData = new File(bamRoot + ".recal.recal_data.csv")
      //add(new CountCovariates(root, recalData, "-OQ"))
      val tableRecal = new TableRecalibrate(bamIn, recalData, recalBam, "-OQ")
      if ( scatter ) {
            tableRecal.intervals = new File("/humgen/gsa-hpprojects/GATK/data/chromosomes.hg18.interval_list")
      	    tableRecal.scatterCount = 25
      }
      add(tableRecal)
      add(new Index(recalBam))
      add(new CountCovariates(recalBam, recalRecalData, "-nt 4"))
      add(new AnalyzeCovariates(recalData, new File(recalData.getPath() + ".analyzeCovariates")))
      add(new AnalyzeCovariates(recalRecalData, new File(recalRecalData.getPath() + ".analyzeCovariates")))
    }
}

runPipeline(unusedArgs(0))

// Populate parameters passed in via -P
setParams

// Run the pipeline
run

def bai(bam: File) = new File(bam + ".bai")

class Index(bamIn: File) extends GatkFunction {
    @Input(doc="foo") var bam = bamIn
    @Output(doc="foo") var bamIndex = bai(bamIn)
    memoryLimit = Some(1)
    override def dotString = "Index: %s".format(bamIn.getName)
    def commandLine = "samtools index %s".format(bam)
}

class CountCovariates(bamIn: File, recalDataIn: File, args: String = "") extends GatkFunction {
    @Input(doc="foo") var bam = bamIn
    @Input(doc="foo") var bamIndex = bai(bamIn)
    @Output(doc="foo") var recalData = recalDataIn
    memoryLimit = Some(4)
    override def dotString = "CountCovariates: %s [args %s]".format(bamIn.getName, args)
    def commandLine = gatkCommandLine("CountCovariates") + args + " -l INFO -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_hg18.rod -I %s --max_reads_at_locus 20000 -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile %s".format(bam, recalData)
}

class TableRecalibrate(bamInArg: File, recalDataIn: File, bamOutArg: File, args: String = "") extends GatkFunction {
    @Input(doc="foo") var bamIn = bamInArg
    @Input(doc="foo") var recalData = recalDataIn
    @Gather(classOf[BamGatherFunction])
    @Output(doc="foo") var bamOut = bamOutArg
    override def dotString = "TableRecalibrate: %s => %s [args %s]".format(bamInArg.getName, bamOutArg.getName, args)
    memoryLimit = Some(2)
    def commandLine = gatkCommandLine("TableRecalibration") + args + " -l INFO -I %s -recalFile %s -outputBam %s".format(bamIn, recalData, bamOut) // bamOut.getPath())
}

class AnalyzeCovariates(recalDataIn: File, outputDir: File) extends GatkFunction {
    @Input(doc="foo") var recalData = recalDataIn
    memoryLimit = Some(4)
    override def dotString = "AnalyzeCovariates: %s".format(recalDataIn.getName)
    def commandLine = "java -Xmx4g -jar /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/dist/AnalyzeCovariates.jar -recalFile %s -outputDir %s -resources /home/radon01/depristo/dev/GenomeAnalysisTK/trunk/R/ -ignoreQ 5 -Rscript /broad/tools/apps/R-2.6.0/bin/Rscript".format(recalData, outputDir)
}
