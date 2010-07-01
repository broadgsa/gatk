import org.broadinstitute.sting.queue.QScript._
// Other imports can be added here

setArgs(args)

val gList = List(30)
val sList = List(0.0001, 0.01)
val dList = List(0.0001, 1000.0)
val bList = List(1.0, 1.3)

for (g: Int <- gList) {
  for (s: Double <- sList) {
    for (d: Double <- dList) {
      for(b: Double <- bList) {

        // Using classes defined below
        val gvc = new GenerateVariantClusters
        val vr = new VariantRecalibrator

        gvc.maxGaussians = g
        gvc.shrinkage = s
        gvc.dirichlet = d
        gvc.clusterFile = new File("g%d_s%.6f_d%.6f_b%.2f.cluster".format(g,s,d,b))
        gvc.jobOutputFile = swapExt(gvc.clusterFile, ".cluster", ".gvc.out")

        vr.clusterFile = gvc.clusterFile
        vr.jobOutputFile = swapExt(vr.clusterFile, ".cluster", ".vr.out")
        vr.backOff = b

        add(gvc, vr)
      }
    }
  }
}

// Populate parameters passed in via -P
setParams

// Run the pipeline
run



// A very basic GATK UnifiedGenotyper
class GenerateVariantClusters extends GatkFunction {
  var maxGaussians: Int = _
  var shrinkage: Double = _
  var dirichlet: Double = _

  @Output
  var clusterFile: File = _

  def commandLine = gatkCommandLine("GenerateVariantClusters") +
    "-B input20,VCF,/broad/shptmp/rpoplin/CEUTSI.chr20.filtered.vcf " +
    "-l INFO -L 20 -an QD -an SB -an HaplotypeScore -an HRun " +
    "-resources /humgen/gsa-scr1/rpoplin/sting_dev_vb/R/ " +
    "-mG %d ".format(maxGaussians) +
    "-shrinkage %.6f ".format(shrinkage) +
    "-dirichlet %.6f ".format(dirichlet) +
    "-clusterFile %s".format(clusterFile)
}

// A basic GATK VariantFiltration
class VariantRecalibrator extends GatkFunction {
  var backOff: Double = _

  @Input
  var clusterFile: File = _

  def commandLine = gatkCommandLine("VariantRecalibrator") +
    "-B input20,VCF,/broad/shptmp/rpoplin/CEUTSI.chr20.filtered.vcf " +
    "-l INFO -L 20 -titv 2.1 " +
    "--ignore_filter HARD_TO_VALIDATE " +
    "-resources /humgen/gsa-scr1/rpoplin/sting_dev_vb/R/ " +
    "-backOff %.2f ".format(backOff) +
    "-clusterFile %s ".format(clusterFile) +
    "-output %s".format(clusterFile)
}
