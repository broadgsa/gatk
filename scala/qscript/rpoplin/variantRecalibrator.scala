import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript

class variantRecalibrator extends QScript {
  @Argument(doc="gatkJarFile")
  var gatkJarFile: File = _
  
  def script = {

val gList = List(30)
val sList = List(0.0001, 0.01)
val dList = List(0.0001, 1000.0)
val bList = List(1.0, 1.3)

for (g: Int <- gList) {
  for (s: Double <- sList) {
    for (d: Double <- dList) {
      for(b: Double <- bList) {

        // Using classes defined by QueueGATKExtensions.jar
        val gvc = new GenerateVariantClusters
        val vr = new VariantRecalibrator

        gvc.jarFile = gatkJarFile
        gvc.rodBind :+= RodBind("input20", "VCF", new File("/broad/shptmp/rpoplin/CEUTSI.chr20.filtered.vcf"))
        gvc.logging_level = "INFO"
        gvc.intervalsString :+= "20"
        gvc.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
        gvc.path_to_resources = "/humgen/gsa-scr1/rpoplin/sting_dev_vb/R/"
        gvc.maxGaussians = Some(g)
        gvc.shrinkage = Some(s)
        gvc.shrinkageFormat = "%.6f"
        gvc.dirichlet = Some(d)
        gvc.dirichletFormat = "%.6f"
        gvc.clusterFile = "g%d_s%.6f_d%.6f_b%.2f.cluster".format(g,s,d,b)
        gvc.jobOutputFile = new File(gvc.clusterFile.stripSuffix(".cluster") + ".gvc.out")

        vr.jarFile = gatkJarFile
        vr.rodBind :+= RodBind("input20", "VCF", new File("/broad/shptmp/rpoplin/CEUTSI.chr20.filtered.vcf"))
        vr.logging_level = "INFO"
        vr.intervalsString :+= "20"
        vr.target_titv = Some(2.1)
        vr.ignore_filter :+= "HARD_TO_VALIDATE"
        vr.path_to_resources = "/humgen/gsa-scr1/rpoplin/sting_dev_vb/R/"
        vr.clusterFile = gvc.clusterFile
        vr.jobOutputFile = new File(vr.clusterFile.stripSuffix(".cluster") + ".vr.out")
        vr.backOff = Some(b)
        vr.backOffFormat = "%.2f"

        add(gvc, vr)
      }
    }
  }
}
  }
}
