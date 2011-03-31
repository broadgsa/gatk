/*
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 3/29/11
 * Time: 5:31 PM
 */
package oneoffs.carneiro;

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class QuickCCTest extends QScript {
  qscript =>

  @Input(doc="path to GenomeAnalysisTK.jar", shortName="gatk", required=true)
  var GATKjar: File = _

  @Input(doc="input BAM file - or list of BAM files", shortName="i", required=true)
  var input: File = _

  @Input(doc="Reference fasta file", shortName="R", required=false)
  var reference: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(doc="dbsnp ROD to use (VCF)", shortName="D", required=false)
  var dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  @Input(shortName="L", required=false)
  var intervals: List[String] = Nil


  val queueLogDir: String = ".qlog/"


  def script = {

    val recal = new File("recal.csv")

    val cc = new CountCovariates()
    cc.reference_sequence = reference
    cc.input_file :+= input
    cc.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    cc.intervalsString = intervals
    cc.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    cc.scatterCount = 4
    cc.recal_file = recal
    cc.memoryLimit = 4

    add(cc);
  }
}