import java.io.File
import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.queue.extensions.gatk.{SelectVariants, RodBind, VariantsToTable}
import org.broadinstitute.sting.queue.QScript

/*
* Created by IntelliJ IDEA.
* User: carneiro
* Date: 4/12/11
* Time: 11:24 AM
*/

class mendelianViolation extends QScript
{

  @Argument(shortName="trio", doc="input trio VCF file", required=true)
  var trio: File = _

  @Argument(shortName="daughter", doc="daughter input VCF file", required=true)
  var daughter: File = _

  @Argument(shortName="family", doc="family string", required=false)
  var family: String = "NA12891+NA12892=NA12878"

  @Argument(shortName="mvq", doc="mendelian violation quality", required=false)
  var mvq: Double = 20

  @Input(doc="path to GenomeAnalysisTK.jar", shortName="gatk", required=false)
  var GATKjar: File = new File("/humgen/gsa-scr1/carneiro/stable/dist/GenomeAnalysisTK.jar")

  def script = {
    val reference = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
    val trioViolations = "trio_violations.vcf"
    val daughterViolations = "daughter_violations.vcf"

    val mv = new SelectVariants()
    mv.rodBind :+= RodBind("variant", "VCF", trio)
    mv.family = family
    mv.reference_sequence = reference
    mv.mvq = mvq
    mv.out = trioViolations
    mv.jarFile = GATKjar
    mv.memoryLimit = 4

    val intersection = new SelectVariants()
    intersection.rodBind :+= RodBind("variant", "VCF", daughter)
    intersection.rodBind :+= RodBind("conc","VCF", trioViolations)
    intersection.reference_sequence = reference
    intersection.conc = "conc"
    intersection.out = daughterViolations
    intersection.jarFile = GATKjar
    intersection.memoryLimit = 4

    add(mv, intersection)
  }
}
