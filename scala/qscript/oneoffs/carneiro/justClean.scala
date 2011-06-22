/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: 3/17/11
 * Time: 11:29 AM
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.queue.extensions.gatk.{IndelRealigner, RealignerTargetCreator, RodBind}
import org.broadinstitute.sting.queue.QScript


class justClean extends QScript {

  @Input(doc="path to GenomeAnalysisTK.jar", shortName="gatk", required=true)
  var GATKjar: File = _

  @Input(doc="input BAM file - or list of BAM files", shortName="i", required=true)
  var input: File = _

  @Input(doc="Reference fasta file", shortName="R", required=false)
  var reference: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(doc="dbsnp ROD to use (VCF)", shortName="D", required=false)
  var dbSNP: File = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_132_b37.leftAligned.vcf")

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", shortName="indels", required=false)
  var indels: File = new File("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf")


  val queueLogDir: String = ".qlog/"


  def script = {

    println(GATKjar)

    val outBam = swapExt(input, ".bam", ".clean.bam")
    val tIntervals = swapExt(input, ".bam", ".all_indels.intervals")

    val target = new RealignerTargetCreator()
    target.input_file :+= input
    target.out = tIntervals
    target.reference_sequence = reference
    target.mismatchFraction = 0.0
    target.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    target.rodBind :+= RodBind("indels", "VCF", indels)
    target.memoryLimit = 6
    target.jobName = queueLogDir + tIntervals + ".atarget"
    target.jarFile = GATKjar
    target.scatterCount = 84
    


    val clean = new IndelRealigner()
    clean.input_file :+= input
    clean.targetIntervals = tIntervals
    clean.out = outBam
    clean.reference_sequence = reference
    clean.rodBind :+= RodBind("dbsnp", "VCF", dbSNP)
    clean.rodBind :+= RodBind("indels", "VCF", indels)
    clean.doNotUseSW = false
    clean.jobName = queueLogDir + outBam + ".clean"
    clean.jarFile = GATKjar
    clean.memoryLimit = 8
    clean.scatterCount = 84

    add(target, clean);
  }
}
