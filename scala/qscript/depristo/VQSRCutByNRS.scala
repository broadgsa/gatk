import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.apache.commons.io.FilenameUtils;
import scala.io.Source._

class recalibrate extends QScript {
  // @Input(doc="bamIn", shortName="I", required=true)
  // var bamList: File = _
  
  @Argument(doc="gatk jar file")
  var gatkJarFile: File = _

  // @Argument(shortName = "R", doc="ref")
  // var referenceFile: File = _

  @Argument(fullName = "prefix", doc="Prefix argument", required=false)
  var prefix: String = ""

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK { logging_level = "INFO"; jarFile = gatkJarFile; memoryLimit = Some(4) }

class Target(val name: String, val reference: File, val rodName: String, val VCF: File, val intervals: Option[String], val titvTarget: Double) {
    def clusterFile = new File(name + ".clusters")
}

val hg18 = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
val b36 = new File("/humgen/1kg/reference/human_b36_both.fasta")
val hg19 = new File("/seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19.fasta")

val HiSeq = new Target("NA12878.HiSeq", hg18, "hg18", new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/HiSeq.WGS.cleaned.ug.snpfiltered.indelfiltered.vcf"), None, 2.07)
val WEx = new Target("NA12878.WEx", hg18, "hg18", new File("/home/radon01/depristo/work/oneOffProjects/1000GenomesProcessingPaper/wgs.v13/GA2.WEx.cleaned.ug.snpfiltered.indelfiltered.vcf"), Some("~/localData/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.interval_list"), 2.6)
val LowPassN60 = new Target("lowpass.N60", b36, "b36", new File("lowpass.N60.chr20.filtered.vcf"), Some("20"), 2.3)
val LowPassAugust = new Target("ALL.august.v3", hg19, "b37", new File("ALL.august.v3.chr20.filtered.vcf"), Some("20"), 2.3)
val TGPWExFH = new Target("1000G.WEx.FH", hg19, "b37", new File("/humgen/gsa-pipeline/PQ7LC/all_batches_v006/Plate_1/SnpCalls/Barcoded_1000G_WEx_Plate_1.cleaned.annotated.handfiltered.vcf"), Some("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"), 3.0)
val TGPWExGdA = new Target("1000G.WEx.GdA", hg19, "b37", new File("/humgen/gsa-scr1/delangel/NewUG/calls/AugustRelease.filtered_Q50_QD5.0_SB0.0.allSamples.SNPs_hg19.WEx_UG_newUG_MQC.vcf"), Some("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"), 3.0)

val targets = List(HiSeq, WEx, LowPassN60, LowPassAugust, TGPWExFH, TGPWExGdA)

def script = {
    for (target <- targets) {
      add(new GenerateVariantClusters(target) )
      add(new VariantRecalibratorTiTv(target) )
      add(new VariantRecalibratorNRS(target) )
    }
}

def bai(bam: File) = new File(bam + ".bai")

val FiltersToIgnore = List("DPFilter", "ABFilter", "ESPStandard", "QualByDepth", "StrandBias", "HomopolymerRun")

class GenerateVariantClusters(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.GenerateVariantClusters with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.DBSNP = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_" + t.rodName + ".rod")
    this.rodBind :+= RodBind("hapmap", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr." + t.rodName + "_fwd.vcf")
    this.rodBind :+= RodBind("input", "VCF", t.VCF)
    this.clusterFile = t.clusterFile
    this.use_annotation ++= List("QD", "SB", "HaplotypeScore", "HRun")
    this.analysisName = t.name + "_Cluster"
    if ( t.intervals != None ) this.intervalsString ++= List(t.intervals.get)
    this.qual = Some(300)
    this.std = Some(3.5)
    this.mG = Some(16) // v2 calls
    // ignores
    this.ignoreFilter ++= FiltersToIgnore
}


class VariantRecalibratorBase(t: Target) extends org.broadinstitute.sting.queue.extensions.gatk.VariantRecalibrator with UNIVERSAL_GATK_ARGS {
    this.reference_sequence = t.reference
    this.DBSNP = new File("/humgen/gsa-hpprojects/GATK/data/dbsnp_129_" + t.rodName + ".rod")
    this.rodBind :+= RodBind("hapmap", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr." + t.rodName + "_fwd.vcf")
    this.rodBind :+= RodBind("truth", "VCF", "/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr." + t.rodName + "_fwd.vcf")
    this.rodBind :+= RodBind("input", "VCF", t.VCF)
    this.clusterFile = t.clusterFile
    this.analysisName = t.name+"_VR"
    if ( t.intervals != None ) this.intervalsString ++= List(t.intervals.get)
    this.ignoreFilter ++= FiltersToIgnore
    this.ignoreFilter ++= List("HARD_TO_VALIDATE")
    this.priorDBSNP = Some(2.0)
    this.priorHapMap = Some(2.0)
    this.target_titv = t.titvTarget
}

class VariantRecalibratorTiTv(t: Target) extends VariantRecalibratorBase(t) {
    this.tranche ++= List("0.1", "1.0", "10.0", "100.0")
    this.out = new File(t.name + ".titv.recalibrated.vcf")
    this.tranchesFile = new File(t.name + ".titv.tranches")
}

class VariantRecalibratorNRS(t: Target) extends VariantRecalibratorBase(t) {
    this.sm = Some(org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibrator.SelectionMetricType.TRUTH_SENSITIVITY)
    this.tranche ++= List("50", "25", "10", "5", "2", "1", "0.5", "0.1")
    this.out = new File(t.name + ".ts.recalibrated.vcf")
    this.tranchesFile = new File(t.name + ".ts.tranches")
}
}
