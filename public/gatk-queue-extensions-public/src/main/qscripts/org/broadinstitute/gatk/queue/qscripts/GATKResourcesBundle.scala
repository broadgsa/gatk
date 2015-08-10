/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException
import org.broadinstitute.gatk.queue.extensions.picard.PicardBamFunction
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction

class GATKResourcesBundle extends QScript {
  // todo -- update to released version when things stabilize
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("dist/GenomeAnalysisTK.jar")

  @Argument(doc="liftOverPerl", required=false)
  var liftOverPerl: File = new File("./public/perl/liftOverVCF.pl")

  @Argument(shortName = "ver", doc="The GIT version of this release", required=true)
  var BUNDLE_VERSION: String = _

  @Argument(shortName = "bundleDir", doc="Path to root where resource files will be placed", required=false)
  val BUNDLE_ROOT = new File("/humgen/gsa-hpprojects/GATK/bundle")

  @Argument(shortName = "downloadDir", doc="Path to root where resource files will be placed for users to download", required=false)
  val DOWNLOAD_ROOT = new File("/humgen/gsa-scr1/pub/bundle")

  @Argument(shortName = "test", doc="", required=false)
  val TEST = false

  @Argument(shortName = "phase2", doc="", required=false)
  val DO_DOWNLOAD = false

  val SITES_EXT: String = "sites"

  def BUNDLE_DIR: File = BUNDLE_ROOT + "/" + BUNDLE_VERSION
  def DOWNLOAD_DIR: File = DOWNLOAD_ROOT + "/" + BUNDLE_VERSION

  // REFERENCES
  class Reference( val name: String, val file: File ) { }
  var hg19: Reference = _
  var b37: Reference = _
  var hg18: Reference = _
  var b36: Reference = _
  var exampleFASTA: Reference = _
  var refs: List[Reference] = _

  class Resource(val file: File, val name: String, val ref: Reference, val useName: Boolean = true, val makeSites: Boolean = true, val makeCallsIfBam: Boolean = true ) {
    def destname(target: Reference): String = {
      if ( useName )
        return name + "." + target.name + "." + getExtension(file)
      else
        return file.getName
    }
  }

  def liftover(in: File, inRef: Reference, out: File, outRef: Reference): CommandLineFunction = {
    //Console.printf("liftover(%s => %s)%n", inRef.name, outRef.name)
    (inRef.name, outRef.name) match {
      case ("b37", "hg19") =>
        return new LiftOverPerl(in, out, new File("public/chainFiles/b37tohg19.chain"), inRef, outRef)
      case ("b37", "hg18") =>
        return new LiftOverPerl(in, out, new File("public/chainFiles/b37tohg18.chain"), inRef, outRef)
      case ("b37", "b36") =>
        return new LiftOverPerl(in, out, new File("public/chainFiles/b37tob36.chain"), inRef, outRef)
      case _ => return null
    }
  }

  def isVCF(file: File) = file.getName.endsWith(".vcf")
  def isBAM(file: File) = file.getName.endsWith(".bam")
  def isOUT(file: File) = file.getName.endsWith(".out")
  def isFASTA(file: File) = file.getName.endsWith(".fasta")
  def isIntervalList(file: File) = file.getName.endsWith(".interval_list")

  var RESOURCES: List[Resource] = Nil
  def addResource(comp: Resource) { RESOURCES = comp :: RESOURCES }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.memoryLimit = 2
  }

  def initializeTestDataFiles() = {
    //
    // Standard evaluation files for indel
    //
    b37 = new Reference("b37", new File("/Users/depristo/Desktop/broadLocal/localData/human_g1k_v37.fasta"))
    hg18 = new Reference("hg18", new File("/Users/depristo/Desktop/broadLocal/localData/Homo_sapiens_assembly18.fasta"))
    exampleFASTA = new Reference("exampleFASTA", new File("public/testdata/exampleFASTA.fasta"))
    refs = List(b37, hg18, exampleFASTA)

    val DATAROOT = "/Users/depristo/Desktop/broadLocal/localData/"
    //addResource(new Resource(DATAROOT + "human_g1k_v37.fasta", "human_g1k_v37.fasta", b37, false))
    addResource(new Resource(DATAROOT + "1000G.snp.validation.b37.vcf", "1000G.snp.validation", b37))
    addResource(new Resource(DATAROOT + "dbsnp_132_b37.vcf", "dbsnp_132", b37, true, false))

    addResource(new Resource(exampleFASTA.file, "exampleFASTA", exampleFASTA, false))
    addResource(new Resource("public/testdata/exampleBAM.bam", "exampleBAM", exampleFASTA, false, false, false))
  }

  def initializeStandardDataFiles() = {
    //
    // references
    //
    hg19 = new Reference("hg19", new File("/humgen/gsa-hpprojects/GATK/data/ucsc.hg19/ucsc.hg19.fasta"))
    b37 = new Reference("b37", new File("/humgen/1kg/reference/human_g1k_v37.fasta"))
    hg18 = new Reference("hg18", new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"))
    b36 = new Reference("b36", new File("/humgen/1kg/reference/human_b36_both.fasta"))
    exampleFASTA = new Reference("exampleFASTA", new File("public/testdata/exampleFASTA.fasta"))
    refs = List(hg19, b37, hg18, b36, exampleFASTA)

    addResource(new Resource(b37.file, "", b37, false))
    addResource(new Resource(b36.file, "", b36, false))
    addResource(new Resource(hg19.file, "", hg19, false))
    addResource(new Resource(hg18.file, "", hg18, false))

    //
    // The b37_decoy reference
    //
    addResource(new Resource("/humgen/1kg/reference/human_g1k_v37_decoy.fasta",
          "IGNORE", b37, false, false))

    //
    // standard VCF files.  Will be lifted to each reference
    //
    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_138_b37.leftAligned.vcf",
      "dbsnp_138", b37, true, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_sites_2141_samples.b37.vcf",
      "1000G_omni2.5", b37, true, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/sites_r27_nr.b37_fwd.vcf",
      "hapmap_3.3", b37, true, false))

    addResource(new Resource("/humgen/1kg/DCC/ftp/technical/working/20120312_phase1_v2_indel_cleaned_sites_list/ALL.wgs.phase1_release_v2.20101123.official_indel_calls.20120312.sites.vcf",
      "1000G_phase1.indels", b37, true, false))

    addResource(new Resource("/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.highQuality.vcf",
      "1000G_phase1.snps.high_confidence", b37, true, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Unvalidated/GoldStandardIndel/gold.standard.indel.MillsAnd1000G.b37.vcf",
      "Mills_and_1000G_gold_standard.indels", b37, true, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Affymetrix_Axiom/Axiom_Exome_Plus.sites_only.all_populations.poly.vcf",
      "Axiom_Exome_Plus.sites_only.all_populations.poly", b37, true, false))

    //
    // CEU trio (NA12878,NA12891,NA12892) best practices results
    //

    addResource(new Resource("/humgen/1kg/processing/production_wgs_final/trio/CEU/CEU.wgs.HaplotypeCaller.20131118.snps_indels.high_coverage_pcr_free.genotypes.vcf",
      "CEUTrio.HiSeq.WGS.b37.bestPractices",b37,true,false))

    //
    // NA12878 knowledgebase snapshot
    //

    addResource(new Resource("/humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/snapshots/NA12878.wgs.broad_truth_set.20131119.snps_and_indels.genotypes.vcf",
      "NA12878.knowledgebase.snapshot.20131119",b37,true,false))

    //
    // example call set for documentation guide tutorial
    //
    addResource(new Resource("/humgen/gsa-hpprojects/NA12878Collection/exampleCalls/NA12878.HiSeq.WGS.bwa.cleaned.raw.b37.subset.vcf",
      "NA12878.HiSeq.WGS.bwa.cleaned.raw.subset", b37, true, true))

    //
    // Test BAM file, only for the b37 reference
    //
    addResource(new Resource("/humgen/gsa-hpprojects/NA12878Collection/bams/CEUTrio.HiSeq.WGS.b37.NA12878.bam",
      "IGNORE", b37, false, false))

    //
    // Exome targets file, only for the b37 reference
    //
    addResource(new Resource("/seq/references/HybSelOligos/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list",
      "Broad.human.exome", b37, true, false, false))

    //
    // refGene files specific to each reference
    //
    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/refGene_b37.sorted.txt",
      "refGene", b37, true, false))

    addResource(new Resource("public/chainFiles/hg18tob37.chain", "", hg18, false, false))
    addResource(new Resource("public/chainFiles/b36tob37.chain", "", b36, false, false))

    // todo -- chain files?
    // todo 1000G SNP and indel call sets?

    //
    // exampleFASTA file
    //
    addResource(new Resource(exampleFASTA.file, "exampleFASTA", exampleFASTA, false))
    addResource(new Resource("public/testdata/exampleBAM.bam", "exampleBAM", exampleFASTA, false, false, false))
  }

  def createBundleDirectories(dir: File) = {
    if ( ! dir.exists ) dir.mkdirs()

    for ( ref <- refs ) {
      val refDir = new File(dir + "/" + ref.name)
      if ( ! refDir.exists ) refDir.mkdirs()
    }
  }

  def createCurrentLink(bundleDir: File) = {

    val currentLink = new File(BUNDLE_ROOT + "/current")

    if ( currentLink.exists ) add(new deleteLink(currentLink))

    add(new linkFile(bundleDir, currentLink))
  }

  def script = {
    if ( TEST )
      initializeTestDataFiles();
    else
      initializeStandardDataFiles();

    if ( ! DO_DOWNLOAD ) {
      // create destination directory structure
      createBundleDirectories(BUNDLE_DIR)

      for ( resource: Resource <- RESOURCES ) {
        if ( isFASTA(resource.file) ) {
          copyBundleFasta(resource, resource.ref)
        } else if ( isBAM(resource.file) ) {
          val f = copyBundleFile(resource, resource.ref)
          add(new IndexBAM(f))
	  if ( resource.makeCallsIfBam ) {
            @Output val outvcf: File = swapExt(f.getParent, f, ".bam", ".vcf")
            add(new UG(resource.file, resource.ref.file, outvcf))
	  }
        } else if ( isVCF(resource.file) ) {
          for ( destRef <- refs ) {
            val out = destFile(BUNDLE_DIR, destRef, resource.destname(destRef))
            var continue = true

            // copy or lift over the original vcf
            if ( resource.ref == destRef ) {
              add(new cpFile(resource.file, out))
            } else {
              val clf = liftover(resource.file, resource.ref, out, destRef)
              if ( clf != null ) {
                add(clf)
              } else {
                continue = false
              }
            }

            if ( continue ) {
              add(new IndexVCF(out, destRef.file))

              if ( resource.makeSites ) {
                val sites: Resource = new Resource(swapExt(out.getParent, out, ".vcf", "." + SITES_EXT + ".vcf"), "", destRef, false)
                add(new JustSites(out, sites.file))
                add(new IndexVCF(sites.file, destRef.file))
              }

              if ( resource.name.contains("dbsnp") ) {
                val dbsnp129: Resource = new Resource(swapExt(out.getParent, out, ".vcf", ".excluding_sites_after_129.vcf"), "", destRef, false)
                add(new MakeDBSNP129(out, destRef.file, dbsnp129.file))
                add(new IndexVCF(dbsnp129.file, destRef.file))
              }
            }
          }
        } else if ( isIntervalList(resource.file) ) {
            val out = destFile(BUNDLE_DIR, resource.ref, resource.destname(resource.ref))
            add(new cpFile(resource.file, out))
        } else {
          //throw new ReviewedStingException("Unknown file type: " + resource)
        }
      }

    } else {
      createCurrentLink(BUNDLE_DIR)
      createBundleDirectories(DOWNLOAD_DIR)
      createDownloadsFromBundle(BUNDLE_DIR, DOWNLOAD_DIR)
    }
  }


  def createDownloadsFromBundle(in: File, out: File) {
    Console.printf("Visiting %s%n", in)
    if (! in.getName.startsWith(".")) {
      if ( in.isDirectory ) {
        out.mkdirs

        for ( child: File <- in.listFiles ) {
          createDownloadsFromBundle(child, out + "/" + child.getName)
        }
      } else {
        if ( isBAM(in) ) {
          add(new cpFile(in, out))
          add(new md5sum(out))
        } else if ( !isOUT(in) ) {
          add(new GzipFile(in, out + ".gz"))
          add(new md5sum(out + ".gz"))
        }

      }
    }
  }

  def copyBundleFasta(res: Resource, ref: Reference) {
    val out = destFile(BUNDLE_DIR, ref, res.destname(ref))
    add(new cpFile(res.file, out))

    val oldRefDict = swapExt(res.file.getParent, res.file, ".fasta", ".dict")
    val newRefDict = swapExt(out.getParent, out, ".fasta", ".dict")

    val oldRefFai = swapExt(res.file.getParent, res.file, ".fasta", ".fasta.fai")
    val newRefFai = swapExt(out.getParent, out, ".fasta", ".fasta.fai")

    add(new cpFile(oldRefDict, newRefDict))
    add(new cpFile(oldRefFai, newRefFai))
  }

  def copyBundleFile(res: Resource, ref: Reference): File = {
    val out = destFile(BUNDLE_DIR, ref, res.destname(ref))
    add(new cpFile(res.file, out))
    return out
  }

  def destFile(dir: File, ref: Reference, f: File): File = {
    return destFile(dir, ref, f.getName)
  }

  def destFile(dir: File, ref: Reference, name: String): File = {
    return new File(dir + "/" + ref.name + "/" + name)
  }

  /**
   * A command line (cut) that removes all genotyping information from a file
   */
  class JustSites(@Input(doc="foo") in: File, @Output(doc="foo") out: File) extends CommandLineFunction {
    def commandLine = "cut -f 1-8 %s > %s".format(in, out)
  }

  class GzipFile(@Input val in: File, @Output val out: File) extends CommandLineFunction {
    def commandLine = "gzip -c %s > %s".format(in.getAbsolutePath, out.getAbsolutePath)
  }

  class cpFile(@Input val in: File, @Output val out: File) extends CommandLineFunction {
    def commandLine = "cp %s %s".format(in.getAbsolutePath, out.getAbsolutePath)
  }

  class deleteLink(@Input val in: File) extends CommandLineFunction {
    def commandLine = "rm %s".format(in.getAbsolutePath)
  }

  class linkFile(@Input val in: File, @Output val out: File) extends CommandLineFunction {
    def commandLine = "ln -s %s %s".format(in.getAbsolutePath, out.getAbsolutePath)
  }

  class md5sum(@Input val in: File) extends CommandLineFunction {
    @Output val o: File = new File(in.getAbsolutePath + ".md5")
    def commandLine = "md5sum %s > %s".format(in.getAbsolutePath, o)
  }

  class IndexBAM(bamIn: File) extends SamtoolsIndexFunction {
    bamFile = bamIn
  }

  class IndexVCF(@Input vcf: File, @Input ref: File) extends CountRODs with UNIVERSAL_GATK_ARGS {
    //@Output val vcfIndex: File = swapExt(vcf.getParent, vcf, ".vcf", ".vcf.idx")
    this.rod :+= vcf
    this.reference_sequence = ref
  }

  class UG(@Input bam: File, @Input ref: File, @Input outVCF: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file = List(bam)
    this.reference_sequence = ref
    this.intervalsString ++= List("20");
    this.out = outVCF
  }

  class MakeDBSNP129(@Input dbsnp: File, @Input ref: File, @Output dbsnp129: File) extends SelectVariants with UNIVERSAL_GATK_ARGS {
    this.variant = dbsnp
    this.select ++= List("dbSNPBuildID <= 129")
    this.reference_sequence = ref
    this.out = dbsnp129
  }

  class LiftOverPerl(@Input val in: File, @Output val out: File, @Input val chain: File, oldRef: Reference, newRef: Reference) extends CommandLineFunction {
    this.memoryLimit = 12
    def commandLine = ("%s -vcf %s -chain %s -out %s " +
      "-gatk ./ -newRef %s -oldRef %s -tmp %s").format(liftOverPerl, in.getAbsolutePath, chain,
      out.getAbsolutePath, newRef.file.replace(".fasta", ""),
      oldRef.file.replace(".fasta", ""), jobTempDir)
  }

  def getExtension(f: File): String = {
    val i = f.getName.lastIndexOf('.');
    if (i > 0 && i < f.getName.length() - 1)
      return f.getName.substring(i+1).toLowerCase();
    else
      return "";
  }
}

