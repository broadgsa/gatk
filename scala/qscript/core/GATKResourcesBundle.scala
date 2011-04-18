package core

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException
import org.broadinstitute.sting.queue.extensions.picard.PicardBamFunction
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction

class GATKResourcesBundle extends QScript {
  // todo -- update to released version when things stabilize
  @Argument(doc="gatkJarFile", required=false)
  var gatkJarFile: File = new File("dist/GenomeAnalysisTK.jar")

  @Argument(doc="liftOverPerl", required=false)
  var liftOverPerl: File = new File("./perl/liftOverVCF.pl")

  @Argument(shortName = "bundleDir", doc="Path to root where resource files will be placed", required=false)
  val BUNDLE_DIR = new File("bundle")

  @Argument(shortName = "downloadDir", doc="Path to root where resource files will be placed for users to download", required=false)
  val DOWNLOAD_DIR = new File("downloads")

  @Argument(shortName = "test", doc="", required=false)
  val TEST = false

  @Argument(shortName = "phase2", doc="", required=false)
  val DO_DOWNLOAD = false

  val SITES_EXT: String = "sites"

  class Reference( val name: String, val file: File ) { }

  val b37 = new Reference("b37", new File("/Users/depristo/Desktop/broadLocal/localData/human_g1k_v37.fasta"))
  //val b36 = new Reference("b36", new File("/Users/depristo/Desktop/broadLocal/localData/human_b36_both.fasta"))
  //val hg19 = new Reference("hg19")
  val hg18 = new Reference("hg18", new File("/Users/depristo/Desktop/broadLocal/localData/Homo_sapiens_assembly18.fasta"))
  val exampleFASTA = new Reference("exampleFASTA", new File("testdata/exampleFASTA.fasta"))

  val refs = List(b37, hg18, exampleFASTA)
  //val refs = List(b37, b36, hg19, hg18, exampleFASTA)

  class Resource(val file: File, val name: String, val ref: Reference, val useName: Boolean = true, val makeSites: Boolean = true ) {
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
//      case ("b37", "hg19") =>
//        return new HGRFC2UCSC(in, out)
      case ("b37", "hg18") =>
        return new LiftOverPerl(in, out, new File("chainFiles/b37Tohg18.chain"), inRef, outRef)
//      case ("b37", "b36") =>
//        return new LiftOverPerl(in, out, new File("chainFiles/b37Tob36.chain"), inRef, outRef)
      case _ => return null
    }
  }

  def isVCF(file: File) = file.getName.endsWith(".vcf")
  def isBAM(file: File) = file.getName.endsWith(".bam")
  def isFASTA(file: File) = file.getName.endsWith(".fasta")

  var RESOURCES: List[Resource] = Nil
  def addResource(comp: Resource) { RESOURCES = comp :: RESOURCES }

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.jarFile = gatkJarFile;
    this.memoryLimit = 2
  }

  def initializeTestDataFiles() = {
    //
    // Standard evaluation files for indels
    //
    val DATAROOT = "/Users/depristo/Desktop/broadLocal/localData/"
    addResource(new Resource(DATAROOT + "human_g1k_v37.fasta", "human_g1k_v37.fasta", b37, false))
    addResource(new Resource(DATAROOT + "1000G.snp.validation.b37.vcf", "1000G.snp.validation", b37))
    //addResource(new Resource(DATAROOT + "dbsnp_132_b37.vcf", "dbsnp_132", b37, true, false))

    addResource(new Resource("testdata/exampleFASTA.fasta", "exampleFASTA", exampleFASTA, false))
    addResource(new Resource("testdata/exampleBAM.bam", "exampleBAM", exampleFASTA, false))
  }

  def initializeStandardDataFiles() = {
    // references
    addResource(new Resource("/humgen/1kg/reference/human_g1k_v37.fasta", "human_g1k_v37.fasta", b37, false))
    addResource(new Resource("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta", "Homo_sapiens_assembly18.fasta", hg18, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/dbSNP/dbsnp_132_b37.leftAligned.vcf",
      "dbsnp_132", b37, true, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/Omni2.5_chip/Omni25_genotypes_1525_samples.b37.vcf",
      "1000G_omni2.5", b37, true, true))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.3/genotypes_r27_nr.b37_fwd.vcf",
      "hapmap_3.3", b37, true, true))

    addResource(new Resource("AFR+EUR+ASN+1KG.dindel_august_release_merged_pilot1.20110126.sites.vcf",
      "1000G_indels_for_realignment", b37, true, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/refGene_b37.sorted.txt",
      "refGene", b37, true, false))

    addResource(new Resource("/humgen/gsa-hpprojects/GATK/data/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam",
      "IGNORE", b37, false, false))

    // todo -- chain files?
    // todo 1000G SNP and indel call sets?

    addResource(new Resource("testdata/exampleFASTA.fasta", "exampleFASTA", exampleFASTA, false))
    addResource(new Resource("testdata/exampleBAM.bam", "exampleBAM", exampleFASTA, false))
  }

  def createBundleDirectories(dir: File) = {
    if ( ! dir.exists ) dir.mkdirs()

    for ( ref <- refs ) {
      val refDir = new File(dir + "/" + ref.name)
      if ( ! refDir.exists ) refDir.mkdirs()
    }
  }

  def script = {
    if ( ! DO_DOWNLOAD ) {
      createBundleDirectories(BUNDLE_DIR)
      if ( TEST )
        initializeTestDataFiles();
      else
        initializeStandardDataFiles();

      for ( resource <- RESOURCES ) {
        if ( isFASTA(resource.file) ) {
          val f = copyBundleFile(resource, resource.ref)
          add(new createDictandFAI(f))
        } else if ( isBAM(resource.file) ) {
          val f = copyBundleFile(resource, resource.ref)
          add(new IndexBAM(f))
          @Output val outvcf: File = swapExt(f.getParent, f, ".bam", ".vcf")
          add(new UG(resource.file, resource.ref.file, outvcf))
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
            }
          }
        } else {
          //throw new ReviewedStingException("Unknown file type: " + resource)
        }
      }
    } else {
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
        } else {
          add(new GzipFile(in, out + ".gz"))
          add(new md5sum(out + ".gz"))
        }

      }
    }
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

  class md5sum(@Input val in: File) extends CommandLineFunction {
    @Output val o: File = new File(in.getAbsolutePath + ".md5")
    def commandLine = "md5 %s > %s".format(in.getAbsolutePath, o)
  }

  class IndexBAM(bamIn: File) extends SamtoolsIndexFunction {
    bamFile = bamIn
  }

  class IndexVCF(@Input vcf: File, @Input ref: File) extends CountRod with UNIVERSAL_GATK_ARGS {
    this.rodBind :+= RodBind(vcf.getName, "VCF", vcf)
    this.reference_sequence = ref
  }

  class UG(@Input bam: File, @Input ref: File, @Input outVCF: File) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    this.input_file = List(bam)
    this.reference_sequence = ref
    this.out = outVCF
  }

  class LiftOverPerl(@Input val in: File, @Output val out: File, @Input val chain: File, oldRef: Reference, newRef: Reference) extends CommandLineFunction {
    def commandLine = ("%s -vcf %s -chain %s -out %s " +
      "-gatk ./ -newRef %s -oldRef %s ").format(liftOverPerl, in.getAbsolutePath, chain,
      out.getAbsolutePath, newRef.file.replace(".fasta", ""), oldRef.file.replace(".fasta", ""))
  }

  class HGRFC2UCSC(@Input val in: File, @Output val out: File) extends CommandLineFunction {
    def commandLine = "python ../python/vcf_b36_to_hg18.py %s %s".format(in.getAbsolutePath, out.getAbsolutePath)
  }

  def getExtension(f: File): String = {
    val i = f.getName.lastIndexOf('.');
    if (i > 0 && i < f.getName.length() - 1)
      return f.getName.substring(i+1).toLowerCase();
    else
      return "";
  }

  class createDictandFAI (@Input ref: File) extends FastaStats with UNIVERSAL_GATK_ARGS {
    @Output val outDict: File = swapExt(ref.getParent, ref, ".fasta", ".dict")
    @Output val outFai: File = swapExt(ref.getParent, ref, ".fasta", ".fasta.fai")
    @Output val outStats: File = swapExt(ref.getParent, ref, ".fasta", ".stats")
    this.reference_sequence = ref
    this.out = outStats
  }
}

