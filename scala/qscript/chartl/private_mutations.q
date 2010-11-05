import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils
import org.broadinstitute.sting.queue.extensions.gatk.{UnifiedGenotyper, RodBind, CombineVariants, SelectVariants}
import org.broadinstitute.sting.queue.QScript
import tools.nsc.io.File

class private_mutations extends QScript {
  pm_script => // alias for the script arguments

  val eomi_merged_calls : File = new File("/humgen/gsa-hphome1/chartl/projects/private_mutations/resources/esp_merged.vcf")
  val g1k_exome_calls : File = new File("/humgen/gsa-hphome1/chartl/projects/private_mutations/resources/g1k_exomes.vcf")
  val ceu_pilot_calls : File = new File("/humgen/gsa-hphome1/chartl/projects/private_mutations/resources/CEU.low_coverage.genotypes.vcf")

  val SINGLE_SAMPLE_EXPRESSION = "\"PASS.*AC=1;|PASS.*AC=2;.*1/1|PASS.*AC=2.*1\\|1\""
  val NOVEL_VARIANT_EXPRESSION = "-v \"rs[0-9]\""
  val KNOWN_VARIANT_EXPRESSION = "\"rs[0-9].*PASS\""

  val BASE_DIR = "/humgen/gsa-hphome1/chartl/projects/private_mutations/queue/"
  val LIST_DIR = BASE_DIR + "lists/"
  val VCF_DIR = BASE_DIR + "vcfs/"
  val EVAL_DIR = BASE_DIR + "evals/"
  val SCRATCH_DIR = BASE_DIR + "intermediate/"
  
  /*
   * Steps of analysis.
   *
   * todo -- 1) Calculate number of EOMI calls top-down, cumulative include/exclude that are:
   *  - present in one sample
   *  - not in dbsnp 129
   *  - not in 1000G pilot CEU
   *  - not in 1000G production EUR
   *  - not in 1000G production exome EUR
   *
   * todo -- 2) Distributional aspects of remaining mutations
   *  - Depth of Coverage
   *  - Functional Class
   *
   * todo -- 3) 1000G Exome v Lowpass -- Venn of single-sample variants
   *  - Venn numbers
   *  - Distribution of depth for exome variants missed in low-pass
   *
   * todo -- 4) Estimate relationship between depth and sensitivity to private variation
   *
   * todo -- 5) Calculate % of exome empowered to detect private variation
   *
   * todo -- 6) Identify loci across 1000G Ex and EOMI that are low-powered
   *   - venn of low-power targets
   *   - mapping the intersection
   *   - genes affected by low power
   */


  /*
   * Scripting commands go here
   */
  def script = {
    // setup analysis resources (best done in loop)
    val g1k_lowpass_chr : List[File] = (1 to 22).toList.map( i => new File("/humgen/1kg/processing/allPopulations_wholeGenome_august_release/calls/chr%d/ALL.august.chr%d.recal.vrcut.EUR.vcf".format(i,i)))
    val eomi_annot_set : List[File] = (1 to 7).toList.map( i => new File("/humgen/gsa-hphome1/chartl/projects/private_mutations/resources/esp_set%d.vcf".format(i)))
    val g1k_lowpass_hg19_merged : File = new File(VCF_DIR+"g1k_lowpass_merged.hg19.sites.vcf")
    val eomi_merged_hg19 : File = eomi_merged_calls
    val ceu_hg19 : File = new File(VCF_DIR+"G1K_Pilot.CEU.hg19.vcf")

    combineAndMergeSites(g1k_lowpass_chr,g1k_lowpass_hg19_merged)
    liftoverVCF(ceu_pilot_calls,ceu_hg19)



    // get EOMI ss-only list
    var eomi_single_sample : Vcf2List = new Vcf2List
    eomi_single_sample.in_vcf = eomi_merged_hg19
    eomi_single_sample.filter = SINGLE_SAMPLE_EXPRESSION
    eomi_single_sample.out_list = new File(LIST_DIR+"EOMI_single_sample_variants.sites.list")

    add(eomi_single_sample)

    var eomi_dbsnp : Vcf2List = new Vcf2List
    eomi_dbsnp.in_vcf = eomi_merged_hg19
    eomi_dbsnp.filter = KNOWN_VARIANT_EXPRESSION
    eomi_dbsnp.out_list = new File(LIST_DIR+"EOMI_novel_variants.sites.list")

    add(eomi_dbsnp)

    var ceu_sites : Vcf2List = new Vcf2List
    ceu_sites.in_vcf = ceu_hg19
    ceu_sites.out_list = new File(SCRATCH_DIR+"G1KP_CEU_variants.sites.list")

    add(ceu_sites)

    var g1k_lowpass_sites : Vcf2List = new Vcf2List
    g1k_lowpass_sites.in_vcf = g1k_lowpass_hg19_merged
    g1k_lowpass_sites.out_list = new File(SCRATCH_DIR+"G1K_lowpass_EUR.sites.list")

    add(g1k_lowpass_sites)

    var g1k_exome_sites : Vcf2List = new Vcf2List
    g1k_exome_sites.in_vcf = g1k_exome_calls
    g1k_exome_sites.out_list = new File(SCRATCH_DIR+"G1K_exomes.sites.list")

    add(g1k_exome_sites)

    var remove_dbsnp : RemoveIntersect = new RemoveIntersect
    remove_dbsnp.fileToUnique = eomi_single_sample.out_list
    remove_dbsnp.comparisonFile = eomi_dbsnp.out_list
    remove_dbsnp.outputFile = new File(LIST_DIR+"EOMI_single_sample.nodbsnp.list")

    add(remove_dbsnp)

    var remove_ceu : RemoveIntersect = new RemoveIntersect
    remove_ceu.fileToUnique = remove_dbsnp.outputFile
    remove_ceu.comparisonFile = ceu_sites.out_list
    remove_ceu.outputFile = new File(LIST_DIR+"EOMI_single_sample.nodbsnp.noceu.list")

    add(remove_ceu)

    var remove_lowpass : RemoveIntersect = new RemoveIntersect
    remove_lowpass.fileToUnique = remove_ceu.outputFile
    remove_lowpass.comparisonFile = g1k_lowpass_sites.out_list
    remove_lowpass.outputFile = new File(LIST_DIR+"EOMI_single_sample.nodbsnp.noceu.noeur.list")

    add(remove_lowpass)

    var remove_exome : RemoveIntersect = new RemoveIntersect
    remove_exome.fileToUnique = remove_lowpass.outputFile
    remove_exome.comparisonFile = g1k_exome_sites.out_list
    remove_exome.outputFile = new File(LIST_DIR+"EOMI_single_sample.nodbsnp.noceu.noeur.noexome.list")

    add(remove_exome)

    var subset : SelectVariants = new SelectVariants
    subset.jarFile = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")
    subset.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
    subset.variantVCF = eomi_merged_hg19
    subset.out = new File(VCF_DIR+"EOMI_merged.hg19.private.vcf")
    subset.intervals :+= remove_exome.outputFile

    add(subset)

    var varDepths : GetVariantDepths = new GetVariantDepths
    varDepths.inVCF = subset.out
    varDepths.outFile = new File("/humgen/gsa-hphome1/chartl/projects/private_mutations/results/EOMI.private.variant.depths.txt")

    add(varDepths)

    var lowpass_single : Vcf2List = new Vcf2List
    lowpass_single.in_vcf = g1k_lowpass_hg19_merged
    lowpass_single.filter = SINGLE_SAMPLE_EXPRESSION
    lowpass_single.out_list = new File(LIST_DIR+"g1k_lowpass_single_sample.variants.list")

    add(lowpass_single)

    var exome_single : Vcf2List = new Vcf2List
    exome_single.in_vcf = g1k_exome_calls
    exome_single.filter = SINGLE_SAMPLE_EXPRESSION
    exome_single.out_list = new File(LIST_DIR+"g1k_exome_single_sample.variants.list")

    add(exome_single)

    var exome_extract : ExtractSites = new ExtractSites(g1k_exome_calls,new File(VCF_DIR+"g1k_exome.sites.vcf"))
    add(exome_extract)

    var getVennExS : CombineVariants = new CombineVariants
    getVennExS.rodBind :+= new RodBind("lowpass","VCF",g1k_lowpass_hg19_merged);
    getVennExS.rodBind :+= new RodBind("exome","VCF",exome_extract.outputVCF);
    getVennExS.priority = "exome,lowpass"
    getVennExS.intervals :+= exome_single.out_list
    getVennExS.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
    getVennExS.jarFile = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")
    getVennExS.genotypeMergeOptions = Some(VariantContextUtils.GenotypeMergeType.UNIQUIFY)
    getVennExS.variantMergeOptions = Some(VariantContextUtils.VariantMergeType.UNION)
    getVennExS.out = new File(VCF_DIR + "g1k_exome_plus_lowpass.singlesample.merged.exome.sites.vcf")

    //add(getVennExS)

    var getVennLPS : CombineVariants = new CombineVariants
    getVennLPS.rodBind :+= new RodBind("lowpass","VCF",g1k_lowpass_hg19_merged);
    getVennLPS.rodBind :+= new RodBind("exome","VCF",exome_extract.outputVCF);
    getVennLPS.priority = "exome,lowpass"
    getVennLPS.intervals :+= lowpass_single.out_list
    getVennLPS.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
    getVennLPS.jarFile = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")
    getVennLPS.genotypeMergeOptions = Some(VariantContextUtils.GenotypeMergeType.UNIQUIFY)
    getVennLPS.variantMergeOptions = Some(VariantContextUtils.VariantMergeType.UNION)
    getVennLPS.out = new File(VCF_DIR + "g1k_exome_plus_lowpass.singlesample.merged.lowpass.sites.vcf")

    add(getVennLPS)

    var getG1KOverlap : GetOverlapSamples = new GetOverlapSamples
    getG1KOverlap.vcf1 = g1k_exome_calls
    getG1KOverlap.vcf2 = g1k_lowpass_chr(1) // EUR only
    getG1KOverlap.samples = new File(LIST_DIR+"g1k_EUR_exome_lowpass_overlap.txt")

    add(getG1KOverlap)

    //callOverlaps(getG1KOverlap.samples,exome_single.out_list, g1k_exome_calls)

  }

  /*
   * defined modules go here
   */

  class Vcf2List extends CommandLineFunction {
    @Input(doc="The vcf file to convert to list",required=true)
    var in_vcf: File = _
    @Argument(doc="An egrep-based filter to apply during conversion",required=false)
    var filter: String = "PASS"
    @Output(doc="The list file to write to", required=true)
    var out_list: File = _

    def commandLine = {
      "egrep %s %s | awk '{print $1\":\"$2}' > %s".format(filter,in_vcf.getAbsolutePath,out_list.getAbsolutePath)
    }
  }

  class GetOverlapSamples extends CommandLineFunction {
    @Input(doc="vcf1")
    var vcf1: File = _
    @Input(doc="vcf2")
    var vcf2: File = _
    @Output(doc="sample file")
    var samples: File = _

    def commandLine = {
      "head -n 500 %s %s | grep \\\\#CHR | cut -f10- | tr '\\t' '\\n' | sort | uniq -c | awk '{if ($1 == 2) print $2}' > %s".format(
        vcf1.getAbsolutePath,
        vcf2.getAbsolutePath,
        samples.getAbsolutePath
        )
    }
  }

  class RemoveIntersect extends CommandLineFunction {
    @Input(doc="The vcf file whose unique entries should be retained",required=true)
    var fileToUnique: File = _
    @Input(doc="The vcf file whose entries overlapping the unique file you would like to remove",required=true)
    var comparisonFile: File = _
    @Output(doc="The file containing only the fileToUnique-unique entries", required=true)
    var outputFile: File = _

    def commandLine = {
      var tmpFile: File = java.io.File.createTempFile("removeintersect","tmp")

      "cat %s %s | sort | uniq -c | awk '{if ($1==2) print $2}' > %s ; cat %s %s | sort | uniq -c | awk '{if ($1 == 1) print $2}' > %s".format(
      fileToUnique.getAbsolutePath,comparisonFile.getAbsolutePath, tmpFile.getAbsolutePath,
      tmpFile.getAbsolutePath, fileToUnique.getAbsolutePath, outputFile.getAbsolutePath
        )
    }
  }

  class ExtractSites(inVCF: File, oVCF: File) extends CommandLineFunction {
    @Input(doc="The vcf file to generate a sites only file from",required=true)
    var inputVCF: File = inVCF
    @Output(doc="The sites-only vcf to write to")
    var outputVCF: File = oVCF

    def commandLine = {
      "cut -f1-9 %s > %s".format(inputVCF.getAbsolutePath,outputVCF.getAbsolutePath)
    }
  }

  class ConcatVCF(in: List[File], out: File) extends CommandLineFunction {
    @Input(doc="The files to concatenate",required=true)
    var inputVCFs : List[File] = in
    @Output(doc="The file to write to", required=true)
    var outputVCF: File = out

    def commandLine = {
      var header : File = java.io.File.createTempFile("concatVCF","header.tmp")
      var body_unsorted : File = java.io.File.createTempFile("concatVCF","body.unsorted.tmp")
      var body : File = java.io.File.createTempFile("concatVCF","body.tmp")

      "grep \\\\# %s > %s ; grep -v \\\\# %s | sed 's/.vcf:/\\t/g' | cut -f2- > %s ; perl %s -tmp . %s %s > %s ; cat %s %s > %s".format(
        inputVCFs(0).getAbsolutePath,
        header.getAbsolutePath,
        inputVCFs.foldLeft[String]("")((x1: String, x2: File ) => x1 + " " + x2.getAbsolutePath ),
        body_unsorted.getAbsolutePath,
        "/humgen/gsa-scr1/chartl/sting/perl/sortByRef.pl",
        body_unsorted.getAbsolutePath,
        "/humgen/1kg/reference/human_g1k_v37.fasta.fai",
        body.getAbsolutePath,
        header.getAbsolutePath,
        body.getAbsolutePath,
        outputVCF.getAbsolutePath
        )
    }
  }

  class ExtractSample(vcfIn: File, sample: String, vcfOut: File) extends CommandLineFunction {
    @Input(doc="File from which to extract sample")
    var inVCF: File = vcfIn
    @Argument(doc="The sample to extract")
    var inSample: String = sample
    @Output(doc="The VCF file to write to")
    var outVCF: File = vcfOut

    def commandLine = {
      "head -n 500 %s | grep \\\\#CHR | tr '\\t' '\\n' | grep -n %s | tr ':' '\\t' | cut -f1 | xargs -i cut -f1-9,\\{\\} %s | egrep \"\\\\#|PASS.*0/1|PASS.*1/0|PASS.*1/1\" > %s".format(
        inVCF.getAbsolutePath, inSample, inVCF.getAbsolutePath, outVCF.getAbsolutePath
        )
    }
  }

  class GetVariantDepths extends CommandLineFunction {
    @Input(doc="The VCF to get the per-sample depth at variant calls")
    var inVCF : File = _
    @Output(doc="The depth file to write to")
    var outFile : File = _

    def commandLine = {
      "grep PASS %s | cut -f10- | tr '\\t' '\\n' | egrep \"1/0|0/1|1/1|1\\|0|0\\|1|1\\|1\" | tr ':' '\\t' | awk '{print $3}' > %s".format(
        inVCF.getAbsolutePath, outFile.getAbsolutePath
        )
    }
  }

  def combineAndMergeSites(chrInputs: List[File], output: File) = {
    var sites_only: List[ExtractSites] = chrInputs.map( f => new ExtractSites(f,new File(VCF_DIR+swapExt(f,".vcf",".sites.vcf").getName)))
    for ( s <- sites_only ) {
      add(s)
    }
    var trivialMerge : ConcatVCF = new ConcatVCF(sites_only.map( s => s.outputVCF ), output)
    add(trivialMerge)
  }

  def liftoverVCF(b36vcf: File, b37vcf: File) = {
    
    class Liftover(in: File, out: File) extends CommandLineFunction {
      @Input(doc="foo")
      val input: File = in
      @Output(doc="foo")
      val output : File = out

      val sting : String = "/humgen/gsa-scr1/chartl/sting/"
      val chain : String = "/humgen/gsa-hpprojects/GATK/data/Liftover_Chain_Files/b36ToHg19.broad.over.chain"

      def commandLine = {
        "%s/perl/liftOverVCF.pl -vcf %s -chain %s -out %s -gatk %s -newRef %s -oldRef %s -tmp .".format(
        sting, input, chain, output, sting, "/seq/references/Homo_sapiens_assembly19/v0/Homo_sapiens_assembly19",
          "/humgen/1kg/reference/human_b36_both"
          )
      }

      class RMDupes(in: File, out:File) extends CommandLineFunction {
        @Input(doc="foo")
        var input : File = in
        @Output(doc="foo")
        var output : File = out

        def commandLine = "cat %s | python /humgen/gsa-scr1/chartl/projects/omni/scripts/filterDupes.py > %s".format(
        input.getAbsolutePath,output.getAbsolutePath
          )
      }

      var lift : Liftover = new Liftover(b36vcf,swapExt(b37vcf,".vcf",".undeduped.vcf"))

      add(lift)

      var rmDupes : RMDupes = new RMDupes(lift.output,b37vcf)

      add(rmDupes)
    }
  }

  def callOverlaps(samples: File, intervals: File, wexCalls: File)  = {
    class GetBamList(sample: String, oList: File) extends CommandLineFunction {
      @Argument(doc="foo")
      var inSample = sample
      @Output(doc="foo")
      var outList = oList
      @Input(doc = "foo")
      var waitForMe: File = _

      def commandLine = {
        "grep %s /humgen/1kg/processing/allPopulations_wholeGenome_august_release/bamLists/*.list | tr ':' '\\t' | awk '{print $2}' | sort | uniq > %s".format(
        inSample,
        outList.getAbsolutePath
          )
      }
    }
    for ( s <- scala.io.Source.fromFile(samples).getLines ) {
      var extract : ExtractSample = new ExtractSample(wexCalls,s,new File(VCF_DIR+"g1k_exome.subset.%s.vcf".format(s)))
      add(extract)

      var sites : Vcf2List = new Vcf2List
      sites.in_vcf = extract.outVCF
      sites.out_list = new File(SCRATCH_DIR+"g1k_exome.subset.%s.variants.list".format(s))
      add(sites)

      var bamList : GetBamList = new GetBamList(s,new File(SCRATCH_DIR+"%s.bams.list".format(s)))
      bamList.waitForMe = sites.out_list

      add(bamList)

      var genotype : UnifiedGenotyper = new UnifiedGenotyper
      genotype.jarFile = new File("/humgen/gsa-scr1/chartl/sting/dist/GenomeAnalysisTK.jar")
      genotype.reference_sequence = new File("/humgen/1kg/reference/human_g1k_v37.fasta")
      genotype.intervals :+= sites.out_list
      genotype.out = new File(VCF_DIR+"g1k_lowpass.%s.exome_sites.vcf".format(s))
      genotype.input_file :+= bamList.outList
      genotype.memoryLimit = Some(3)
      genotype.output_all_callable_bases = true

      add(genotype)

    }
  }

}