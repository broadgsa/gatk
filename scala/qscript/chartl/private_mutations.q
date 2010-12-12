import collection.JavaConversions._
import java.io.FileNotFoundException
import org.broadinstitute.sting.datasources.pipeline._
import org.broadinstitute.sting.queue.extensions.gatk.{VariantFiltration, UnifiedGenotyper}
import org.broadinstitute.sting.queue.library.clf.vcf._
import org.broadinstitute.sting.queue.pipeline._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.yaml.YamlUtils

class private_mutations extends QScript {
  @Argument(shortName="yaml",fullName="eomiYaml",doc="Project YAML file",required=true) var eomiYaml: File = _
  @Argument(shortName="sting",fullName="stingDir",doc="path to the Sting directory",required=true) var sting: String = _
  @Argument(shortName="out",fullName="finalVCF",doc="the merged vcf to write to", required=true) var finalMergedVCF : File = _

  var gatkjar : File = _
  def script = {
    gatkjar = new File(sting+"/dist/GenomeAnalysisTK.jar")
    var input_pipeline : Pipeline = YamlUtils.load(classOf[Pipeline],eomiYaml)
    var eomi_pipeline : Pipeline = new Pipeline
    // use only QC-positive samples
    eomi_pipeline.setProject( input_pipeline.getProject )
    eomi_pipeline.setSamples( input_pipeline.getSamples.filter( p => p.getTags.get("QCStatus").equals("PASS")) )
    var vcLib : VariantCalling = new VariantCalling(eomi_pipeline,gatkjar)
    var pmLib : ProjectManagement = new ProjectManagement(sting)

    /*eomi_pipeline.getSamples.foreach( p =>
      if ( ! p.getBamFiles.get("recalibrated").exists) throw new FileNotFoundException(
        p.getBamFiles.get("recalibrated").getAbsolutePath+" does not exist" ))*/

    var batches : List[List[PipelineSample]] = eomi_pipeline.getSamples.toList.grouped(100).toList
    var genotypers : List[UnifiedGenotyper] = batches.map( pl => pl.map( p => p.getBamFiles.get("recalibrated") ) ).zipWithIndex.map(
    b => vcLib.StandardUnifiedGenotyper(b._1,new File(eomi_pipeline.getProject.getName+"_batch%d.raw.vcf".format(1+b._2))))
    addAll(genotypers)

    var handFilters : List[VariantFiltration] = genotypers.map( g => vcLib.StandardHandfilter(g.out,swapExt(g.out,".raw.vcf",".handfiltered.vcf")))

    addAll(handFilters)

    addAll(pmLib.MergeBatches(handFilters.map( _.out), batches.flatten.map( p => p.getBamFiles.get("recalibrated")),
      finalMergedVCF,eomi_pipeline.getProject.getReferenceFile,20))

    var afr_sams : List[PipelineSample] = eomi_pipeline.getSamples.toList.filter( p => p.getTags.get("Population").equals("AFRAMR"))
    var eur_sams : List[PipelineSample] = eomi_pipeline.getSamples.toList.filter( p => p.getTags.get("Population").equals("EURAMR") ||
      p.getTags.get("Population").equals("UNK"))

    var variant_loci : VCFExtractIntervals = new VCFExtractIntervals(finalMergedVCF,swapExt(finalMergedVCF,".vcf",".intervals.list"),false)

    add(variant_loci)

    var extract_afr : VCFExtractSamples = new VCFExtractSamples(finalMergedVCF,swapExt(finalMergedVCF,".vcf",".afr.vcf"),afr_sams.map(_.getId))
    var extract_eur : VCFExtractSamples = new VCFExtractSamples(finalMergedVCF,swapExt(finalMergedVCF,".vcf",".eur+unk.vcf"),eur_sams.map(_.getId))

    add(extract_afr)
    add(extract_eur)
  }
  
}