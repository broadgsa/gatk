import org.broadinstitute.sting.commandline.ArgumentCollection
import org.broadinstitute.sting.datasources.pipeline.Pipeline
import org.broadinstitute.sting.queue.extensions.gatk.CommandLineGATK
import org.broadinstitute.sting.queue.pipeline._
import org.broadinstitute.sting.queue.util.PipelineUtils
import org.broadinstitute.sting.queue.{QException, QScript}
import collection.JavaConversions._
import org.broadinstitute.sting.utils.yaml.YamlUtils

class fullCallingPipelineV2 extends QScript {
  fcp =>

  @ArgumentCollection var pipelineArgs = new PipelineArgumentCollection
  
  private var pipeline: Pipeline = _

  def script = {
    pipelineArgs.verifyArguments
    pipeline = PipelineUtils.loadPipelineFromPAC(pipelineArgs)
    var callingLib: VariantCalling = new VariantCalling(pipeline,fcp.pipelineArgs.gatkJar)
    var cleaningLib: BamProcessing = new BamProcessing(pipeline,fcp.pipelineArgs.gatkJar,fcp.pipelineArgs.picardFixMatesJar)

    val projectBase: String = fcp.pipeline.getProject.getName
    val cleanedBase: String = projectBase + ".cleaned"
    val uncleanedBase: String = projectBase + ".uncleaned"

    // there are commands that use all the bam files
    val recalibratedSamples = fcp.pipeline.getSamples.filter( u => ( u.getBamFiles.contains("recalibrated") || u.getBamFiles.contains("cleaned") ) )

    var bamsToClean: List[(File,File)] = Nil
    var recalBams: List[File] = Nil
    var cleanedBams: List[File] = Nil

    for ( sample <- recalibratedSamples ) {
      val bam = sample.getBamFiles.get("recalibrated")
      recalBams :+= bam
      if (!sample.getBamFiles.contains("cleaned")) {
        sample.getBamFiles.put("cleaned", swapExt(bam,"bam","cleaned.bam"))
        bamsToClean :+= (bam,sample.getBamFiles.get("cleaned"))
      }

      cleanedBams :+= sample.getBamFiles.get("cleaned")
    }

    if ( !fcp.pipelineArgs.skip_cleaning ) {
      addAll(cleaningLib.StandardIndelRealign(bamsToClean,fcp.pipelineArgs.cleaningJobs))
    }

    if (!fcp.pipelineArgs.skip_cleaning ) {
      endToEnd(cleanedBase, cleanedBams, callingLib)
    } else {
      endToEnd(uncleanedBase, recalBams, callingLib)
    }
  }


  def endToEnd(base: String, bamFiles: List[File], lib: VariantCalling) = {
    var recal_vcf = new File(base+"_snps.recal.annotated.tranched.vcf")
    var handfilt_vcf = new File(base+"_snps.handfiltered.annotated.vcf")
    var indel_vcf = new File(base+"_indel_calls.vcf")

    addAll(lib.StandardCallingPipeline(bamFiles,indel_vcf,recal_vcf,handfilt_vcf,fcp.pipelineArgs.target_titv,fcp.pipelineArgs.refseqTable))
  }
}
