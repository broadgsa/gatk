package org.broadinstitute.sting.queue.util

import net.sf.picard.reference.ReferenceSequenceFileFactory
import java.io.File
import org.broadinstitute.sting.utils.GenomeLocParser
import collection.JavaConversions._
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.queue.pipeline.PipelineArgumentCollection
import org.broadinstitute.sting.utils.yaml.YamlUtils
import org.broadinstitute.sting.datasources.pipeline.{PipelineSample, PipelineProject, Pipeline}
import org.broadinstitute.sting.utils.text.XReadLines

class PipelineUtils {

}

object PipelineUtils{

  def smartSplitContigs(reference: File, intervals: File, sets: Int) : List[List[String]] = {
    var genomeLocParser: GenomeLocParser = new GenomeLocParser(ReferenceSequenceFileFactory.getReferenceSequenceFile(reference))
    val targets = IntervalUtils.parseIntervalArguments(genomeLocParser,List(intervals.getAbsolutePath), false)

    // Build up a map of contigs with sizes.
    var contigSizes = Map.empty[String, Long]
    // todo -- make this look like functional code for Q's sake
    //targets.foreach( loc => { contigSizes += loc -> { contigSizes.get(loc.getContig) match { case Some(size) => size + loc.size case None => loc.size } } })

    for (loc <- targets) {
      val contig = loc.getContig
      val contigSize = loc.size
      contigSizes += contig -> {
        contigSizes.get(contig) match {
          case Some(size) => size + contigSize
          case None => contigSize
        }
      }
    }

    // Keep a list of pairs of sizes with lists of contigs.
    var splitContigs = List.empty[(Long, List[String])]
    for ((contigName, contigSize) <- contigSizes) {
      if (splitContigs.size < sets) {
        // If there are fewer than the requested number of sets, just add this contig.
        splitContigs :+= contigSize -> List(contigName)
      } else {
        // If there is already a number of sets
        // sort the contigs to get the smallest one first.
        splitContigs = splitContigs.sortBy{case (size, contigs) => size}
        // Update the pair with the new contig size and name.
        var smallContigs = splitContigs.head
        smallContigs = (smallContigs._1 + contigSize) -> (smallContigs._2 :+ contigName)
        // Re add the pair to the list.
        splitContigs = smallContigs :: splitContigs.tail
      }
    }

    splitContigs.map{case (size, contigs) => contigs}
  }

  def loadPipelineFromPAC(args: PipelineArgumentCollection) : Pipeline = {
    if ( args.yamlFile != null ) {
      return YamlUtils.load(classOf[Pipeline], args.yamlFile)
    } else {
      return loadPipelineFromSpec(args.projectName,args.projectRef,args.projectIntervals,args.projectDBSNP,args.projectBams)
    }
  }

  def loadPipelineFromSpec(name: String, ref: File, ivals: File, dbsnp: File, pBamList: File) : Pipeline = {
    var newPipeline : Pipeline = new Pipeline
    var pipeProject : PipelineProject = new PipelineProject
    var pipeSamples : List[PipelineSample] = ((new XReadLines(pBamList)).readLines).toList.map( bamSpecToSample )

    pipeProject.setName(name)
    pipeProject.setReferenceFile(ref)
    pipeProject.setIntervalList(ivals)
    pipeProject.setDbsnpFile(dbsnp)

    newPipeline.setProject(pipeProject)
    newPipeline.setSamples(pipeSamples)
    
    return newPipeline
  }

  //todo -- find a better name for this function
  def bamSpecToSample(spec: String) : PipelineSample = {
    var sam : PipelineSample = new PipelineSample
    var spStr : Array[String] = spec.split("\\s")
    sam.setId(spStr(0))
    var tagStr : Array[String] = spStr(1).split(",")
    var tagMap : java.util.HashMap[String,String] = new java.util.HashMap[String,String](tagStr.size)
    tagStr.filter( u => ! u.equals("")).foreach( u => tagMap.put(u.split(":")(0),u.split(":")(1)) )
    sam.setTags(tagMap)
    var bamStr : Array[String] = spStr(2).split(",")
    var bamMap : java.util.HashMap[String,File] = new java.util.HashMap[String,File](bamStr.size)
    bamStr.foreach( u => bamMap.put( u.split(":")(0), new File(u.split(":")(1) ) ) )
    sam.setBamFiles(bamMap)

    return sam

  }
}