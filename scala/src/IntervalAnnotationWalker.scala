package org.broadinstitute.sting.scala

import java.io.PrintStream
import org.broadinstitute.sting.gatk.contexts.{ReferenceContext, AlignmentContext}
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker
import collection.JavaConversions._
import collection.immutable.List
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqFeature
import org.broadinstitute.sting.gatk.walkers.{TreeReducible, RefWalker}
import org.broadinstitute.sting.commandline.{Output, Argument}
import org.broadinstitute.sting.utils.{BaseUtils, GenomeLoc}
import collection.mutable.{ListBuffer, HashSet}

class IntervalAnnotationWalker extends RefWalker[AnnotationMetaData,List[IntervalInfoBuilder]] {
  @Argument(doc="Min proportion of bases overlapping between an interval of interest and an annotation interval for annotation to occur",shortName="mpb")
  var minPropOverlap : Double = 0.50 // default to 50%
  @Output var out : PrintStream = _

  override def reduceInit : List[IntervalInfoBuilder] = {
    Nil
  }

  override def map(tracker: RefMetaDataTracker, ref: ReferenceContext, context: AlignmentContext) : AnnotationMetaData = {
    val ivalsOfInterest : List[GATKFeature] = tracker.getGATKFeatureMetaData("interval_list",false).toList
    val refSeqData : List[Object] = tracker.getReferenceMetaData("refseq").toList
    val tableMetaData : List[Object] = tracker.getReferenceMetaData("table",false).toList
    var amd : AnnotationMetaData = new AnnotationMetaData(ref.getLocus,ref.getBase)
    if ( refSeqData != null ) {
      amd.refSeqFeatures = refSeqData.map( (o: Object) => o.asInstanceOf[RefSeqFeature])
    }
    if ( tableMetaData != null ) {
      amd.tableFeatures = tableMetaData.map( (o: Object) => o.asInstanceOf[TableFeature]).map( (u: TableFeature) =>
      new Tuple2(u.getLocation,u.getAllValues.toList.reduceLeft(_ + "," + _)))
    }
    amd.newIvals = ivalsOfInterest.map(_.getLocation)
    return amd
  }

  override def reduce(mapMetaData : AnnotationMetaData, rList : List[IntervalInfoBuilder]) : List[IntervalInfoBuilder] = {
    rList.filter( (u : IntervalInfoBuilder) => u.location.isBefore(mapMetaData.locus)).foreach(u => out.print("%s%n".format(u.done)))
    val newList : List[IntervalInfoBuilder] = rList.filter( ! _.finalized ) :::
      mapMetaData.newIvals.filter( (g: GenomeLoc) => g.getStart == mapMetaData.locus.getStart ).map( u => new IntervalInfoBuilder(u,minPropOverlap) )
    newList.foreach(_.update(mapMetaData))
    return newList
  }

  override def onTraversalDone(finalList : List[IntervalInfoBuilder]) = {
    finalList.foreach(u => out.print("%s%n".format(u.done)))
  }

}

class AnnotationMetaData(loc: GenomeLoc, ref: Byte) {
  val locus : GenomeLoc = loc
  var newIvals : List[GenomeLoc] = Nil
  var refSeqFeatures : List[RefSeqFeature] = Nil
  var tableFeatures : List[(GenomeLoc,String)] = Nil
  val refBase : Byte = ref
}

class IntervalInfoBuilder(loc : GenomeLoc, minProp : Double) {

  val OVERLAP_PROP : Double = minProp

  var metaData : HashSet[String] = new HashSet[String]
  var seenTranscripts : HashSet[String] = new HashSet[String]
  var baseContent : ListBuffer[Byte] = new ListBuffer[Byte]
  var location : GenomeLoc = loc
  var gcContent : Double = _
  var entropy : Double = _
  var finalized : Boolean = false


  def update(mmd : AnnotationMetaData) = {
    baseContent += mmd.refBase
    mmd.refSeqFeatures.filter( p => ! seenTranscripts.contains(p.getTranscriptUniqueGeneName)).view.foreach(u => addRefSeq(u))
    mmd.tableFeatures.filter( (t : (GenomeLoc,String)) =>
      ! metaData.contains(t._2) && (t._1.intersect(location).size.asInstanceOf[Double]/location.size() >= OVERLAP_PROP) ).foreach( t => metaData.add(t._2))
  }

  def addRefSeq( rs: RefSeqFeature ) : Unit = {
    seenTranscripts.add(rs.getTranscriptUniqueGeneName)
    val exons : List[(GenomeLoc,Int)] = rs.getExons.toList.zipWithIndex.filter( p => ((p._1.intersect(location).size+0.0)/location.size) >= minProp )
    if ( exons.size > 0 ) {
      metaData.add("%s[%s]".format(rs.getTranscriptUniqueGeneName,exons.map( p => "exon%d".format(p._2) ).reduceLeft(_ + "," + _)))
    } else {
      if ( location.isBefore(rs.getExons.get(0)) || location.isPast(rs.getExons.get(rs.getNumExons-1)) ) {
        metaData.add("%s[%s]".format(rs.getTranscriptUniqueGeneName,"UTR"))
      } else {
        metaData.add("%s[%s]".format(rs.getTranscriptUniqueGeneName,"Intron"))
      }
    }

  }

  def done : String = {
    finalized = true
    def isGC(b : Byte) : Int = if ( BaseUtils.gIndex == b || BaseUtils.cIndex == b ) { 1 } else { 0 }
    gcContent = baseContent.foldLeft[Int](0)( (a,b) => a + isGC(b)).asInstanceOf[Double]/location.size()
    entropy = 0.0 // todo -- implement me
    val meta : String = metaData.reduceLeft(_ + "\t" + _)
    return "%s\t%d\t%d\t%.2f\t%.2f\t%s".format(location.getContig,location.getStart,location.getStop,gcContent,entropy,meta)
  }

}