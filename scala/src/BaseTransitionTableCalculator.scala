package org.broadinstitute.sting.scala

import gatk.walkers.genotyper.{UnifiedGenotyper, GenotypeCall}
import java.io.File
import net.sf.samtools.SAMRecord
import org.broadinstitute.sting.gatk.walkers.LocusWalker
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker
import org.broadinstitute.sting.gatk.contexts.ReferenceContext
import org.broadinstitute.sting.gatk.contexts.AlignmentContext
import org.broadinstitute.sting.utils.Pair
import org.broadinstitute.sting.utils.genotype.GenotypeMetaData
import utils._
import cmdLine.Argument

class TransitionTable() {
  var table: Array[Array[Double]] = new Array[Array[Double]](4,4)

  def makeKey(b: Char): Int = {
    return BaseUtils.simpleBaseToBaseIndex(b)
  }

  def fromKey(i: Int): Char = {
    return BaseUtils.baseIndexToSimpleBase(i)
  }

  def add(ref: Char, subst: Char) = {
    table(makeKey(ref))(makeKey(subst)) += 1
  }

  def get(i: Int, j: Int): Double = {
    return table(i)(j)
  }

  def set(i: Int, j: Int, d: Double) = {
    //printf("Setting %d / %d => %f%n", i, j, d)
    table(i)(j) = d
  }

  def mapEntries[A](f: ((Int, Int, Double) => A)): List[A] = {
    return for { i <- List.range(0, 4); j <- List.range(0,4) } yield f(i, j, get(i,j))
  }

  def header(): String = {
    return Utils.join("\t", mapEntries((i,j,x) => "P(%c|true=%c)" format (fromKey(j), fromKey(i))).toArray)
  }

  override def toString(): String = {
    return Utils.join("\t", mapEntries((i,j,x) => "%.2f" format x).toArray)
  }

  def sum(): Double = {
     return sumList(table.toList map (x => sumList(x.toList)))
  }

  def sumList(xs: List[Double]): Double = {
      return (0.0 :: xs) reduceLeft {(x, y) => x + y}
  }

  def getRateMatrix(): TransitionTable = {
    var sums = (table.toList map (x => sumList(x.toList))).toArray
    var t = new TransitionTable()

    //println(sums.toList)
    this mapEntries ((i,j,x) => t.set(i, j, x / sums(i)))

    return t
  }
}

class BaseTransitionTableCalculator extends LocusWalker[Unit,Int] {
  private var MIN_MAPPING_QUALITY = 30
  private var MIN_BASE_QUALITY = 20
  private var MIN_LOD = 5
  private var PRINT_FREQ = 1000000
  private var CARE_ABOUT_STRAND = true

  private val SSG = new UnifiedGenotyper()
  private val table1 = new TransitionTable()
  private val tableFWD = new TransitionTable()
  private val tableREV = new TransitionTable()
  private val table2 = new TransitionTable()

  @Argument() { val shortName = "v", val doc = "Print out verbose output", val required = false }
  var VERBOSE: Boolean = false

  override def initialize() {
    SSG.initialize();
  }

  override def map(tracker: RefMetaDataTracker, ref: ReferenceContext, context: AlignmentContext): Unit = {
    val callPair = SSG.map(tracker, ref, context)
    //val call = callPair.getFirst().get(0)
    val pileup = new ReadBackedPileup(ref.getBase(), context)

    def hasNoNs(): Boolean = {
      return ! (pileup.getBases() exists ('N' ==))
    }

    if ( isDefinitelyHomRef(callPair) && hasNoNs() ) {
      var (refBases, nonRefBases) = splitBases(ref.getBase, pileup.getBases)

      //printf("nonRefBases %s%n", nonRefBases)
      if ( nonRefBases.length > 0 ) {
        var (nonRefReads, offsets) = getNonRefReads(ref.getBase, pileup)
        var nonRefRead: SAMRecord = nonRefReads.head

        nonRefBases.length match {
          case 1 =>
            addRead(nonRefRead, offsets.head, nonRefBases.head, ref, table1)
            addRead(nonRefRead, offsets.head, nonRefBases.head, ref, if ( nonRefRead.getReadNegativeStrandFlag ) tableREV else tableFWD )
          case 2 =>
            addRead(nonRefRead, offsets.head, nonRefBases.head, ref, table2)
            addRead(nonRefReads.tail.head, offsets.tail.head, nonRefBases.tail.head, ref, table2)
          case x => 
        }

        if ( VERBOSE && pileup.size > 30 && nonRefReads.length < 4 )
          printNonRefBases(ref, nonRefReads, offsets, context.getReads.size())
      }
    }
  }

  implicit def listToJavaList[T](l: List[T]) = java.util.Arrays.asList(l.toArray: _*)

  def printNonRefBases(ref: ReferenceContext, nonRefReads: List[SAMRecord], offsets: List[Integer], depth: Integer): Unit = {
    out.println(new ReadBackedPileup(ref.getLocus(), ref.getBase(), listToJavaList(nonRefReads), offsets) + " " + depth)
  }

  def addRead(nonRefRead: SAMRecord, offset: Integer, nonRefBase: Char, ref: ReferenceContext, table: TransitionTable): Unit = {
    if ( goodRead(nonRefRead, offset) ) {
      //printf("Including site:%n  call=%s%n  nonRefBases=%s%n  refBases=%s%n  nonRefRead=%s%n  offset=%s%n", call, nonRefBases, refBases, nonRefRead.format, offsets.head)
      //println(call.getReadDepth, call.getReferencebase, nonRefBases, refBases, nonRefRead.getMappingQuality)
      if ( CARE_ABOUT_STRAND && nonRefRead.getReadNegativeStrandFlag() ) {
        // it's on the reverse strand
        val refComp: Char = BaseUtils.simpleComplement(ref.getBase)
        val nonRefComp: Char = BaseUtils.simpleComplement(nonRefBase)

        //printf("Negative: %c/%c is actually %c/%c%n", ref.getBase, nonRefBases.head, refComp, nonRefComp)
        table.add(refComp, nonRefComp)
      } else {
        table.add(ref.getBase, nonRefBase)
      }
    }
  }

  def getNonRefReads(ref: Char, pileup: ReadBackedPileup): (List[SAMRecord], List[Integer]) = {
    def getBase(i: Int): Char = {
      var offset = pileup.getOffsets().get(i)
      return pileup.getReads().get(i).getReadString().charAt(offset.asInstanceOf[Int])
    }

    var subset = for { i <- List.range(0, pileup.size) if getBase(i) != ref } yield i
    var reads = subset map (i => pileup.getReads().get(i))
    var offsets = subset map (i => pileup.getOffsets().get(i))

    //println(reads, offsets, subset)
    return (reads, offsets)
  }


  def hasFewNonRefBases(refBase: Char, pileup: String): List[Char] = {
    splitBases(refBase, pileup) match {
      case (refs, nonRefs) =>
        //printf("nrf=%s, rb=%s%n", nonRefs.mkString, refs.mkString)
        return nonRefs
    }
  }

  def goodRead( read: SAMRecord, offset: Integer ): Boolean = {
    return read.getMappingQuality() > MIN_MAPPING_QUALITY && read.getBaseQualities()(offset.intValue) > MIN_BASE_QUALITY
  }

  def splitBases(refBase: Char, pileup: String): (List[Char], List[Char]) = {
    val refBases = pileup.toList filter (refBase ==)
    val nonRefBases = pileup.toList filter (refBase !=)
    return (refBases, nonRefBases)
  }
  
  def isDefinitelyHomRef(call: Pair[java.util.List[GenotypeCall],GenotypeMetaData]): Boolean = {
    if ( call == null )
       return false;

    return ! call.getFirst().get(0).isVariant() && call.getFirst().get(0).getNegLog10PError() > MIN_LOD
  }

  def reduceInit(): Int = {
    return 0;
  }

  def reduce(m: Unit, r: Int): Int = {
    return r;
  }

  override def onTraversalDone(r: Int) = {
    out.println("-------------------- FINAL RESULT --------------------")
    out.println("type\t" + table1.header)
    def print1(t: String, x: TransitionTable) = {
      out.println(t + "\t" + x)
    }

    print1("1-mismatch", table1)
    print1("2-mismatch", table2)
    print1("1-mismatch-fwd", tableFWD)
    print1("1-mismatch-rev", tableREV)
  }
}
