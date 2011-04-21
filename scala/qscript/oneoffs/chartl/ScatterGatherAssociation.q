import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import org.broadinstitute.sting.queue.extensions.gatk.{IntervalScatterFunction, CommandLineGATK}
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.utils.text.XReadLines
import collection.JavaConversions._

class ScatterGatherAssociation extends QScript {

  @Argument(fullName="gatkJar",shortName="gatk",doc="Path to the GATK jarfile",required=true)
  var gatkJar : File = _
  @Argument(fullName="metaData",shortName="SM",doc="Sample meta data",required=true)
  var metaData : File = _
  @Argument(fullName="bamList",shortName="I",doc="list of bam files (single .list file)",required=true)
  var bamList : File = _
  @Argument(fullName="outputBase",shortName="o",doc="Base for output files",required=true)
  var outBase : String = _
  @Argument(fullName="noBedGraph",shortName="nbg",doc="Don't use bedgraph format",required=false)
  var dontUseBedGraph : Boolean = false
  @Argument(fullName="reference",shortName="R",doc="Reference file, if not hg19",required=false)
  var referenceFile : File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
  @Argument(fullName="intervals",shortName="L",doc="Interval list, if not whole-exome 1.1",required=false)
  var intervalsFile : File = new File("/seq/references/HybSelOligos/whole_exome_agilent_1.1_refseq_plus_3_boosters/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list")
  @Argument(fullName="memoryLimit",shortName="M",doc="Memory limit for SG jobs",required=false)
  var memLimit : Int = 4
  @Argument(fullName="scatterJobs",shortName="SJ",doc="Number of scatter jobs",required=false)
  var scatterJobs : Int = 125

  val ASSOCIATION_TESTS = List("BaseQualityScore","InsertSizeDistribution","MappingQuality0",
    "MateMappingQuality","MateOtherContig","MateSameStrand","MateUnmapped","MismatchRate",
    "ProperPairs","ReadClipping","ReadIndels","ReadMappingQuality","ReferenceMismatches",
    "SampleDepth")

  class RegionalAssociationSG(base : String, ext : String) extends CommandLineGATK with ScatterGatherableFunction{
    this.analysis_type = "RegionalAssociation"

    @Argument(doc="useBed")
    var useBed : Boolean = true

    // the rest are output files implicitly constructed by the multiplexer

    @Output(doc="bqs")
    @Gather(classOf[SimpleTextGatherFunction])
    var bqs : File = new File(String.format("%s.%s.%s", base, "BaseQualityScore", ext))
    /*
    @Output(doc="isd")
    @Gather(classOf[SimpleTextGatherFunction])
    var isd : File = new File(String.format("%s.%s.%s",base,"InsertSizeDistribution",ext))
    @Output(doc="mq0")
    @Gather(classOf[SimpleTextGatherFunction])
    var mq0 : File = new File(String.format("%s.%s.%s",base,"MappingQuality0",ext))
    @Output(doc="mmq")
    @Gather(classOf[SimpleTextGatherFunction])
    var mmq : File = new File(String.format("%s.%s.%s",base,"MateMappingQuality",ext))
    @Output(doc="moc")
    @Gather(classOf[SimpleTextGatherFunction])
    var moc : File = new File(String.format("%s.%s.%s",base,"MateOtherContig",ext))
    @Output(doc="mss")
    @Gather(classOf[SimpleTextGatherFunction])
    var mss : File = new File(String.format("%s.%s.%s",base,"MateSameStrand",ext))
    /@Output(doc="mu")
    @Gather(classOf[SimpleTextGatherFunction])
    var mu : File = new File(String.format("%s.%s.%s",base,"MateUnmapped",ext))
    @Output(doc="mmr")
    @Gather(classOf[SimpleTextGatherFunction])
    var mmr : File = new File(String.format("%s.%s.%s",base,"MismatchRate",ext))
    @Output(doc="pp")
    @Gather(classOf[SimpleTextGatherFunction])
    var pp : File = new File(String.format("%s.%s.%s",base,"ProperPairs",ext))
    @Output(doc="rc")
    @Gather(classOf[SimpleTextGatherFunction])
    var rc : File = new File(String.format("%s.%s.%s",base,"ReadClipping",ext))
    @Output(doc="ri")
    @Gather(classOf[SimpleTextGatherFunction])
    var ri : File = new File(String.format("%s.%s.%s",base,"ReadIndels",ext))
    @Output(doc="rmq")
    @Gather(classOf[SimpleTextGatherFunction])
    var rmq : File = new File(String.format("%s.%s.%s",base,"ReadMappingQuality",ext))
    @Output(doc="rm")
    @Gather(classOf[SimpleTextGatherFunction])
    var rm : File = new File(String.format("%s.%s.%s",base,"ReferenceMismatches",ext))
    @Output(doc="sd")
    @Gather(classOf[SimpleTextGatherFunction])
    var sd : File = new File(String.format("%s.%s.%s",base,"SampleDepth",ext))
    @Output(doc="rai")
    @Gather(classOf[SimpleTextGatherFunction])
    var rli : File = new File(String.format("%s.%s.%s",base,"ReadsAberrantInsertSize",ext))
    @Output(doc="rwi")
    @Gather(classOf[SimpleTextGatherFunction])
    var rwi : File = new File(String.format("%s.%s.%s",base,"ReadsWithIndels",ext))
    */

    override def commandLine = {
      var bedStr : String = ""
      if ( useBed ) {
        bedStr = " -bg "
      }
      super.commandLine + " -AT ALL -o %s%s".format(base,bedStr)
    }
  }

  def script = {

    var ext : String = ""
    if ( dontUseBedGraph ) {
      ext = "tdf"
    } else {
      ext = "bedgraph"
    }

    var association = new RegionalAssociationSG(outBase,ext)
    association.useBed = ! dontUseBedGraph
    association.sample_metadata :+= metaData
    association.intervals :+= intervalsFile
    association.reference_sequence = referenceFile
    association.jarFile = gatkJar
    association.input_file ++= asScalaIterable((new XReadLines(bamList)).readLines).map(u => new File(u)).toList
    association.scatterCount = scatterJobs
    association.memoryLimit = Some(memLimit)
    association.scatterClass = classOf[IntervalScatterFunction]

    add(association)
  }
}
