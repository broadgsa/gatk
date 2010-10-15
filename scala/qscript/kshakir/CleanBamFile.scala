import org.broadinstitute.sting.queue.extensions.firehose.ImportSingleValueFunction
import org.broadinstitute.sting.queue.extensions.picard.PicardBamJarFunction
import org.broadinstitute.sting.queue.extensions.samtools.SamtoolsIndexFunction
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._

class CleanBamFile extends QScript {
  qscript =>

  @Argument(doc="gatk jar", shortName="gatk")
  var gatkJar: File = _

  @Argument(doc="samtools binary", shortName="samtools")
  var samtoolsBinary: String = _

  @Argument(doc="fix mates jar", shortName="fixMates")
  var fixMatesJar: File = _

  @Input(doc="Script that can merge text files, for example Sting/shell/mergeText.sh.", shortName="MTS")
  var mergeTextScript: String = _

  @Argument(doc="base name for output files", shortName="base")
  var baseName: String = _

  @Input(doc="reference genome", shortName="R")
  var referenceFile: File = _

  @Input(doc="recalibrated bam", shortName="I")
  var recalibratedBam: File = _

  @Argument(doc="read group blacklist conversion script that can convert firehose outputs to a GATK blacklist file.", shortName="RGBLS")
  var readGroupBlackListScript: String = _

  @Argument(doc="read group blacklist", shortName="RGBL", required=false)
  var readGroupBlackList: String = _

  @Argument(doc="intervals", shortName="L")
  var intervals: File = _

  @Argument(doc="Script that can split the interval file by contig, for example Sting/python/splitIntervalsByContig.py.", shortName="RTCSS")
  var targetCreatorScatterScript: String = _

  @Argument(doc="RealignerTargetCreator scatter count.  " +
          "Best if it is either 1 or the number of contigs in the interval list.  " +
          "If used the compute farm must also be used.", shortName="RTCSC")
  var targetCreatorScatterCount = 0

  @Argument(doc="Script that can split the intervals evenly, for example Sting/shell/splitIntervals.sh.", shortName="IRSS")
  var indelRealignerScatterScript: String = _

  @Argument(doc="IndelRealigner scatter count.", shortName="IRSC")
  var indelRealignerScatterCount = 0

  @Input(doc="dbsnp file", shortName="D")
  var dbsnpFile: File = _

  @Argument(doc="firehose import jar", shortName="importJar")
  var firehoseImportJar: File = _

  @Argument(doc="short job queue", shortName="shortQueue", required=false)
  var shortJobQueue: String = _

  @Argument(doc="firehose host", shortName="FHHost")
  var firehoseHost: String = _

  @Argument(doc="firehose port", shortName="FHPort")
  var firehosePort: Int = _

  @Argument(doc="firehose domain", shortName="FHDomain")
  var firehoseDomain: String = _

  @Argument(doc="clean bam firehose security token", shortName="FHToken")
  var firehoseSecurityToken: String = _

  @Argument(doc="clean bam firehose entity type", shortName="bamFHEType")
  var bamFirehoseEntityType: String = _

  @Argument(doc="clean bam firehose entity id", shortName="bamFHEID")
  var bamFirehoseEntityID: String = _

  @Argument(doc="clean bam firehose annotation type name", shortName="bamFHAnn")
  var bamFirehoseAnnotationTypeName: String = _

  trait GATKCommonArgs extends CommandLineGATK {
    this.jarFile = qscript.gatkJar
    this.reference_sequence = qscript.referenceFile
    this.intervals = qscript.intervals
    this.input_file :+= recalibratedBam
  }

  def baseFile(suffix: String) = new File(baseName + suffix)

  def script = {
    val blacklistConverter = new CommandLineFunction {
      @Output(doc="blacklist file") var blacklistFile: File = _
      def commandLine = readGroupBlackListScript + " " + blacklistFile + " \"" + readGroupBlackList + "\""
    }

    if (readGroupBlackList != null) {
      blacklistConverter.blacklistFile = baseFile(".blacklist.txt")
      add(blacklistConverter)
    }

    // -T RealignerTargetCreator -I <input.bam> -R <reference.genome> <interval.list> <blacklist.file> -o <base.name>.merged.intervals
    val targetCreator = new RealignerTargetCreator with GATKCommonArgs
    targetCreator.memoryLimit = Some(2)
    targetCreator.read_group_black_list :+= blacklistConverter.blacklistFile
    targetCreator.out = baseFile(".merged.intervals")
    targetCreator.scatterCount = targetCreatorScatterCount
    targetCreator.setupScatterFunction = {
      case (scatter: IntervalScatterFunction, _) =>
        scatter.splitIntervalsScript = targetCreatorScatterScript
    }
    targetCreator.setupGatherFunction = {
      case (gather: SimpleTextGatherFunction, _) =>
        gather.mergeTextScript = mergeTextScript
    }

    // -T IndelRealigner -I <input.bam> -R <reference.genome> <blacklist.file> -stats <base.name>.indel.stats
    // -O <base.name>.unfixed.cleaned.bam -maxInRam 200000 -targetIntervals <merged.intervals> -D <dbsnp.file>
    val realigner = new IndelRealigner with GATKCommonArgs
    realigner.memoryLimit = Some(4)
    realigner.read_group_black_list :+= blacklistConverter.blacklistFile
    realigner.statisticsFileForDebugging = baseFile(".indel.stats")
    realigner.maxReadsInRam = Some(200000)
    realigner.targetIntervals = targetCreator.out
    realigner.DBSNP = dbsnpFile
    realigner.scatterCount = indelRealignerScatterCount

    var fixedBam: File = null

    if (realigner.scatterCount > 1) {
      realigner.out = baseFile(".cleaned.bam")
      // While gathering run fix mates.
      realigner.setupScatterFunction = {
        case (scatter: IntervalScatterFunction, _) =>
          scatter.splitIntervalsScript = indelRealignerScatterScript
      }
      realigner.setupGatherFunction = {
        case (gather: BamGatherFunction, _) =>
          gather.memoryLimit = Some(4)
          gather.jarFile = fixMatesJar
          // Don't pass this AS=true to fix mates!
          gather.assumeSorted = None
        case (gather: SimpleTextGatherFunction, _) =>
          gather.mergeTextScript = mergeTextScript
      }

      fixedBam = realigner.out
    } else {
      realigner.out = baseFile(".unfixed.cleaned.bam")

      // Explicitly run fix mates if the function won't be scattered.
      var fixMates = new PicardBamJarFunction {
        // Declare inputs/outputs for dependency tracking.
        @Input(doc="unfixed bam") var unfixed: File = _
        @Output(doc="fixed bam") var fixed: File = _
        def inputBams = List(unfixed)
        def outputBam = fixed
      }
      fixMates.memoryLimit = Some(4)
      fixMates.jarFile = fixMatesJar
      fixMates.unfixed = realigner.out
      fixMates.fixed = baseFile(".cleaned.bam")

      fixedBam = fixMates.fixed

      // Add the fix mates explicitly
      add(fixMates)
    }

    val bamIndex = new SamtoolsIndexFunction
    bamIndex.samtools = samtoolsBinary
    bamIndex.bamFile = fixedBam
    bamIndex.bamFileIndex = swapExt(fixedBam, "bam", "bam.bai")

    val importer = new ImportSingleValueFunction {
      /** Files that this job should wait on before running. */
      @Input(doc="Explicit job dependencies", required=false)
      var jobDependencies: List[File] = Nil
    }
    importer.jobQueue = shortJobQueue
    importer.jarFile = firehoseImportJar
    importer.host = firehoseHost
    importer.port = firehosePort
    importer.domain = firehoseDomain
    importer.securityToken = firehoseSecurityToken
    importer.entityType = bamFirehoseEntityType
    importer.entityID = bamFirehoseEntityID
    importer.annotationTypeName = bamFirehoseAnnotationTypeName
    importer.importValue = fixedBam
    importer.jobDependencies :+= bamIndex.bamFileIndex

    add(targetCreator, realigner, bamIndex, importer)
  }
}
