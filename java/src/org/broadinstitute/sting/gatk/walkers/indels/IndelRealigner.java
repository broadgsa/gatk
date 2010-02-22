package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.IntervalMergingRule;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;

import java.util.*;
import java.io.File;
import java.io.FileWriter;

/**
 * Performs local realignment of reads based on misalignments due to the presence of indels.
 * Unlike most mappers, this walker uses the full alignment context to determine whether an
 * appropriate alternate reference (i.e. indel) exists and updates SAMRecords accordingly.
 */
public class IndelRealigner extends ReadWalker<Integer, Integer> {

    @Argument(fullName="targetIntervals", shortName="targetIntervals", doc="intervals file output from RealignerTargetCreator", required=true)
    protected String intervalsFile = null;

    @Argument(fullName="LODThresholdForCleaning", shortName="LOD", doc="LOD threshold above which the cleaner will clean", required=false)
    protected double LOD_THRESHOLD = 5.0;

    @Argument(fullName="entropyThreshold", shortName="entropy", doc="percentage of mismatches at a locus to be considered having high entropy", required=false)
    protected double MISMATCH_THRESHOLD = 0.15;

    @Argument(fullName="output", shortName="O", required=false, doc="Output bam (or directory if using --NwayOutput)")
    protected String baseWriterFilename = null;

    @Argument(fullName="NWayOutput", shortName="nway", required=false, doc="Should the reads be emitted in a separate bam file for each one of the input bams? [default:no]")
    protected boolean NWAY_OUTPUT = false;

    @Argument(fullName="outputSuffix", shortName="suffix", required=false, doc="Suffix to append to output bams (when using --NwayOutput) [default:'.cleaned']")
    protected String outputSuffix = "cleaned";

    @Argument(fullName="bam_compression", shortName="compress", required=false, doc="Compression level to use for output bams [default:5]")
    protected Integer compressionLevel = 5;

    public enum RealignerSortingStrategy {
        NO_SORT,
        ON_DISK,
        IN_MEMORY
    }

    @Argument(fullName="sortStrategy", shortName="sort", required=false, doc="What type of sorting strategy should we use?  Options include NO_SORT, ON_DISK, and IN_MEMORY.  Sorting in memory is much faster than on disk but should be used with care - with too much coverage or with long reads it might generate failures [default:ON_DISK]")
    protected RealignerSortingStrategy SORTING_STRATEGY = RealignerSortingStrategy.ON_DISK;

    // ADVANCED OPTIONS FOLLOW

    @Argument(fullName="outputIndels", shortName="indels", required=false, doc="Output file (text) for the indels found")
    protected String OUT_INDELS = null;

    @Argument(fullName="statisticsFile", shortName="stats", doc="print out statistics (what does or doesn't get cleaned)", required=false)
    protected String OUT_STATS = null;

    @Argument(fullName="SNPsFile", shortName="snps", doc="print out whether mismatching columns do or don't get cleaned out", required=false)
    protected String OUT_SNPS = null;

    @Argument(fullName="maxConsensuses", shortName="maxConsensuses", doc="max alternate consensuses to try (necessary to improve performance in deep coverage)", required=false)
    protected int MAX_CONSENSUSES = 30;

    @Argument(fullName="maxReadsForConsensuses", shortName="greedy", doc="max reads used for finding the alternate consensuses (necessary to improve performance in deep coverage)", required=false)
    protected int MAX_READS_FOR_CONSENSUSES = 120;

    @Argument(fullName="maxReadsForRealignment", shortName="maxReads", doc="max reads allowed at an interval for realignment", required=false)
    protected int MAX_READS = 20000;

    @Argument(fullName="writerWindowSize", shortName="writerWindowSize", doc="the window over which the writer will store reads when --sortInMemory is enabled", required=false)
    protected int SORTING_WRITER_WINDOW = 100;


    // the intervals input by the user
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;

    // the reads and known indels that fall into the current interval
    private final ReadBin readsToClean = new ReadBin();
    private final ArrayList<SAMRecord> readsNotToClean = new ArrayList<SAMRecord>();
    private final ArrayList<VariationRod> knownIndelsToTry = new ArrayList<VariationRod>();

    // the wrapper around the SAM writer
    private Map<String, SAMFileWriter> writers = null;

    // random number generator
    private static final long RANDOM_SEED = 1252863495;
    private static final Random generator = new Random(RANDOM_SEED);

    private static final int MAX_QUAL = 99;

    // fraction of mismatches that need to no longer mismatch for a column to be considered cleaned
    private static final double MISMATCH_COLUMN_CLEANED_FRACTION = 0.75;

    private static final double SW_MATCH = 30.0;      // 1.0;
    private static final double SW_MISMATCH = -10.0;  //-1.0/3.0;
    private static final double SW_GAP = -10.0;       //-1.0-1.0/3.0;
    private static final double SW_GAP_EXTEND = -2.0; //-1.0/.0;

    // other output files
    private FileWriter indelOutput = null;
    private FileWriter statsOutput = null;
    private FileWriter snpsOutput = null;

    public void initialize() {

        if ( LOD_THRESHOLD < 0.0 )
            throw new RuntimeException("LOD threshold cannot be a negative number");
        if ( MISMATCH_THRESHOLD <= 0.0 || MISMATCH_THRESHOLD > 1.0 )
            throw new RuntimeException("Entropy threshold must be a fraction between 0 and 1");

        // read in the intervals for cleaning
        List<GenomeLoc> locs = GenomeAnalysisEngine.parseIntervalRegion(Arrays.asList(intervalsFile), IntervalMergingRule.OVERLAPPING_ONLY);
        intervals = GenomeLocSortedSet.createSetFromList(locs).iterator();
        currentInterval = intervals.hasNext() ? intervals.next() : null;

        // set up the output writer(s)
        if ( baseWriterFilename != null ) {
            writers = new HashMap<String, SAMFileWriter>();
            Map<File, Set<String>> readGroupMap = getToolkit().getFileToReadGroupIdMapping();
            SAMFileWriterFactory factory = new SAMFileWriterFactory();

            if ( NWAY_OUTPUT ) {
                Map<File, SAMFileReader> readerMap = getToolkit().getFileToReaderMapping();
                for ( File file : readerMap.keySet() ) {
                    SAMFileHeader header = readerMap.get(file).getFileHeader();
                    if ( SORTING_STRATEGY == RealignerSortingStrategy.NO_SORT )
                        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
                    String newFileName = file.getName().substring(0, file.getName().length()-3) + outputSuffix + ".bam";
                    SAMFileWriter writer = factory.makeBAMWriter(header, SORTING_STRATEGY == RealignerSortingStrategy.IN_MEMORY, new File(baseWriterFilename, newFileName), compressionLevel);
                    if ( SORTING_STRATEGY == RealignerSortingStrategy.IN_MEMORY )
                        writer = new SortingSAMFileWriter(writer, SORTING_WRITER_WINDOW);
                    for ( String rg : readGroupMap.get(file) )
                        writers.put(rg, writer);
                }
            } else {
                SAMFileHeader header = getToolkit().getSAMFileHeader();
                if ( SORTING_STRATEGY == RealignerSortingStrategy.NO_SORT )
                    header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
                SAMFileWriter writer = factory.makeBAMWriter(header, SORTING_STRATEGY == RealignerSortingStrategy.IN_MEMORY, new File(baseWriterFilename), compressionLevel);
                if ( SORTING_STRATEGY == RealignerSortingStrategy.IN_MEMORY )
                    writer = new SortingSAMFileWriter(writer, SORTING_WRITER_WINDOW);
                for ( Set<String> set : readGroupMap.values() ) {
                    for ( String rg : set )
                        writers.put(rg, writer);                                     
                }
            }
        }

        if ( OUT_INDELS != null ) {
            try {
                indelOutput = new FileWriter(new File(OUT_INDELS));
            } catch (Exception e) {
                logger.warn("Failed to create output file "+ OUT_INDELS+". Indel output will be suppressed");
                err.println(e.getMessage());
                indelOutput = null;
            }
        }
        if ( OUT_STATS != null ) {
            try {
                statsOutput = new FileWriter(new File(OUT_STATS));
            } catch (Exception e) {
                logger.warn("Failed to create output file "+ OUT_STATS+". Cleaning stats output will be suppressed");
                err.println(e.getMessage());
                statsOutput = null;
            }
        }
        if ( OUT_SNPS != null ) {
            try {
                snpsOutput = new FileWriter(new File(OUT_SNPS));
            } catch (Exception e) {
                logger.warn("Failed to create output file "+ OUT_SNPS+". Cleaning snps output will be suppressed");
                err.println(e.getMessage());
                snpsOutput = null;
            }
        }
    }

    private void emit(final SAMRecord read) {
        if ( writers != null ) {
            SAMReadGroupRecord readGroup = read.getReadGroup();
            if ( readGroup == null || readGroup.getReadGroupId() == null ) {
                if ( writers.size() > 1 )
                    throw new StingException("There are multiple output writers but read " + read.toString() + " has no read group");
                writers.values().iterator().next().addAlignment(read);
            } else {
                writers.get(readGroup.getReadGroupId()).addAlignment(read);
            }
        }
    }

    private void emit(final List<SAMRecord> reads) {
        if ( writers == null )
            return;

        // break out the reads into sets for their respective writers
        Map<SAMFileWriter, Set<SAMRecord>> bins = new HashMap<SAMFileWriter, Set<SAMRecord>>();
        for ( SAMFileWriter writer : writers.values() )
            bins.put(writer, new HashSet<SAMRecord>());

        for ( SAMRecord read : reads ) {
            SAMReadGroupRecord readGroup = read.getReadGroup();
            if ( readGroup == null || readGroup.getReadGroupId() == null ) {
                if ( writers.size() > 1 )
                    throw new StingException("There are multiple output writers but read " + read.toString() + " has no read group");
                bins.get(writers.values().iterator().next()).add(read);
            } else {
                bins.get(writers.get(readGroup.getReadGroupId())).add(read);
            }
        }

        for ( Map.Entry<SAMFileWriter, Set<SAMRecord>> entry : bins.entrySet() ) {
            if ( SORTING_STRATEGY == RealignerSortingStrategy.IN_MEMORY ) {
                // we can be efficient in this case by batching the reads all together
                ((SortingSAMFileWriter)entry.getKey()).addAlignments(entry.getValue());
            } else {
                for ( SAMRecord read : entry.getValue() )
                    entry.getKey().addAlignment(read);
            }
        }
    }

    public Integer map(char[] ref, SAMRecord read) {
        if ( currentInterval == null ) {
            emit(read);
            return 0;
        }

        GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);
        // hack to get around unmapped reads having screwy locations
        if ( readLoc.getStop() == 0 )
            readLoc = GenomeLocParser.createGenomeLoc(readLoc.getContig(), readLoc.getStart(), readLoc.getStart());

        if ( readLoc.isBefore(currentInterval) || Utils.is454Read(read) ) {
            // TODO -- it would be nice if we could use indels from 454 reads as alternate consenses
            emit(read);
            return 0;
        }
        else if ( readLoc.overlapsP(currentInterval) ) {
            if ( read.getReadUnmappedFlag() ||
                 read.getDuplicateReadFlag() ||
                 read.getNotPrimaryAlignmentFlag() ||
                 read.getMappingQuality() == 0 ||
                 read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ) {
                readsNotToClean.add(read);
            } else {
                readsToClean.add(read, ref);
            }

            if ( readsToClean.size() + readsNotToClean.size() >= MAX_READS ) {
                // merge the two sets for emission
                readsNotToClean.addAll(readsToClean.getReads());
                emit(readsNotToClean);
                readsToClean.clear();
                readsNotToClean.clear();
                currentInterval = intervals.hasNext() ? intervals.next() : null;
            }
        }
        else {  // the read is past the current interval

            clean(readsToClean);
            knownIndelsToTry.clear();

            // merge the two sets for emission
            readsNotToClean.addAll(readsToClean.getReads());
            emit(readsNotToClean);
            readsToClean.clear();
            readsNotToClean.clear();

            do {
                currentInterval = intervals.hasNext() ? intervals.next() : null;
            } while ( currentInterval != null && currentInterval.isBefore(readLoc) );

            // call back into map now that the state has been updated
            map(ref, read);
        }

        return 0;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        if ( readsToClean.size() > 0 || readsNotToClean.size() > 0 ) {
            clean(readsToClean);
            knownIndelsToTry.clear();

            // merge the two sets for emission
            readsNotToClean.addAll(readsToClean.getReads());
            emit(readsNotToClean);
        }

        if ( writers != null ) {
            HashSet<SAMFileWriter> uniqueWriters = new HashSet<SAMFileWriter>(writers.values());
            for ( SAMFileWriter writer : uniqueWriters ) {
                writer.close();
                if ( SORTING_STRATEGY == RealignerSortingStrategy.IN_MEMORY )
                    ((SortingSAMFileWriter)writer).getBaseWriter().close();
            }
        }

        if ( OUT_INDELS != null ) {
            try {
                indelOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_INDELS+" gracefully. Data may be corrupt.");
            }
        }
        if ( OUT_STATS != null ) {
            try {
                statsOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_STATS+" gracefully. Data may be corrupt.");
            }
        }
        if ( OUT_SNPS != null ) {
            try {
                snpsOutput.close();
            } catch (Exception e) {
                logger.error("Failed to close "+OUT_SNPS+" gracefully. Data may be corrupt.");
            }
        }
    }

    private static int mismatchQualitySumIgnoreCigar(final AlignedRead aRead, final byte[] refSeq, int refIndex, int quitAboveThisValue) {
        final byte[] readSeq = aRead.getRead().getReadBases();
        final byte[] quals = aRead.getRead().getBaseQualities();
        int sum = 0;
        for (int readIndex = 0 ; readIndex < readSeq.length ; refIndex++, readIndex++ ) {
            if ( refIndex >= refSeq.length ) {
                sum += MAX_QUAL;
                // optimization: once we pass the threshold, stop calculating
                if ( sum > quitAboveThisValue )
                    return sum;
            } else {
                byte refChr = refSeq[refIndex];
                byte readChr = readSeq[readIndex];
                if ( !BaseUtils.isRegularBase((char)readChr) || !BaseUtils.isRegularBase((char)refChr) )
                    continue; // do not count Ns/Xs/etc ?
                if ( readChr != refChr ) {
                    sum += (int)quals[readIndex];
                    // optimization: once we pass the threshold, stop calculating
                    if ( sum > quitAboveThisValue )
                        return sum;
                }
            }
        }
        return sum;
    }

    private static boolean readIsClipped(final SAMRecord read) {
        final Cigar c = read.getCigar();
        final int n = c.numCigarElements();
        if ( c.getCigarElement(n-1).getOperator() == CigarOperator.S ||
             c.getCigarElement(0).getOperator() == CigarOperator.S) return true;
        return false;
    }

    private void clean(ReadBin readsToClean) {

        final List<SAMRecord> reads = readsToClean.getReads();
        if ( reads.size() == 0 )
            return;

        final byte[] reference = readsToClean.getRereference();
        final long leftmostIndex = readsToClean.getLocation().getStart();

        final ArrayList<SAMRecord> refReads = new ArrayList<SAMRecord>();                   // reads that perfectly match ref
        final ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();               // reads that don't perfectly match
        final LinkedList<AlignedRead> altAlignmentsToTest = new LinkedList<AlignedRead>();  // should we try to make an alt consensus from the read?
        final ArrayList<AlignedRead> leftMovedIndels = new ArrayList<AlignedRead>();
        final Set<Consensus> altConsenses = new LinkedHashSet<Consensus>();                 // list of alt consenses
        int totalMismatchSum = 0;

        // if there are any known indels for this region, get them
        for ( VariationRod knownIndel : knownIndelsToTry ) {
            String indelStr = knownIndel.isInsertion() ? knownIndel.getAlternateAlleleList().get(0) : Utils.dupString('-', knownIndel.getAlleleList().get(0).length());
            Consensus c = createAlternateConsensus((int)(knownIndel.getLocation().getStart() - leftmostIndex), reference, indelStr, knownIndel.isDeletion());
            if ( c != null )
                altConsenses.add(c);
        }

        // decide which reads potentially need to be cleaned
        for ( final SAMRecord read : reads ) {

            //            if ( debugOn ) {
            //                System.out.println(read.getReadName()+" "+read.getCigarString()+" "+read.getAlignmentStart()+"-"+read.getAlignmentEnd());
            //                System.out.println(reference.substring((int)(read.getAlignmentStart()-leftmostIndex),(int)(read.getAlignmentEnd()-leftmostIndex)));
            //                System.out.println(read.getReadString());
            //            }

            // we currently can not deal with clipped reads correctly (or screwy record)
            if ( read.getCigar().numCigarElements() == 0 || readIsClipped(read) ) {
                refReads.add(read);
                continue;
            }

            final AlignedRead aRead = new AlignedRead(read);

            // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
            int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
            if ( numBlocks == 2 ) {
                Cigar newCigar = indelRealignment(read.getCigar(), reference, read.getReadBases(), read.getAlignmentStart()-(int)leftmostIndex, 0);
                if ( aRead.setCigar(newCigar) ) {
                    leftMovedIndels.add(aRead);
                }
            }

            final int mismatchScore = mismatchQualitySumIgnoreCigar(aRead, reference, read.getAlignmentStart()-(int)leftmostIndex, Integer.MAX_VALUE);
            //            if ( debugOn ) System.out.println("mismatchScore="+mismatchScore);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( mismatchScore > 0 ) {
                altReads.add(aRead);
                if ( !read.getDuplicateReadFlag() )
                    totalMismatchSum += mismatchScore;
                aRead.setMismatchScoreToReference(mismatchScore);
                // if it has an indel, let's see if that's the best consensus
                if ( numBlocks == 2 )  {
                    Consensus c = createAlternateConsensus(aRead.getAlignmentStart() - (int)leftmostIndex, aRead.getCigar(), reference, aRead.getRead().getReadBases());
                    if ( c == null ) {} //System.out.println("ERROR: Failed to create alt consensus for read "+aRead.getRead().getReadName());
                    else altConsenses.add(c);
                }
                else {
                    //                    if ( debugOn ) System.out.println("Going to test...");
                    altAlignmentsToTest.add(aRead);
                }
            }
            // otherwise, we can emit it as is
            else {
                // if ( debugOn ) System.out.println("Emitting as is...");
                refReads.add(read);
            }
        }

        // choose alternate consensuses
        if ( altAlignmentsToTest.size() <= MAX_READS_FOR_CONSENSUSES ) {
            for ( AlignedRead aRead : altAlignmentsToTest ) {
                // do a pairwise alignment against the reference
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getRead().getReadBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
                Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, aRead.getRead().getReadBases());
                if ( c != null ) {
                    //                    if ( debugOn ) System.out.println("NEW consensus generated by SW: "+c.str ) ;
                    altConsenses.add(c);
                } else {
                    //   if ( debugOn ) System.out.println("FAILED to create Alt consensus from SW");
                }
            }
        } else {
            // choose alternate consenses randomly
            int readsSeen = 0;
            while ( readsSeen++ < MAX_READS_FOR_CONSENSUSES && altConsenses.size() <= MAX_CONSENSUSES) {
                int index = generator.nextInt(altAlignmentsToTest.size());
                AlignedRead aRead = altAlignmentsToTest.remove(index);
                // do a pairwise alignment against the reference
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getRead().getReadBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
                Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, aRead.getRead().getReadBases());
                if ( c != null )
                    altConsenses.add(c);
            }
        }

        Consensus bestConsensus = null;
        Iterator<Consensus> iter = altConsenses.iterator();

        // if ( debugOn ) System.out.println("------\nChecking consenses...\n--------\n");

        while ( iter.hasNext() ) {
            Consensus consensus = iter.next();

            // if ( debugOn ) System.out.println("Consensus: "+consensus.str);

            for ( int j = 0; j < altReads.size(); j++ ) {
                AlignedRead toTest = altReads.get(j);
                Pair<Integer, Integer> altAlignment = findBestOffset(consensus.str, toTest, (int)leftmostIndex);

                // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                int myScore = altAlignment.second;
                if ( myScore >= toTest.getMismatchScoreToReference() )
                    myScore = toTest.getMismatchScoreToReference();
                // keep track of reads that align better to the alternate consensus.
                // By pushing alignments with equal scores to the alternate, it means we'll over-call (het -> hom non ref) but are less likely to under-call (het -> ref, het non ref -> het)
                else
                    consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.first));

                //logger.debug(consensus.str +  " vs. " + toTest.getRead().getReadString() + " => " + myScore + " - " + altAlignment.first);
                if ( !toTest.getRead().getDuplicateReadFlag() )
                    consensus.mismatchSum += myScore;

                // optimization: once the mismatch sum is higher than the best consensus, quit since this one can't win
                //  THIS MUST BE DISABLED IF WE DECIDE TO ALLOW MORE THAN ONE ALTERNATE CONSENSUS!
                if ( bestConsensus != null && consensus.mismatchSum > bestConsensus.mismatchSum )
                    break;
            }

            //logger.debug(consensus.str +  " " + consensus.mismatchSum);
            if ( bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                // we do not need this alt consensus, release memory right away!!
                if ( bestConsensus != null )
                    bestConsensus.readIndexes.clear();
                bestConsensus = consensus;
                //logger.debug(consensus.str +  " " + consensus.mismatchSum);
            } else {
                // we do not need this alt consensus, release memory right away!!
                consensus.readIndexes.clear();
            }
        }

        // if the best alternate consensus has a smaller sum of quality score mismatches (more than
        // the LOD threshold), and it didn't just move around the mismatching columns, then clean!
        final double improvement = (bestConsensus == null ? -1 : ((double)(totalMismatchSum - bestConsensus.mismatchSum))/10.0);
        if ( improvement >= LOD_THRESHOLD ) {

            bestConsensus.cigar = indelRealignment(bestConsensus.cigar, reference, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference);

           // start cleaning the appropriate reads
            for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes ) {
                AlignedRead aRead = altReads.get(indexPair.first);
                updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.second, aRead, (int)leftmostIndex);
            }
            if ( !alternateReducesEntropy(altReads, reference, leftmostIndex) ) {
                if ( statsOutput != null ) {
                    try {
                        statsOutput.write(readsToClean.getLocation().toString());
                        statsOutput.write("\tFAIL (bad indel)\t"); // if improvement > LOD_THRESHOLD *BUT* entropy is not reduced (SNPs still exist)
                        statsOutput.write(Double.toString(improvement));
                        statsOutput.write("\n");
                        statsOutput.flush();
                    } catch (Exception e) {}
                }
            } else {
                //logger.debug("CLEAN: " + AlignmentUtils.cigarToString(bestConsensus.cigar) + " " + bestConsensus.str );
                if ( indelOutput != null && bestConsensus.cigar.numCigarElements() > 1 ) {
                    // NOTE: indels are printed out in the format specified for the low-coverage pilot1
                    //  indel calls (tab-delimited): chr position size type sequence
                    StringBuilder str = new StringBuilder();
                    str.append(reads.get(0).getReferenceName());
                    int position = bestConsensus.positionOnReference + bestConsensus.cigar.getCigarElement(0).getLength();
                    str.append("\t" + (leftmostIndex + position - 1));
                    CigarElement ce = bestConsensus.cigar.getCigarElement(1);
                    str.append("\t" + ce.getLength() + "\t" + ce.getOperator() + "\t");
                    int length = ce.getLength();
                    if ( ce.getOperator() == CigarOperator.D ) {
                        for ( int i = 0; i < length; i++)
                            str.append((char)reference[position+i]);
                    } else {
                        for ( int i = 0; i < length; i++)
                            str.append((char)bestConsensus.str[position+i]);
                    }
                    str.append("\t" + (((double)(totalMismatchSum - bestConsensus.mismatchSum))/10.0) + "\n");
                    try {
                        indelOutput.write(str.toString());
                        indelOutput.flush();
                    } catch (Exception e) {}
                }
                if ( statsOutput != null ) {
                    try {
                        statsOutput.write(readsToClean.getLocation().toString());
                        statsOutput.write("\tCLEAN"); // if improvement > LOD_THRESHOLD *AND* entropy is reduced
                        if ( bestConsensus.cigar.numCigarElements() > 1 )
                            statsOutput.write(" (found indel)");
                        statsOutput.write("\t");
                        statsOutput.write(Double.toString(improvement));
                        statsOutput.write("\n");
                        statsOutput.flush();
                    } catch (Exception e) {}
                }

                // We need to update the mapping quality score of the cleaned reads;
                // however we don't have enough info to use the proper MAQ scoring system.
                // For now, we'll use a heuristic:
                // the mapping quality score is improved by the LOD difference in mismatching
                // bases between the reference and alternate consensus (divided by 10)

                // finish cleaning the appropriate reads
                for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes ) {
                    final AlignedRead aRead = altReads.get(indexPair.first);
                    if ( aRead.finalizeUpdate() ) {
                        aRead.getRead().setMappingQuality(Math.min(aRead.getRead().getMappingQuality() + (int)(improvement/10.0), 255));
                        aRead.getRead().setAttribute("NM", AlignmentUtils.numMismatches(aRead.getRead(), reference, aRead.getRead().getAlignmentStart()-(int)leftmostIndex));
                    }
                }
            }

            // END IF ( improvement >= LOD_THRESHOLD )

        } else if ( statsOutput != null ) {
            try {
                statsOutput.write(readsToClean.getLocation().toString());
                statsOutput.write("\tFAIL\t"); // if improvement < LOD_THRESHOLD
                statsOutput.write(Double.toString(improvement));
                statsOutput.write("\n");
                statsOutput.flush();
            } catch (Exception e) {}
        }
    }

    private Consensus createAlternateConsensus(final int indexOnRef, final Cigar c, final byte[] reference, final byte[] readStr) {
        if ( indexOnRef < 0 )
            return null;

        // create the new consensus
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < indexOnRef; i++)
            sb.append((char)reference[i]);
        //logger.debug("CIGAR = " + AlignmentUtils.cigarToString(c));

        int indelCount = 0;
        int altIdx = 0;
        int refIdx = indexOnRef;
        boolean ok_flag = true;
        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
            CigarElement ce = c.getCigarElement(i);
            int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
            case D:
                indelCount++;
                refIdx += elementLength;
                break;
            case M:
                if ( reference.length < refIdx + elementLength )
                    ok_flag = false;
                else  {
                    for (int j = 0; j < elementLength; j++)
                        sb.append((char)reference[refIdx+j]);
                }
                refIdx += elementLength;
                altIdx += elementLength;
                break;
            case I:
                for (int j = 0; j < elementLength; j++)
                    sb.append((char)readStr[altIdx + j]);
                altIdx += elementLength;
                indelCount++;
                break;
            }
        }
        // make sure that there is at most only a single indel and it aligns appropriately!
        if ( !ok_flag || indelCount != 1 || reference.length < refIdx )
            return null;

        for (int i = refIdx; i < reference.length; i++)
            sb.append((char)reference[i]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, c, indexOnRef);
    }

    private Consensus createAlternateConsensus(final int indexOnRef, final byte[] reference, final String indelStr, final boolean isDeletion) {
        if ( indexOnRef < 0 )
            return null;

        // create the new consensus
        StringBuilder sb = new StringBuilder();
        Cigar cigar = new Cigar();
        int refIdx;

        for (refIdx = 0; refIdx < indexOnRef; refIdx++)
            sb.append((char)reference[refIdx]);
        if ( indexOnRef > 0 )
            cigar.add(new CigarElement(indexOnRef, CigarOperator.M));

        if ( isDeletion ) {
            refIdx += indelStr.length();
            cigar.add(new CigarElement(indelStr.length(), CigarOperator.D));
        }
        else {
            sb.append(indelStr);
            cigar.add(new CigarElement(indelStr.length(), CigarOperator.I));
        }

        if ( reference.length - refIdx > 0 )
            cigar.add(new CigarElement(reference.length - refIdx, CigarOperator.M));        
        for (; refIdx < reference.length; refIdx++)
            sb.append((char)reference[refIdx]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, cigar, indexOnRef);
    }

    private Pair<Integer, Integer> findBestOffset(final byte[] ref, final AlignedRead read, final int leftmostIndex) {

        // optimization: try the most likely alignment first (to get a low score to beat)
        int originalAlignment = read.getOriginalAlignmentStart() - leftmostIndex;
        int bestScore = mismatchQualitySumIgnoreCigar(read, ref, originalAlignment, Integer.MAX_VALUE);
        int bestIndex = originalAlignment;

        // optimization: we can't get better than 0, so we can quit now
        if ( bestScore == 0 )
            return new Pair<Integer, Integer>(bestIndex, 0);

        // optimization: the correct alignment shouldn't be too far from the original one (or else the read wouldn't have aligned in the first place)
        for ( int i = 0; i < originalAlignment; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        final int maxPossibleStart = ref.length - read.getReadLength();
        for ( int i = originalAlignment + 1; i <= maxPossibleStart; i++ ) {
            int score = mismatchQualitySumIgnoreCigar(read, ref, i, bestScore);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
            // optimization: we can't get better than 0, so we can quit now
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
        }

        return new Pair<Integer, Integer>(bestIndex, bestScore);
    }


    private void updateRead(final Cigar altCigar, final int altPosOnRef, final int myPosOnAlt, final AlignedRead aRead, final int leftmostIndex) {
        Cigar readCigar = new Cigar();

        // special case: there is no indel
        if ( altCigar.getCigarElements().size() == 1 ) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.setCigar(readCigar);
            return;
        }

        CigarElement altCE1 = altCigar.getCigarElement(0);
        CigarElement altCE2 = altCigar.getCigarElement(1);

        int leadingMatchingBlockLength = 0; // length of the leading M element or 0 if the leading element is I

        CigarElement indelCE;
        if ( altCE1.getOperator() == CigarOperator.I  ) {
            indelCE=altCE1;
            if ( altCE2.getOperator() != CigarOperator.M  )
                throw new StingException("When first element of the alt consensus is I, the second one must be M. Actual: "+altCigar.toString());
        }
        else {
            if ( altCE1.getOperator() != CigarOperator.M  )
                throw new StingException("First element of the alt consensus cigar must be M or I. Actual: "+altCigar.toString());
            if ( altCE2.getOperator() == CigarOperator.I  || altCE2.getOperator() == CigarOperator.D ) indelCE=altCE2;
            else
                throw new StingException("When first element of the alt consensus is M, the second one must be I or D. Actual: "+altCigar.toString());
            leadingMatchingBlockLength = altCE1.getLength();
        }

        // the easiest thing to do is to take each case separately
        int endOfFirstBlock = altPosOnRef + leadingMatchingBlockLength;
        boolean sawAlignmentStart = false;

        // for reads starting before the indel
        if ( myPosOnAlt < endOfFirstBlock) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt);
            sawAlignmentStart = true;

            // for reads ending before the indel
            if ( myPosOnAlt + aRead.getReadLength() <= endOfFirstBlock) {
                readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
                aRead.setCigar(readCigar);
                return;
            }
            readCigar.add(new CigarElement(endOfFirstBlock - myPosOnAlt, CigarOperator.M));
        }

        int indelOffsetOnRef = 0, indelOffsetOnRead = 0;
        // forward along the indel
        if ( indelCE.getOperator() == CigarOperator.I ) {
            // for reads that end in an insertion
            if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + indelCE.getLength() ) {
                readCigar.add(new CigarElement(myPosOnAlt + aRead.getReadLength() - endOfFirstBlock, CigarOperator.I));
                aRead.setCigar(readCigar);
                return;
            }

            // for reads that start in an insertion
            if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + indelCE.getLength() ) {
                aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(indelCE.getLength() - (myPosOnAlt - endOfFirstBlock), CigarOperator.I));
                indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
                sawAlignmentStart = true;
            } else if ( sawAlignmentStart ) {
                readCigar.add(indelCE);
                indelOffsetOnRead = indelCE.getLength();
            }
        } else if ( indelCE.getOperator() == CigarOperator.D ) {
            if ( sawAlignmentStart )
                readCigar.add(indelCE);
            indelOffsetOnRef = indelCE.getLength();
        }

        // for reads that start after the indel
        if ( !sawAlignmentStart ) {
            aRead.setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.setCigar(readCigar);
            return;
        }

        int readRemaining = aRead.getReadLength();
        for ( CigarElement ce : readCigar.getCigarElements() ) {
            if ( ce.getOperator() != CigarOperator.D )
                readRemaining -= ce.getLength();
        }
        if ( readRemaining > 0 )
            readCigar.add(new CigarElement(readRemaining, CigarOperator.M));
        aRead.setCigar(readCigar);
    }

    private boolean alternateReducesEntropy(final List<AlignedRead> reads, final byte[] reference, final long leftmostIndex) {
        final int[] originalMismatchBases = new int[reference.length];
        final int[] cleanedMismatchBases = new int[reference.length];
        final int[] totalOriginalBases = new int[reference.length];
        final int[] totalCleanedBases = new int[reference.length];

        // set to 1 to prevent dividing by zero
        for ( int i=0; i < reference.length; i++ )
            originalMismatchBases[i] = totalOriginalBases[i] = cleanedMismatchBases[i] = totalCleanedBases[i] = 1;

        for (int i=0; i < reads.size(); i++) {
            final AlignedRead read = reads.get(i);
            if ( read.getRead().getAlignmentBlocks().size() > 1 )
                 continue;

            int refIdx = read.getOriginalAlignmentStart() - (int)leftmostIndex;
            final byte[] readStr = read.getRead().getReadBases();
            final byte[] quals = read.getRead().getBaseQualities();

            for (int j=0; j < readStr.length; j++, refIdx++ ) {
                if ( refIdx < 0 || refIdx >= reference.length ) {
                    //System.out.println( "Read: "+read.getRead().getReadName() + "; length = " + readStr.length() );
                    //System.out.println( "Ref left: "+ leftmostIndex +"; ref length=" + reference.length() + "; read alignment start: "+read.getOriginalAlignmentStart() );
                    break;
                }
                totalOriginalBases[refIdx] += quals[j];
                if ( readStr[j] != reference[refIdx] )
                    originalMismatchBases[refIdx] += quals[j];
            }

            // reset and now do the calculation based on the cleaning
            refIdx = read.getAlignmentStart() - (int)leftmostIndex;
            int altIdx = 0;
            Cigar c = read.getCigar();
            for (int j = 0 ; j < c.numCigarElements() ; j++) {
                CigarElement ce = c.getCigarElement(j);
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case M:
                        for (int k = 0 ; k < elementLength ; k++, refIdx++, altIdx++ ) {
                            if ( refIdx >= reference.length )
                                break;
                            totalCleanedBases[refIdx] += quals[altIdx];
                            if ( readStr[altIdx] != reference[refIdx] )
                                cleanedMismatchBases[refIdx] += quals[altIdx];
                        }
                        break;
                    case I:
                        altIdx += elementLength;
                        break;
                    case D:
                        refIdx += elementLength;
                        break;
                }

            }
        }

        int originalMismatchColumns = 0, cleanedMismatchColumns = 0;
        StringBuilder sb = new StringBuilder();
        for ( int i=0; i < reference.length; i++ ) {
            if ( cleanedMismatchBases[i] == originalMismatchBases[i] )
                continue;
            boolean didMismatch = false, stillMismatches = false;
            if ( originalMismatchBases[i] > totalOriginalBases[i] * MISMATCH_THRESHOLD )  {
                didMismatch = true;
                originalMismatchColumns++;
                if ( ((double)cleanedMismatchBases[i] / (double)totalCleanedBases[i]) > ((double)originalMismatchBases[i] / (double)totalOriginalBases[i]) * (1.0 - MISMATCH_COLUMN_CLEANED_FRACTION) ) {
                    stillMismatches = true;
                    cleanedMismatchColumns++;
                }
            } else if ( cleanedMismatchBases[i] > totalCleanedBases[i] * MISMATCH_THRESHOLD ) {
                cleanedMismatchColumns++;
            }
            if ( snpsOutput != null ) {
                    if ( didMismatch ) {
                        sb.append(reads.get(0).getRead().getReferenceName() + ":");
                        sb.append(((int)leftmostIndex + i));
                        if ( stillMismatches )
                            sb.append(" SAME_SNP\n");
                        else
                            sb.append(" NOT_SNP\n");
                    }
            }
        }

        //logger.debug("Original mismatch columns = " + originalMismatchColumns + "; cleaned mismatch columns = " + cleanedMismatchColumns);

        final boolean reduces = (originalMismatchColumns == 0 || cleanedMismatchColumns < originalMismatchColumns);
        if ( reduces && snpsOutput != null ) {
            try {
                snpsOutput.write(sb.toString());
                snpsOutput.flush();
            } catch (Exception e) {}
        }
        return reduces;
    }

    /** Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
     * starting at 0-based position <code>refIndex</code> on the <code>refSeq</code> and specified by its <code>cigar</code>.
     * The last argument <code>readIndex</code> specifies 0-based position on the read where the alignment described by the
     * <code>cigar</code> starts. Usually cigars specify alignments of the whole read to the ref, so that readIndex is normally 0.
     * Use non-zero readIndex only when the alignment cigar represents alignment of a part of the read. The refIndex in this case
     * should be the position where the alignment of that part of the read starts at. In other words, both refIndex and readIndex are
     * always the positions where the cigar starts on the ref and on the read, respectively.
     *
     * If the alignment has an indel, then this method attempts moving this indel left across a stretch of repetitive bases. For instance, if the original cigar
     * specifies that (any) one AT  is deleted from a repeat sequence TATATATA, the output cigar will always mark the leftmost AT
     * as deleted. If there is no indel in the original cigar, or the indel position is determined unambiguously (i.e. inserted/deleted sequence
     * is not repeated), the original cigar is returned.
     * @param cigar structure of the original alignment
     * @param refSeq reference sequence the read is aligned to
     * @param readSeq read sequence
     * @param refIndex 0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @return a cigar, in which indel is guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    private Cigar indelRealignment(Cigar cigar, final byte[] refSeq, final byte[] readSeq, final int refIndex, final int readIndex) {
        if ( cigar.numCigarElements() < 2 ) return cigar; // no indels, nothing to do

        final CigarElement ce1 = cigar.getCigarElement(0);
        final CigarElement ce2 = cigar.getCigarElement(1);

        // we currently can not handle clipped reads; alternatively, if the alignment starts from insertion, there
        // is no place on the read to move that insertion further left; so we are done:
        if ( ce1.getOperator() != CigarOperator.M ) return cigar;

        int difference = 0; // we can move indel 'difference' bases left
        final int indel_length = ce2.getLength();

        int period = 0; // period of the inserted/deleted sequence
        int indelIndexOnRef = refIndex+ce1.getLength() ; // position of the indel on the REF (first deleted base or first base after insertion)
        int indelIndexOnRead = readIndex+ce1.getLength(); // position of the indel on the READ (first insterted base, of first base after deletion)

        byte[] indelString = new byte[ce2.getLength()];  // inserted or deleted sequence
        if ( ce2.getOperator() == CigarOperator.D )
            System.arraycopy(refSeq, indelIndexOnRef, indelString, 0, ce2.getLength());
        else if ( ce2.getOperator() == CigarOperator.I )
            System.arraycopy(readSeq, indelIndexOnRead, indelString, 0, ce2.getLength());
        else
            // we can get here if there is soft clipping done at the beginning of the read
            // for now, we'll just punt the issue and not try to realign these
            return cigar;

        // now we have to check all WHOLE periods of the indel sequence:
        //  for instance, if
        //   REF:   AGCTATATATAGCC
        //   READ:   GCTAT***TAGCC
        // the deleted sequence ATA does have period of 2, but deletion obviously can not be
        // shifted left by 2 bases (length 3 does not contain whole number of periods of 2);
        // however if 4 bases are deleted:
        //   REF:   AGCTATATATAGCC
        //   READ:   GCTA****TAGCC
        // the length 4 is a multiple of the period of 2, and indeed deletion site can be moved left by 2 bases!
        //  Also, we will always have to check the length of the indel sequence itself (trivial period). If the smallest
        // period is 1 (which means that indel sequence is a homo-nucleotide sequence), we obviously do not have to check
        // any other periods.

        // NOTE: we treat both insertions and deletions in the same way below: we always check if the indel sequence
        // repeats itsels on the REF (never on the read!), even for insertions: if we see TA inserted and REF has, e.g., CATATA prior to the insertion
        // position, we will move insertion left, to the position right after CA. This way, while moving the indel across the repeat
        // on the ref, we can theoretically move it across a non-repeat on the read if the latter has a mismtach.

        while ( period < indel_length ) { // we will always get at least trivial period = indelStringLength

                period = BaseUtils.sequencePeriod(indelString, period+1);

                if ( indel_length % period != 0 ) continue; // if indel sequence length is not a multiple of the period, it's not gonna work

                int newIndex = indelIndexOnRef;

                while ( newIndex >= period ) { // let's see if there is a repeat, i.e. if we could also say that same bases at lower position are deleted

                    // lets check if bases [newIndex-period,newIndex) immediately preceding the indel on the ref
                    // are the same as the currently checked period of the inserted sequence:

                    boolean match = true;

                    for ( int testRefPos = newIndex - period, indelPos = 0 ; testRefPos < newIndex; testRefPos++, indelPos++) {
                        byte indelChr = indelString[indelPos];
                        if ( refSeq[testRefPos] != indelChr || !BaseUtils.isRegularBase((char)indelChr) ) {
                            match = false;
                            break;
                        }
                    }
                    if ( match )
                        newIndex -= period; // yes, they are the same, we can move indel farther left by at least period bases, go check if we can do more...
                    else break; // oops, no match, can not push indel farther left
                }

                final int newDifference = indelIndexOnRef - newIndex;
                if ( newDifference > difference ) difference = newDifference; // deletion should be moved 'difference' bases left

                if ( period == 1 ) break; // we do not have to check all periods of homonucleotide sequences, we already
                                          // got maximum possible shift after checking period=1 above.
        }

        //        if ( ce2.getLength() >= 2 )
        //            System.out.println("-----------------------------------\n  FROM:\n"+AlignmentUtils.alignmentToString(cigar,readSeq,refSeq,refIndex, (readIsConsensusSequence?refIndex:0)));


        if ( difference > 0 ) {

            // The following if() statement: this should've never happened, unless the alignment is really screwed up.
            // A real life example:
            //
            //   ref:    TTTTTTTTTTTTTTTTTT******TTTTTACTTATAGAAGAAAT...
            //  read:       GTCTTTTTTTTTTTTTTTTTTTTTTTACTTATAGAAGAAAT...
            //
            //  i.e. the alignment claims 6 T's to be inserted. The alignment is clearly malformed/non-conforming since we could
            // have just 3 T's inserted (so that the beginning of the read maps right onto the beginning of the
            // reference fragment shown): that would leave us with same 2 mismatches at the beginning of the read
            // (G and C) but lower gap penalty. Note that this has nothing to do with the alignment being "right" or "wrong"
            // with respect to where on the DNA the read actually came from. It is the assumptions of *how* the alignments are
            // built and represented that are broken here. While it is unclear how the alignment shown above could be generated
            // in the first place, we are not in the business of fixing incorrect alignments in this method; all we are
            // trying to do is to left-adjust correct ones. So if something like that happens, we refuse to change the cigar
            // and bail out.
            if ( ce1.getLength()-difference < 0 ) return cigar;

            Cigar newCigar = new Cigar();
            // do not add leading M cigar element if its length is zero (i.e. if we managed to left-shift the
            // insertion all the way to the read start):
            if ( ce1.getLength() - difference > 0 )
                newCigar.add(new CigarElement(ce1.getLength()-difference, CigarOperator.M));
            newCigar.add(ce2);  // add the indel, now it's left shifted since we decreased the number of preceding matching bases

            if ( cigar.numCigarElements() > 2 ) {
                // if we got something following the indel element:

                if ( cigar.getCigarElement(2).getOperator() == CigarOperator.M  ) {
                    // if indel was followed by matching bases (that's the most common situation),
                    // increase the length of the matching section after the indel by the amount of left shift
                    // (matching bases that were on the left are now *after* the indel; we have also checked at the beginning
                    // that the first cigar element was also M):
                    newCigar.add(new CigarElement(cigar.getCigarElement(2).getLength()+difference, CigarOperator.M));
                } else {
                    // if the element after the indel was not M, we have to add just the matching bases that were on the left
                    // and now appear after the indel after we performed the shift. Then add the original element that followed the indel.
                    newCigar.add(new CigarElement(difference, CigarOperator.M));
                    newCigar.add(new CigarElement(cigar.getCigarElement(2).getLength(),cigar.getCigarElement(2).getOperator()));
                }
                // now add remaining (unchanged) cigar elements, if any:
                for ( int i = 3 ; i < cigar.numCigarElements() ; i++ )  {
                    newCigar.add(new CigarElement(cigar.getCigarElement(i).getLength(),cigar.getCigarElement(i).getOperator()));
                }
            }

            //logger.debug("Realigning indel: " + AlignmentUtils.cigarToString(cigar) + " to " + AlignmentUtils.cigarToString(newCigar));
            cigar = newCigar;

        }
        return cigar;
    }

    private class AlignedRead {
        private final SAMRecord read;
        private Cigar newCigar = null;
        private int newStart = -1;
        private int mismatchScoreToReference;

        public AlignedRead(SAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public SAMRecord getRead() {
               return read;
        }

        public int getReadLength() {
               return read.getReadLength();
        }

        public Cigar getCigar() {
            return (newCigar != null ? newCigar : read.getCigar());
        }

        // tentatively sets the new Cigar, but it needs to be confirmed later
        // returns true if the new cigar is a valid change (i.e. not same as original and doesn't remove indel)
        public boolean setCigar(Cigar cigar) {
            if ( getCigar().equals(cigar) )
                return false;

            String str = AlignmentUtils.cigarToString(cigar);
            if ( !str.contains("D") && !str.contains("I") )
                return false;

            newCigar = cigar;
            return true;
        }

        // tentatively sets the new start, but it needs to be confirmed later
        public void setAlignmentStart(int start) {
            newStart = start;
        }

        public int getAlignmentStart() {
            return (newStart != -1 ? newStart : read.getAlignmentStart());
        }

        public int getOriginalAlignmentStart() {
            return read.getAlignmentStart();
        }

        // finalizes the changes made.
        // returns true if this record actually changes, false otherwise
        public boolean finalizeUpdate() {
            // if we haven't made any changes, don't do anything
            if ( newCigar == null )
                return false;
            if ( newStart == -1 )
                newStart = read.getAlignmentStart();

            // if it's a paired end read, we need to update the insert size
            if ( read.getReadPairedFlag() ) {
                int insertSize = read.getInferredInsertSize();
                if ( insertSize > 0 ) {
                    read.setCigar(newCigar);
                    read.setInferredInsertSize(insertSize + read.getAlignmentStart() - newStart);
                    read.setAlignmentStart(newStart);
                } else {
                    // note that the correct order of actions is crucial here
                    int oldEnd = read.getAlignmentEnd();
                    read.setCigar(newCigar);
                    read.setAlignmentStart(newStart);
                    read.setInferredInsertSize(insertSize + oldEnd - read.getAlignmentEnd());
                }
            } else {
                read.setCigar(newCigar);
                read.setAlignmentStart(newStart);
            }
            return true;
        }

        public void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }

        public int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }
    }

    private class Consensus {
        public final byte[] str;
        public final ArrayList<Pair<Integer, Integer>> readIndexes;
        public final int positionOnReference;
        public int mismatchSum;
        public Cigar cigar;

        public Consensus(byte[] str, Cigar cigar, int positionOnReference) {
            this.str = str;
            this.cigar = cigar;
            this.positionOnReference = positionOnReference;
            mismatchSum = 0;
            readIndexes = new ArrayList<Pair<Integer, Integer>>();
        }

        public boolean equals(Object o) {
            return ( this == o || (o instanceof Consensus && this.str.equals(((Consensus)o).str)) );
        }

        public boolean equals(Consensus c) {
            return ( this == c || this.str.equals(c.str) );
        }
    }

    private class ReadBin {

        private final ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        private char[] reference = null;
        private GenomeLoc loc = null;

        public ReadBin() { }

        public void add(SAMRecord read, char[] ref) {
            reads.add(read);

            // set up the reference
            if ( reference == null ) {
                reference = ref;
                loc = GenomeLocParser.createGenomeLoc(read);
            } else {
                long lastPosWithRefBase = loc.getStart() + reference.length -1;
                int neededBases = (int)(read.getAlignmentEnd() - lastPosWithRefBase);
                if ( neededBases > ref.length )
                    throw new StingException("Read " + read.getReadName() + " does not overlap the previous read in this interval; please ensure that you are using the same input bam that was used in the RealignerTargetCreator step");
                if ( neededBases > 0 ) {
                    char[] newReference = new char[reference.length + neededBases];
                    System.arraycopy(reference, 0, newReference, 0, reference.length);
                    System.arraycopy(ref, ref.length-neededBases, newReference, reference.length, neededBases);
                    reference = newReference;
                    loc = GenomeLocParser.createGenomeLoc(loc.getContigIndex(), loc.getStart(), loc.getStop()+neededBases);
                }
            }
        }

        public List<SAMRecord> getReads() { return reads; }

        public byte[] getRereference() {
            // upper case it
            for ( int i = 0; i < reference.length; i++ )
                reference[i] = Character.toUpperCase(reference[i]);

            // convert it to a byte array
            byte[] refArray = new byte[reference.length];
            StringUtil.charsToBytes(reference, 0, reference.length, refArray, 0);
            return refArray;
        }

        public GenomeLoc getLocation() { return loc; }

        public int size() { return reads.size(); }

        public void clear() {
            reads.clear();
            reference = null;
            loc = null;
        }

    }
}