/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.indels;

import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.interval.IntervalFileMergingIterator;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.commandline.Argument;

import java.io.File;
import java.io.FileWriter;
import java.util.*;

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

    @Argument(fullName="output", shortName="O", required=false, doc="Output bam")
    protected String writerFilename = null;

    @Argument(fullName="bam_compression", shortName="compress", required=false, doc="Compression level to use for output bams [default:5]")
    protected Integer compressionLevel = 5;

    @Argument(fullName="maxReadsInRam", shortName="maxInRam", doc="max reads allowed to be kept in memory at a time by the SAMFileWriter. "+
                "If too low, the tool may run out of system file descriptors needed to perform sorting; if too high, the tool may run out of memory.", required=false)
    protected int MAX_RECORDS_IN_RAM = 500000;

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

    @Argument(fullName="maxReadsForRealignment", shortName="maxReads", doc="max reads allowed at an interval for realignment; "+
                       "if this value is exceeded, realignment is not attempted and the reads are passed to the output file(s) as-is", required=false)
    protected int MAX_READS = 20000;

    @Argument(fullName="sortInCoordinateOrderEvenThoughItIsHighlyUnsafe", required=false, doc="Should we sort the final bam in coordinate order even though it will be malformed because mate pairs of realigned reads will contain inaccurate information?")
    protected boolean SORT_IN_COORDINATE_ORDER = false;

    @Argument(fullName="realignReadsWithBadMates", required=false, doc="Should we try to realign paired-end reads whose mates map to other chromosomes?")
    protected boolean REALIGN_BADLY_MATED_READS = false;

    @Argument(fullName="no_pg_tag", shortName="noPG", required=false, doc="Don't output the usual PG tag in the realigned bam file header. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.")
    protected boolean NO_PG_TAG = false;

    @Argument(fullName="targetIntervalsAreNotSorted", shortName="targetNotSorted", required=false, doc="This tool assumes that the target interval list is sorted; if the list turns out to be unsorted, it will throw an exception.  Use this argument when your interval list is not sorted to instruct the Realigner to first sort it in memory.")
    protected boolean TARGET_NOT_SORTED = false;

    // the intervals input by the user
    private Iterator<GenomeLoc> intervals = null;

    // the current interval in the list
    private GenomeLoc currentInterval = null;

    // the reads and known indels that fall into the current interval
    private final ReadBin readsToClean = new ReadBin();
    private final ArrayList<SAMRecord> readsNotToClean = new ArrayList<SAMRecord>();
    private final IdentityHashMap<Object, VariantContext> knownIndelsToTry = new IdentityHashMap<Object, VariantContext>();

    // the wrapper around the SAM writer
    private SAMFileWriter writer = null;

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

        if ( !TARGET_NOT_SORTED && IntervalUtils.isIntervalFile(intervalsFile)) {
            // prepare to read intervals one-by-one, as needed (assuming they are sorted). 
            intervals = new IntervalFileMergingIterator( new java.io.File(intervalsFile), IntervalMergingRule.OVERLAPPING_ONLY );
        } else {
            // read in the whole list of intervals for cleaning
            GenomeLocSortedSet locs = IntervalUtils.sortAndMergeIntervals(IntervalUtils.parseIntervalArguments(Arrays.asList(intervalsFile)), IntervalMergingRule.OVERLAPPING_ONLY);
            intervals = locs.iterator();
        }
        currentInterval = intervals.hasNext() ? intervals.next() : null;

        // set up the output writer(s)
        if ( writerFilename != null ) {
            SAMFileWriterFactory factory = new SAMFileWriterFactory();
            factory.setMaxRecordsInRam(MAX_RECORDS_IN_RAM);

            SAMFileHeader header = getToolkit().getSAMFileHeader();
            File file = new File(writerFilename);
            writer = makeWriter(factory, header, file);
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

    private SAMFileWriter makeWriter(SAMFileWriterFactory factory, SAMFileHeader header, File file) {
        if ( SORT_IN_COORDINATE_ORDER )
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            header.setSortOrder(SAMFileHeader.SortOrder.queryname);

        if ( !NO_PG_TAG ) {
            final SAMProgramRecord programRecord = new SAMProgramRecord("GATK IndelRealigner");
            final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
            programRecord.setProgramVersion(headerInfo.getString("org.broadinstitute.sting.gatk.version"));
            header.addProgramRecord( programRecord );
        }

        SAMFileWriter writer = factory.makeBAMWriter(header, false, file, compressionLevel);

        return writer;
    }

    private void emit(final SAMRecord read) {
        if ( writer != null )
            writer.addAlignment(read);
    }

    private void emit(final List<SAMRecord> reads) {
        if ( writer != null ) {
            for ( SAMRecord read : reads )
                writer.addAlignment(read);
        }
    }

    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        if ( currentInterval == null ) {
            emit(read);
            return 0;
        }

        // edge case: when the last target interval abuts the end of the genome, we'll get one of the
        //   unmapped reads while the currentInterval still isn't null.  We need to trigger the cleaning
        //   at this point without trying to create a GenomeLoc.
        if ( read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX ) {
            cleanAndCallMap(ref, read, metaDataTracker, null);
            return 0;
        }

        GenomeLoc readLoc = GenomeLocParser.createGenomeLoc(read);
        // hack to get around unmapped reads having screwy locations
        if ( readLoc.getStop() == 0 )
            readLoc = GenomeLocParser.createGenomeLoc(readLoc.getContigIndex(), readLoc.getStart(), readLoc.getStart());

        if ( readLoc.isBefore(currentInterval) || ReadUtils.is454Read(read) ) {
            // TODO -- it would be nice if we could use indels from 454 reads as alternate consenses
            emit(read);
            return 0;
        }
        else if ( readLoc.overlapsP(currentInterval) ) {
            if ( doNotTryToClean(read) ) {
                readsNotToClean.add(read);
            } else {
                readsToClean.add(read, ref.getBases());
                // add the rods to the list of known variants
                populateKnownIndels(metaDataTracker, null);
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
            cleanAndCallMap(ref, read, metaDataTracker, readLoc);
        }

        return 0;
    }

    private boolean doNotTryToClean(SAMRecord read) {
        return read.getReadUnmappedFlag() ||
                read.getNotPrimaryAlignmentFlag() ||
                read.getMappingQuality() == 0 ||
                read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START ||
                (!REALIGN_BADLY_MATED_READS && read.getReadPairedFlag() && !read.getMateUnmappedFlag() && read.getMateReferenceIndex() != read.getReferenceIndex());
    }

    private void cleanAndCallMap(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker, GenomeLoc readLoc) {
        clean(readsToClean);
        knownIndelsToTry.clear();

        // merge the two sets for emission
        readsNotToClean.addAll(readsToClean.getReads());
        emit(readsNotToClean);
        readsToClean.clear();
        readsNotToClean.clear();

        do {
            currentInterval = intervals.hasNext() ? intervals.next() : null;
        } while ( currentInterval != null && (readLoc == null || currentInterval.isBefore(readLoc)) );

        // call back into map now that the state has been updated
        map(ref, read, metaDataTracker);
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

        if ( writer != null )
            writer.close();

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

    private void populateKnownIndels(ReadMetaDataTracker metaDataTracker, ReferenceContext ref) {
        for ( Collection<GATKFeature> rods : metaDataTracker.getContigOffsetMapping().values() ) {
            Iterator<GATKFeature> rodIter = rods.iterator();
            while ( rodIter.hasNext() ) {
                Object rod = rodIter.next().getUnderlyingObject();
                if ( knownIndelsToTry.containsKey(rod) )
                    continue;
                if ( VariantContextAdaptors.canBeConvertedToVariantContext(rod))
                    knownIndelsToTry.put(rod, VariantContextAdaptors.toVariantContext("", rod, ref));
                else
                    knownIndelsToTry.put(rod, null);
            }
        }
    }

    private static int mismatchQualitySumIgnoreCigar(final AlignedRead aRead, final byte[] refSeq, int refIndex, int quitAboveThisValue) {
        final byte[] readSeq = aRead.getReadBases();
        final byte[] quals = aRead.getBaseQualities();
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
                if ( !BaseUtils.isRegularBase(readChr) || !BaseUtils.isRegularBase(refChr) )
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

    private void clean(ReadBin readsToClean) {

        final List<SAMRecord> reads = readsToClean.getReads();
        if ( reads.size() == 0 )
            return;

        final byte[] reference = readsToClean.getRereference();
        final long leftmostIndex = readsToClean.getLocation().getStart();

        final ArrayList<SAMRecord> refReads = new ArrayList<SAMRecord>();                   // reads that perfectly match ref
        final ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();               // reads that don't perfectly match
        final LinkedList<AlignedRead> altAlignmentsToTest = new LinkedList<AlignedRead>();  // should we try to make an alt consensus from the read?
        final Set<Consensus> altConsenses = new LinkedHashSet<Consensus>();                 // list of alt consenses
        long totalRawMismatchSum = 0;

        // if there are any known indels for this region, get them
        for ( VariantContext knownIndel : knownIndelsToTry.values() ) {
            if ( knownIndel == null || !knownIndel.isIndel() )
                continue;
            byte[] indelStr = knownIndel.isInsertion() ? knownIndel.getAlternateAllele(0).getBases() : Utils.dupBytes((byte)'-', knownIndel.getReference().length());
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

            // we can not deal with screwy records
            if ( read.getCigar().numCigarElements() == 0 ) {
                refReads.add(read);
                continue;
            }

            final AlignedRead aRead = new AlignedRead(read);

            // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
            int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
            if ( numBlocks == 2 ) {
                Cigar newCigar = AlignmentUtils.leftAlignIndel(unclipCigar(read.getCigar()), reference, read.getReadBases(), read.getAlignmentStart()-(int)leftmostIndex, 0);
                aRead.setCigar(newCigar, false);
            }

            final int startOnRef = read.getAlignmentStart()-(int)leftmostIndex;
            final int rawMismatchScore = mismatchQualitySumIgnoreCigar(aRead, reference, startOnRef, Integer.MAX_VALUE);
            //            if ( debugOn ) System.out.println("mismatchScore="+mismatchScore);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( rawMismatchScore > 0 ) {
                altReads.add(aRead);
                //logger.debug("Adding " + aRead.getRead().getReadName() + " with raw mismatch score " + rawMismatchScore + " to non-ref reads");

                if ( !read.getDuplicateReadFlag() )
                    totalRawMismatchSum += rawMismatchScore;
                aRead.setMismatchScoreToReference(rawMismatchScore);
                aRead.setAlignerMismatchScore(AlignmentUtils.mismatchingQualities(aRead.getRead(), reference, startOnRef));

                // if it has an indel, let's see if that's the best consensus
                if ( numBlocks == 2 )  {
                    Consensus c = createAlternateConsensus(startOnRef, aRead.getCigar(), reference, aRead.getReadBases());
                    if ( c != null )
                        altConsenses.add(c);

                }
                else {
                    //                    if ( debugOn ) System.out.println("Going to test...");
                    altAlignmentsToTest.add(aRead);
                }
            }
            // otherwise, we can emit it as is
            else {
                // if ( debugOn ) System.out.println("Emitting as is...");
                //logger.debug("Adding " + aRead.getRead().getReadName() + " with raw mismatch score " + rawMismatchScore + " to ref reads");
                refReads.add(read);
            }
        }

        // choose alternate consensuses
        if ( altAlignmentsToTest.size() <= MAX_READS_FOR_CONSENSUSES ) {
            for ( AlignedRead aRead : altAlignmentsToTest ) {
                // do a pairwise alignment against the reference
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getReadBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
                Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, aRead.getReadBases());
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
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getReadBases(), SW_MATCH, SW_MISMATCH, SW_GAP, SW_GAP_EXTEND);
                Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, aRead.getReadBases());
                if ( c != null )
                    altConsenses.add(c);
            }
        }

        Consensus bestConsensus = null;
        Iterator<Consensus> iter = altConsenses.iterator();

        // if ( debugOn ) System.out.println("------\nChecking consenses...\n--------\n");

        while ( iter.hasNext() ) {
            Consensus consensus = iter.next();
            //logger.debug("Trying new consensus: " + consensus.cigar + " " + new String(consensus.str));

//            if ( DEBUG ) {
//                System.out.println("Checking consensus with alignment at "+consensus.positionOnReference+" cigar "+consensus.cigar);
//                System.out.println(new String(consensus.str));
//                int z = 0;
//                for ( ; z < consensus.positionOnReference; z++ )  System.out.print('.');
//                for ( z=0 ; z < consensus.cigar.getCigarElement(0).getLength() ; z++ ) System.out.print('.');
//                if ( consensus.cigar.getCigarElement(1).getOperator() == CigarOperator.I ) for ( z= 0; z < consensus.cigar.getCigarElement(1).getLength(); z++ )  System.out.print('I');
//                System.out.println();
//            }

            // if ( debugOn ) System.out.println("Consensus: "+consensus.str);

            for ( int j = 0; j < altReads.size(); j++ ) {
                AlignedRead toTest = altReads.get(j);
                Pair<Integer, Integer> altAlignment = findBestOffset(consensus.str, toTest, (int)leftmostIndex);

                // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                int myScore = altAlignment.second;

                if ( myScore > toTest.getAlignerMismatchScore() || myScore >= toTest.getMismatchScoreToReference() )
                    myScore = toTest.getMismatchScoreToReference();
                // keep track of reads that align better to the alternate consensus.
                // By pushing alignments with equal scores to the alternate, it means we'll over-call (het -> hom non ref) but are less likely to under-call (het -> ref, het non ref -> het)
                else
                    consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.first));

                //logger.debug(consensus.cigar +  " vs. " + toTest.getRead().getReadName() + "-" + toTest.getRead().getReadString() + " => " + myScore + " vs. " + toTest.getMismatchScoreToReference());
                if ( !toTest.getRead().getDuplicateReadFlag() )
                    consensus.mismatchSum += myScore;

                // optimization: once the mismatch sum is higher than the best consensus, quit since this one can't win
                //  THIS MUST BE DISABLED IF WE DECIDE TO ALLOW MORE THAN ONE ALTERNATE CONSENSUS!
                if ( bestConsensus != null && consensus.mismatchSum > bestConsensus.mismatchSum )
                    break;
            }

            //logger.debug("Mismatch sum of new consensus: " + consensus.mismatchSum);
            if ( bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                // we do not need this alt consensus, release memory right away!!
                if ( bestConsensus != null )
                    bestConsensus.readIndexes.clear();
                bestConsensus = consensus;
                //logger.debug("New consensus " + bestConsensus.cigar +  " is now best consensus");
            } else {
                // we do not need this alt consensus, release memory right away!!
                consensus.readIndexes.clear();
            }
        }

        // if:
        // 1) the best alternate consensus has a smaller sum of quality score mismatches than the aligned version of the reads,
        // 2) beats the LOD threshold for the sum of quality score mismatches of the raw version of the reads,
        // 3) didn't just move around the mismatching columns (i.e. it actually reduces entropy), 
        // then clean!
        final double improvement = (bestConsensus == null ? -1 : ((double)(totalRawMismatchSum - bestConsensus.mismatchSum))/10.0);
        if ( improvement >= LOD_THRESHOLD ) {

            bestConsensus.cigar = AlignmentUtils.leftAlignIndel(bestConsensus.cigar, reference, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference);

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
                //logger.debug("CLEAN: " + bestConsensus.cigar + " " + bestConsensus.str.toString() + " " + bestConsensus.cigar.numCigarElements() );
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
                    str.append("\t" + (((double)(totalRawMismatchSum - bestConsensus.mismatchSum))/10.0) + "\n");
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
                statsOutput.write(String.format("%s\tFAIL\t%.1f\t%d%n",
                        readsToClean.getLocation().toString(), improvement));
                statsOutput.flush();
            } catch (Exception e) {}
        }
    }

    // create a Consensus from cigar/read strings which originate somewhere on the reference
    private Consensus createAlternateConsensus(final int indexOnRef, final Cigar c, final byte[] reference, final byte[] readStr) {
        if ( indexOnRef < 0 )
            return null;

        // create the new consensus
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(c.numCigarElements()-1);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < indexOnRef; i++)
            sb.append((char)reference[i]);

        int indelCount = 0;
        int altIdx = 0;
        int refIdx = indexOnRef;
        boolean ok_flag = true;
        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
            CigarElement ce = c.getCigarElement(i);
            int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
            case D:
                refIdx += elementLength;
                indelCount++;
                elements.add(ce);
                break;
            case M:
                altIdx += elementLength;
            case N:
                if ( reference.length < refIdx + elementLength )
                    ok_flag = false;
                else  {
                    for (int j = 0; j < elementLength; j++)
                        sb.append((char)reference[refIdx+j]);
                }
                refIdx += elementLength;
                elements.add(new CigarElement(elementLength, CigarOperator.M));
                break;
            case I:
                for (int j = 0; j < elementLength; j++) {
                    if ( ! BaseUtils.isRegularBase(readStr[altIdx+j]) ) {
                        // Insertions with N's in them cause real problems sometimes; it's better to drop them altogether
                        ok_flag=false;
                        break;
                    }
                    sb.append((char)readStr[altIdx + j]);
                }
                altIdx += elementLength;
                indelCount++;
                elements.add(ce);
                break;
            case S:
            default:
                break;
            }
        }
        // make sure that there is at most only a single indel and it aligns appropriately!
        if ( !ok_flag || indelCount != 1 || reference.length < refIdx )
            return null;

        for (int i = refIdx; i < reference.length; i++)
            sb.append((char)reference[i]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, new Cigar(elements), indexOnRef);
    }

    // create a Consensus from just the indel string that falls on the reference
    private Consensus createAlternateConsensus(final int indexOnRef, final byte[] reference, final byte[] indelStr, final boolean isDeletion) {
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
            refIdx += indelStr.length;
            cigar.add(new CigarElement(indelStr.length, CigarOperator.D));
        }
        else {
            for ( byte b : indelStr )
                sb.append((char)b);
            cigar.add(new CigarElement(indelStr.length, CigarOperator.I));
        }

        if ( reference.length - refIdx > 0 )
            cigar.add(new CigarElement(reference.length - refIdx, CigarOperator.M));        
        for (; refIdx < reference.length; refIdx++)
            sb.append((char)reference[refIdx]);
        byte[] altConsensus =  StringUtil.stringToBytes(sb.toString()); // alternative consensus sequence we just built from the current read

        return new Consensus(altConsensus, cigar, 0);
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
            if ( altCE2.getOperator() == CigarOperator.I  || altCE2.getOperator() == CigarOperator.D )
                indelCE=altCE2;
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

        int readRemaining = aRead.getReadBases().length;
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
            final byte[] readStr = read.getReadBases();
            final byte[] quals = read.getBaseQualities();

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
                    case S:
                    default:
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

    private static Cigar unclipCigar(Cigar cigar) {
        ArrayList<CigarElement> elements = new ArrayList<CigarElement>(cigar.numCigarElements());
        for ( CigarElement ce : cigar.getCigarElements() ) {
            if ( !isClipOperator(ce.getOperator()) )
                elements.add(ce);
        }
        return new Cigar(elements);
    }

    private static boolean isClipOperator(CigarOperator op) {
        return op == CigarOperator.S || op == CigarOperator.H || op == CigarOperator.P;
    }

    private class AlignedRead {
        private final SAMRecord read;
        private byte[] readBases = null;
        private byte[] baseQuals = null;
        private Cigar newCigar = null;
        private int newStart = -1;
        private int mismatchScoreToReference = 0;
        private long alignerMismatchScore = 0;

        public AlignedRead(SAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public SAMRecord getRead() {
               return read;
        }

        public int getReadLength() {
            return readBases != null ? readBases.length : read.getReadLength();
        }

        public byte[] getReadBases() {
            if ( readBases == null )
                getUnclippedBases();
            return readBases;
        }

        public byte[] getBaseQualities() {
            if ( baseQuals == null )
                getUnclippedBases();
            return baseQuals;
        }

        // pull out the bases that aren't clipped out
        private void getUnclippedBases() {
            readBases = new byte[getReadLength()];
            baseQuals = new byte[getReadLength()];
            byte[] actualReadBases = read.getReadBases();
            byte[] actualBaseQuals = read.getBaseQualities();
            int fromIndex = 0, toIndex = 0;

            for ( CigarElement ce : read.getCigar().getCigarElements() ) {
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case S:
                        fromIndex += elementLength;
                        break;
                    case M:
                    case I:
                        System.arraycopy(actualReadBases, fromIndex, readBases, toIndex, elementLength);
                        System.arraycopy(actualBaseQuals, fromIndex, baseQuals, toIndex, elementLength);
                        fromIndex += elementLength;
                        toIndex += elementLength;
                    default:
                        break;
                }
            }

            // if we got clipped, trim the array
            if ( fromIndex != toIndex ) {
                byte[] trimmedRB = new byte[toIndex];
                byte[] trimmedBQ = new byte[toIndex];
                System.arraycopy(readBases, 0, trimmedRB, 0, toIndex);
                System.arraycopy(baseQuals, 0, trimmedBQ, 0, toIndex);
                readBases = trimmedRB;
                baseQuals = trimmedBQ;
            }
        }

        public Cigar getCigar() {
            return (newCigar != null ? newCigar : read.getCigar());
        }

        public void setCigar(Cigar cigar) {
            setCigar(cigar, true);
        }

        // tentatively sets the new Cigar, but it needs to be confirmed later
        public void setCigar(Cigar cigar, boolean fixClippedCigar) {
            if ( fixClippedCigar && getReadBases().length < read.getReadLength() )
                cigar = reclipCigar(cigar);

            // no change?
            if ( getCigar().equals(cigar) )
                return;

            // no indel?
            String str = cigar.toString();
            if ( !str.contains("D") && !str.contains("I") )
                return;

            newCigar = cigar;
        }

        // pull out the bases that aren't clipped out
        private Cigar reclipCigar(Cigar cigar) {
            ArrayList<CigarElement> elements = new ArrayList<CigarElement>();

            int i = 0;
            int n = read.getCigar().numCigarElements();
            while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                elements.add(read.getCigar().getCigarElement(i++));

            elements.addAll(cigar.getCigarElements());

            i++;
            while ( i < n && !isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                i++;

            while ( i < n && isClipOperator(read.getCigar().getCigarElement(i).getOperator()) )
                elements.add(read.getCigar().getCigarElement(i++));

            return new Cigar(elements);
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

        public void setAlignerMismatchScore(long score) {
            alignerMismatchScore = score;
        }

        public long getAlignerMismatchScore() {
            return alignerMismatchScore;
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

        @Override
        public boolean equals(Object o) {
            return ( this == o || (o instanceof Consensus && Arrays.equals(this.str,(((Consensus)o).str)) ) );
        }

        public boolean equals(Consensus c) {
            return ( this == c || Arrays.equals(this.str,c.str) ) ;
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(this.str);
        }
    }

    private class ReadBin {

        private final ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        private byte[] reference = null;
        private GenomeLoc loc = null;

        public ReadBin() { }

        public void add(SAMRecord read, byte[] ref) {
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
                    byte[] newReference = new byte[reference.length + neededBases];
                    System.arraycopy(reference, 0, newReference, 0, reference.length);
                    System.arraycopy(ref, ref.length-neededBases, newReference, reference.length, neededBases);
                    reference = newReference;
                    loc = GenomeLocParser.createGenomeLoc(loc.getContigIndex(), loc.getStart(), loc.getStop()+neededBases);
                }
            }
        }

        public List<SAMRecord> getReads() { return reads; }

        public byte[] getRereference() {
            return reference;
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