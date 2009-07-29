package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWindowWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.filters.Platform454Filter;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import net.sf.samtools.*;

import java.util.*;
import java.io.File;
import java.io.FileWriter;
 
@WalkerName("IntervalCleaner")
@ReadFilters({Platform454Filter.class, ZeroMappingQualityReadFilter.class})
public class  IntervalCleanerWalker extends LocusWindowWalker<Integer, Integer> {
    @Argument(fullName="allow454Reads", shortName="454", doc="process 454 reads", required=false)
    boolean allow454 = false;
    @Argument(fullName="OutputCleaned", shortName="O", required=false, doc="Output file (sam or bam) for improved (realigned) reads")
    SAMFileWriter writer = null;
    @Argument(fullName="OutputIndels", shortName="indels", required=false, doc="Output file (text) for the indels found")
    String OUT_INDELS = null;
    @Argument(fullName="OutputCleanedReadsOnly", shortName="cleanedOnly", doc="print out cleaned reads only (otherwise, all reads within the intervals)", required=false)
    boolean cleanedReadsOnly = false;
    @Argument(fullName="statisticsFile", shortName="stats", doc="print out statistics (what does or doesn't get cleaned)", required=false)
    String OUT_STATS = null;
    @Argument(fullName="SNPsFile", shortName="snps", doc="print out whether mismatching columns do or don't get cleaned out", required=false)
    String OUT_SNPS = null;
    @Argument(fullName="LODThresholdForCleaning", shortName="LOD", doc="LOD threshold above which the cleaner will clean", required=false)
    double LOD_THRESHOLD = 5.0;
    @Argument(fullName="EntropyThreshold", shortName="entropy", doc="percentage of mismatches at a locus to be considered having high entropy", required=false)
    double MISMATCH_THRESHOLD = 0.15;
    @Argument(fullName="maxConsensuses", shortName="maxConsensuses", doc="max alternate consensuses to try (necessary to improve performance in deep coverage)", required=false)
    int MAX_CONSENSUSES = 30;
    @Argument(fullName="maxReadsForConsensuses", shortName="greedy", doc="max reads used for finding the alternate consensuses (necessary to improve performance in deep coverage)", required=false)
    int MAX_READS_FOR_CONSENSUSES = 120;

    public static final int MAX_QUAL = 99;

    // fraction of mismatches that need to no longer mismatch for a column to be considered cleaned
    private static final double MISMATCH_COLUMN_CLEANED_FRACTION = 0.75;

    private FileWriter indelOutput = null;
    private FileWriter statsOutput = null;
    private FileWriter snpsOutput = null;
    Random generator;

    // we need to sort the reads ourselves because SAM headers get messed up and claim to be "unsorted" sometimes
    private TreeSet<ComparableSAMRecord> readsToWrite = null;
    private TreeSet<ComparableSAMRecord> nextSetOfReadsToWrite = null;

    public void initialize() {

        if ( LOD_THRESHOLD < 0.0 )
            throw new RuntimeException("LOD threshold cannot be a negative number");
        if ( MISMATCH_THRESHOLD <= 0.0 || MISMATCH_THRESHOLD > 1.0 )
            throw new RuntimeException("Entropy threshold must be a fraction between 0 and 1");

        SAMFileHeader header = getToolkit().getEngine().getSAMHeader();
        if ( writer != null ) {
            readsToWrite = new TreeSet<ComparableSAMRecord>();
        }

        logger.info("Writing into output BAM file");
        logger.info("Temporary space used: "+System.getProperty("java.io.tmpdir"));
        generator = new Random();

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

    public Integer map(RefMetaDataTracker tracker, String ref, LocusContext context) {
        List<SAMRecord> reads = context.getReads();
        ArrayList<SAMRecord> goodReads = new ArrayList<SAMRecord>();
        for ( SAMRecord read : reads ) {
            if ( !read.getReadUnmappedFlag() &&
                 !read.getNotPrimaryAlignmentFlag() &&
                 !read.getDuplicateReadFlag() &&
                 read.getMappingQuality() != 0 &&           
                 read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START &&
		         (allow454 || !Utils.is454Read(read)) )
                goodReads.add(read);
            else if ( writer != null && !cleanedReadsOnly )
                readsToWrite.add(new ComparableSAMRecord(read));
        }

        clean(goodReads, ref, context.getLocation());
        //bruteForceClean(goodReads, ref, context.getLocation().getStart());
        //testCleanWithDeletion();
        //testCleanWithInsertion();

        if ( writer != null ) {
            // Although we can guarantee that reads will be emitted in order WITHIN an interval
            // (since we sort them ourselves), we can't guarantee it BETWEEN intervals.  So,
            // we need to keep track of the PREVIOUS interval's reads: if they don't overlap
            // with those from this interval then we can emit them; otherwise, we merge them.           
            if ( nextSetOfReadsToWrite != null ) {
                if ( readsToWrite.size() > 0 && nextSetOfReadsToWrite.size() > 0 &&
                     readsToWrite.first().getRecord().getAlignmentStart() < nextSetOfReadsToWrite.last().getRecord().getAlignmentStart() ) {
                    nextSetOfReadsToWrite.addAll(readsToWrite);
                } else {
                    Iterator<ComparableSAMRecord> iter = nextSetOfReadsToWrite.iterator();
                    while ( iter.hasNext() )
                        writer.addAlignment(iter.next().getRecord());
                    nextSetOfReadsToWrite = new TreeSet<ComparableSAMRecord>(readsToWrite);
                }
            } else {
                nextSetOfReadsToWrite = new TreeSet<ComparableSAMRecord>(readsToWrite);
            }
            readsToWrite.clear();
        }
        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        if ( nextSetOfReadsToWrite != null ) {
            Iterator<ComparableSAMRecord> iter = nextSetOfReadsToWrite.iterator();
            while ( iter.hasNext() )
                writer.addAlignment(iter.next().getRecord());
        }
        if ( writer != null ) {
            writer.close();
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
        out.println("Saw " + result + " intervals");
    }

 
    private static int mismatchQualitySumIgnoreCigar(AlignedRead aRead, String refSeq, int refIndex) {
        String readSeq = aRead.getReadString();
        String quals = aRead.getBaseQualityString();
        int sum = 0;
        for (int readIndex = 0 ; readIndex < readSeq.length() ; refIndex++, readIndex++ ) {
            if ( refIndex >= refSeq.length() )
                sum += MAX_QUAL;
            else {
                char refChr = refSeq.charAt(refIndex);
                char readChr = readSeq.charAt(readIndex);
                if ( BaseUtils.simpleBaseToBaseIndex(readChr) == -1 ||
                     BaseUtils.simpleBaseToBaseIndex(refChr)  == -1 )
                    continue; // do not count Ns/Xs/etc ?
                if ( Character.toUpperCase(readChr) != Character.toUpperCase(refChr) )
                    sum += (int)quals.charAt(readIndex) - 33;        
            }
        }
        return sum;
    }

    private static boolean readIsClipped(SAMRecord read) {
        final Cigar c = read.getCigar();
        final int n = c.numCigarElements();
        if ( c.getCigarElement(n-1).getOperator() == CigarOperator.S ||
             c.getCigarElement(0).getOperator() == CigarOperator.S) return true;
        return false;
    }

    private void clean(List<SAMRecord> reads, String reference, GenomeLoc interval) {

        long leftmostIndex = interval.getStart();
        ArrayList<SAMRecord> refReads = new ArrayList<SAMRecord>();                   // reads that perfectly match ref
        ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();               // reads that don't perfectly match
        LinkedList<AlignedRead> altAlignmentsToTest = new LinkedList<AlignedRead>();  // should we try to make an alt consensus from the read?
        ArrayList<AlignedRead> leftMovedIndels = new ArrayList<AlignedRead>();
        HashSet<Consensus> altConsenses = new HashSet<Consensus>();                   // list of alt consenses
        int totalMismatchSum = 0;

        // decide which reads potentially need to be cleaned
        for ( SAMRecord read : reads ) {

            // we currently can not deal with clipped reads correctly (or screwy record)
            if ( read.getCigar().numCigarElements() == 0 || readIsClipped(read) ) {
                refReads.add(read);
                continue;
            }

            AlignedRead aRead = new AlignedRead(read);

            // first, move existing indels (for 1 indel reads only) to leftmost position within identical sequence
            int numBlocks = AlignmentUtils.getNumAlignmentBlocks(read);
            if ( numBlocks == 2 ) {
                Cigar newCigar = indelRealignment(read.getCigar(), reference, read.getReadString(), read.getAlignmentStart()-(int)leftmostIndex, 0);
                if ( aRead.setCigar(newCigar) )
                    leftMovedIndels.add(aRead);
            }

            int mismatchScore = mismatchQualitySumIgnoreCigar(aRead, reference, read.getAlignmentStart()-(int)leftmostIndex);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( mismatchScore > 0 ) {
                altReads.add(aRead);
                totalMismatchSum += mismatchScore;
                aRead.setMismatchScoreToReference(mismatchScore);
                // if it has an indel, let's see if that's the best consensus
                if ( numBlocks == 2 )
                    altConsenses.add(createAlternateConsensus(aRead.getAlignmentStart() - (int)leftmostIndex, aRead.getCigar(), reference, aRead.getReadString()));
                else
                    altAlignmentsToTest.add(aRead);
            }
            // otherwise, we can emit it as is
            else {
                refReads.add(read);
            }
        }

        // choose alternate consensuses
        if ( altAlignmentsToTest.size() <= MAX_READS_FOR_CONSENSUSES ) {
            for ( AlignedRead aRead : altAlignmentsToTest ) {
                // do a pairwise alignment against the reference
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getReadString());
                Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, aRead.getReadString());
                if ( c != null)
                    altConsenses.add(c);
            }
        } else {
            // choose alternate consenses randomly
            int readsSeen = 0;
            while ( readsSeen++ < MAX_READS_FOR_CONSENSUSES && altConsenses.size() <= MAX_CONSENSUSES) {
                int index = generator.nextInt(altAlignmentsToTest.size());
                AlignedRead aRead = altAlignmentsToTest.remove(index);
                // do a pairwise alignment against the reference
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getReadString());
                Consensus c = createAlternateConsensus(swConsensus.getAlignmentStart2wrt1(), swConsensus.getCigar(), reference, aRead.getReadString());
                if ( c != null)
                    altConsenses.add(c);
            }
        }

        Consensus bestConsensus = null;
        Iterator<Consensus> iter = altConsenses.iterator();
        while ( iter.hasNext() ) {
            Consensus consensus = iter.next();

            for ( int j = 0; j < altReads.size(); j++ ) {
                AlignedRead toTest = altReads.get(j);
                Pair<Integer, Integer> altAlignment = findBestOffset(consensus.str, toTest);

                // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                int myScore = altAlignment.second;
                if ( myScore >= toTest.getMismatchScoreToReference() )
                    myScore = toTest.getMismatchScoreToReference();
                // keep track of reads that align better to the alternate consensus.
                // By pushing alignments with equal scores to the alternate, it means we'll over-call (het -> hom non ref) but are less likely to under-call (het -> ref, het non ref -> het)
                else
                    consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.first));

                logger.debug(consensus.str +  " vs. " + toTest.getReadString() + " => " + myScore + " - " + altAlignment.first);
                consensus.mismatchSum += myScore;
            }

            logger.debug(consensus.str +  " " + consensus.mismatchSum);
            if ( bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                bestConsensus = consensus;
                logger.debug(consensus.str +  " " + consensus.mismatchSum);
            }
        }

        // if the best alternate consensus has a smaller sum of quality score mismatches (more than
        // the LOD threshold), and it didn't just move around the mismatching columns, then clean!
        double improvement = (bestConsensus == null ? -1 : ((double)(totalMismatchSum - bestConsensus.mismatchSum))/10.0);
        if ( improvement >= LOD_THRESHOLD ) {

            bestConsensus.cigar = indelRealignment(bestConsensus.cigar, reference, bestConsensus.str, bestConsensus.positionOnReference, bestConsensus.positionOnReference);

           // start cleaning the appropriate reads
            for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes ) {
                AlignedRead aRead = altReads.get(indexPair.first);
                updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.second, aRead, (int)leftmostIndex);
            }
            if( !alternateReducesEntropy(altReads, reference, leftmostIndex) ) {
                if ( statsOutput != null ) {
                    try {
                        statsOutput.write(interval.toString());
                        statsOutput.write("\tFAIL (bad indel)\t"); // if improvement > LOD_THRESHOLD *BUT* entropy is not reduced (SNPs still exist)
                        statsOutput.write(Double.toString(improvement));
                        statsOutput.write("\n");
                        statsOutput.flush();
                    } catch (Exception e) {}
                }
            } else {
                logger.debug("CLEAN: " + AlignmentUtils.cigarToString(bestConsensus.cigar) + " " + bestConsensus.str );
                if ( indelOutput != null && bestConsensus.cigar.numCigarElements() > 1 ) {
                    // NOTE: indels are printed out in the format specified for the low-coverage pilot1
                    //  indel calls (tab-delimited): chr position size type sequence
                    StringBuffer str = new StringBuffer();
                    str.append(reads.get(0).getReferenceName());
                    int position = bestConsensus.positionOnReference + bestConsensus.cigar.getCigarElement(0).getLength();
                    str.append("\t" + (leftmostIndex + position - 1));
                    CigarElement ce = bestConsensus.cigar.getCigarElement(1);
                    str.append("\t" + ce.getLength() + "\t" + ce.getOperator() + "\t");
                    if ( ce.getOperator() == CigarOperator.D )
                        str.append(reference.substring(position, position+ce.getLength()));
                    else
                        str.append(bestConsensus.str.substring(position, position+ce.getLength()));
                    str.append("\t" + (((double)(totalMismatchSum - bestConsensus.mismatchSum))/10.0) + "\n");
                    try {
                        indelOutput.write(str.toString());
                        indelOutput.flush();
                    } catch (Exception e) {}
                }
                if ( statsOutput != null ) {
                    try {
                        statsOutput.write(interval.toString());
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
                    AlignedRead aRead = altReads.get(indexPair.first);
                    if ( aRead.finalizeUpdate() ) {
                        aRead.getRead().setMappingQuality(Math.min(aRead.getRead().getMappingQuality() + (int)(improvement/10.0), 255));
                        aRead.getRead().setAttribute("NM", AlignmentUtils.numMismatches(aRead.getRead(), reference, aRead.getRead().getAlignmentStart()-(int)leftmostIndex));
                    }
                }
            }

            // END IF ( improvement >= LOD_THRESHOLD )

        } else if ( statsOutput != null ) { 
            try {
                statsOutput.write(interval.toString());
                statsOutput.write("\tFAIL\t"); // if improvement < LOD_THRESHOLD
                statsOutput.write(Double.toString(improvement));
                statsOutput.write("\n");
                statsOutput.flush();
            } catch (Exception e) {}
        }

        // write them out
        if ( writer != null ) {
            if ( !cleanedReadsOnly ) {
                for ( SAMRecord rec : refReads )
                    readsToWrite.add(new ComparableSAMRecord(rec));
            }
            for ( AlignedRead aRec : leftMovedIndels )
                aRec.finalizeUpdate();
            for ( AlignedRead aRec : altReads ) {
                if ( !cleanedReadsOnly || aRec.wasUpdated() )
                    readsToWrite.add(new ComparableSAMRecord(aRec.getRead()));
            }
        }
    }

    private Consensus createAlternateConsensus(int indexOnRef, Cigar c, String reference, String readStr) {
        if ( indexOnRef < 0 )
            return null;

        // create the new consensus
        StringBuffer sb = new StringBuffer();
        sb.append(reference.substring(0, indexOnRef));
        logger.debug("CIGAR = " + AlignmentUtils.cigarToString(c));

        int indelCount = 0;
        int altIdx = 0;
        int refIdx = indexOnRef;
        boolean ok_flag = true;
        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
            CigarElement ce = c.getCigarElement(i);
            switch( ce.getOperator() ) {
            case D:
                indelCount++;
                refIdx += ce.getLength();
                break;
            case M:
                if ( reference.length() < refIdx+ce.getLength() )
                    ok_flag = false;
                else
                    sb.append(reference.substring(refIdx, refIdx+ce.getLength()));
                refIdx += ce.getLength();
                altIdx += ce.getLength();
                break;
            case I:
                sb.append(readStr.substring(altIdx, altIdx+ce.getLength()));
                altIdx += ce.getLength();
                indelCount++;
                break;
            }
        }
        // make sure that there is at most only a single indel and it aligns appropriately!
        if ( !ok_flag || indelCount != 1 || reference.length() < refIdx )
            return null;

        sb.append(reference.substring(refIdx));
        String altConsensus =  sb.toString(); // alternative consensus sequence we just built from the cuurent read

        return new Consensus(altConsensus, c, indexOnRef);
    }

    private Pair<Integer, Integer> findBestOffset(String ref, AlignedRead read) {
        int attempts = ref.length() - read.getReadLength() + 1;
        int bestScore = mismatchQualitySumIgnoreCigar(read, ref, 0);
        int bestIndex = 0;
        for ( int i = 1; i < attempts; i++ ) {
            // we can't get better than 0!
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
            int score = mismatchQualitySumIgnoreCigar(read, ref, i);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
        }
        return new Pair<Integer, Integer>(bestIndex, bestScore);
    }

    
    private void updateRead(Cigar altCigar, int altPosOnRef, int myPosOnAlt, AlignedRead aRead, int leftmostIndex) {
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

        // the easiest thing to do is to take each case separately
        int endOfFirstBlock = altPosOnRef + altCE1.getLength();
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
        if ( altCE2.getOperator() == CigarOperator.I ) {
            // for reads that end in an insertion
            if ( myPosOnAlt + aRead.getReadLength() < endOfFirstBlock + altCE2.getLength() ) {
                readCigar.add(new CigarElement(myPosOnAlt + aRead.getReadLength() - endOfFirstBlock, CigarOperator.I));
                aRead.setCigar(readCigar);
                return;
            }

            // for reads that start in an insertion
            if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + altCE2.getLength() ) {
                aRead.setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(altCE2.getLength() - (myPosOnAlt - endOfFirstBlock), CigarOperator.I));
                indelOffsetOnRead = myPosOnAlt - endOfFirstBlock;
                sawAlignmentStart = true;
            } else if ( sawAlignmentStart ) {
                readCigar.add(altCE2);
                indelOffsetOnRead = altCE2.getLength();
            }
        } else if ( altCE2.getOperator() == CigarOperator.D ) {
            if ( sawAlignmentStart )
                readCigar.add(altCE2);
            indelOffsetOnRef = altCE2.getLength();
        } else {
            throw new RuntimeException("Operator of middle block is not I or D: " + altCE2.getOperator());
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

    private boolean alternateReducesEntropy(List<AlignedRead> reads, String reference, long leftmostIndex) {
        int[] originalMismatchBases = new int[reference.length()];
        int[] cleanedMismatchBases = new int[reference.length()];
        int[] totalBases = new int[reference.length()];
        for ( int i=0; i < reference.length(); i++ )
            originalMismatchBases[i] = totalBases[i] = 0;

        for (int i=0; i < reads.size(); i++) {
            AlignedRead read = reads.get(i);
            if ( read.getRead().getAlignmentBlocks().size() > 1 )
                 continue;

            int refIdx = read.getOriginalAlignmentStart() - (int)leftmostIndex;
            String readStr = read.getReadString();
            String qualStr = read.getBaseQualityString();

            for (int j=0; j < readStr.length(); j++, refIdx++ ) {
                //                if ( refIdx < 0 || refIdx >= reference.length() ) {
                //                    System.out.println( "Read: "+read.getRead().getReadName() + "; length = " + readStr.length() );
                //                    System.out.println( "Ref left: "+ leftmostIndex +"; ref length=" + reference.length() + "; read alignment start: "+read.getOriginalAlignmentStart() ); 
                //                }
                totalBases[refIdx] += (int)qualStr.charAt(j) - 33;
                if ( Character.toUpperCase(readStr.charAt(j)) != Character.toUpperCase(reference.charAt(refIdx)) )
                    originalMismatchBases[refIdx] += (int)qualStr.charAt(j) - 33;
            }

            // reset and now do the calculation based on the cleaning
            refIdx = read.getAlignmentStart() - (int)leftmostIndex;
            int altIdx = 0;
            Cigar c = read.getCigar();
            for (int j = 0 ; j < c.numCigarElements() ; j++) {
                CigarElement ce = c.getCigarElement(j);
                switch ( ce.getOperator() ) {
                    case M:
                        for (int k = 0 ; k < ce.getLength() ; k++, refIdx++, altIdx++ ) {
                            if ( refIdx < reference.length() && Character.toUpperCase(readStr.charAt(altIdx)) != Character.toUpperCase(reference.charAt(refIdx)) )
                                cleanedMismatchBases[refIdx] += (int)qualStr.charAt(altIdx) - 33;
                        }
                        break;
                    case I:
                        altIdx += ce.getLength();
                        break;
                    case D:
                        refIdx += ce.getLength();
                        break;
                }

            }
        }

        int originalMismatchColumns = 0, cleanedMismatchColumns = 0;
        StringBuffer sb = new StringBuffer();
        for ( int i=0; i < reference.length(); i++ ) {
            if ( cleanedMismatchBases[i] == originalMismatchBases[i] )
                continue;
            boolean didMismatch = false, stillMismatches = false;
            if ( originalMismatchBases[i] > totalBases[i] * MISMATCH_THRESHOLD )  {
                didMismatch = true;
                originalMismatchColumns++;
                if ( cleanedMismatchBases[i] > originalMismatchBases[i] * (1.0 - MISMATCH_COLUMN_CLEANED_FRACTION) ) {
                    stillMismatches = true;
                    cleanedMismatchColumns++;
                }
            } else if ( cleanedMismatchBases[i] > totalBases[i] * MISMATCH_THRESHOLD ) {
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
                
        logger.debug("Original mismatch columns = " + originalMismatchColumns + "; cleaned mismatch columns = " + cleanedMismatchColumns);

        boolean reduces = (originalMismatchColumns == 0 || cleanedMismatchColumns < originalMismatchColumns);
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
    private Cigar indelRealignment(Cigar cigar, String refSeq, String readSeq, int refIndex, int readIndex) {
        if ( cigar.numCigarElements() < 2 ) return cigar; // no indels, nothing to do
        
        CigarElement ce1 = cigar.getCigarElement(0);
        CigarElement ce2 = cigar.getCigarElement(1);

        int difference = 0; // we can move indel 'difference' bases left
        final int indel_length = ce2.getLength();

        String indelString; // inserted or deleted sequence
        int period = 0; // period of the inserted/deleted sequence
        int indelIndexOnRef = refIndex+ce1.getLength() ; // position of the indel on the REF (first deleted base or first base after insertion)
        int indelIndexOnRead = readIndex+ce1.getLength(); // position of the indel on the READ (first insterted base, of first base after deletion)

        if ( ce2.getOperator() == CigarOperator.D )
            indelString = refSeq.substring(indelIndexOnRef, indelIndexOnRef+ce2.getLength()).toUpperCase(); // deleted bases
        else if ( ce2.getOperator() == CigarOperator.I )
            indelString = readSeq.substring(indelIndexOnRead, indelIndexOnRead+ce2.getLength()).toUpperCase(); // get the inserted bases
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
                        char indelChr = indelString.charAt(indelPos);
                        if ( Character.toUpperCase(refSeq.charAt(testRefPos)) != indelChr || BaseUtils.simpleBaseToBaseIndex(indelChr) == -1 ) {
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
            Cigar newCigar = new Cigar();
            newCigar.add(new CigarElement(ce1.getLength()-difference, CigarOperator.M));
            newCigar.add(ce2);
            newCigar.add(new CigarElement(cigar.getCigarElement(2).getLength()+difference, CigarOperator.M));
            //            System.out.println("  FROM:\n"+AlignmentUtils.alignmentToString(cigar,readSeq,refSeq,refIndex));
            //            if ( ce2.getLength() >=2 )
            //                System.out.println("  REALIGNED TO:\n"+AlignmentUtils.alignmentToString(newCigar,readSeq,refSeq,refIndex,readIndex)+"\n");

            logger.debug("Realigning indel: " + AlignmentUtils.cigarToString(cigar) + " to " + AlignmentUtils.cigarToString(newCigar));
            cigar = newCigar;

        }
        return cigar;
    }

    private class AlignedRead {
        private SAMRecord read;
        private Cigar newCigar = null;
        private int newStart = -1;
        private int mismatchScoreToReference;
        private boolean updated = false;

        public AlignedRead(SAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public SAMRecord getRead() {
               return read;
        }

        public String getReadString() {
               return read.getReadString();
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
            updated = true;
            return true;
        }

        public boolean wasUpdated() {
            return updated;
        }

        public String getBaseQualityString() {
            return read.getBaseQualityString();
        }

        public void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }

        public int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }
    }

    private class Consensus {
        public String str;
        public int mismatchSum;
        public int positionOnReference;
        public Cigar cigar;
        public ArrayList<Pair<Integer, Integer>> readIndexes;

        public Consensus(String str, Cigar cigar, int positionOnReference) {
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

    private void testCleanWithInsertion() {
        String reference = "AAAAAACCCCCCAAAAAA";
        // the alternate reference is: "AAAAAACCCTTCCCAAAAAA";
        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        SAMFileHeader header = getToolkit().getEngine().getSAMHeader();
        SAMRecord r1 = new SAMRecord(header);
        r1.setReadName("1");
        r1.setReadString("AACCCCCC");
        r1.setAlignmentStart(4);
        r1.setBaseQualityString("BBBBBBBB");
        SAMRecord r2 = new SAMRecord(header);
        r2.setReadName("2");
        r2.setReadString("AAAACCCT");
        r2.setAlignmentStart(2);
        r2.setBaseQualityString("BBBBBBBB");
        SAMRecord r3 = new SAMRecord(header);
        r3.setReadName("3");
        r3.setReadString("CTTC");
        r3.setAlignmentStart(10);
        r3.setBaseQualityString("BBBB");
        SAMRecord r4 = new SAMRecord(header);
        r4.setReadName("4");
        r4.setReadString("TCCCAA");
        r4.setAlignmentStart(8);
        r4.setBaseQualityString("BBBBBB");
        SAMRecord r5 = new SAMRecord(header);
        r5.setReadName("5");
        r5.setReadString("AAAGAACC");
        r5.setAlignmentStart(0);
        r5.setBaseQualityString("BBBBBBBB");
        SAMRecord r6 = new SAMRecord(header);
        r6.setReadName("6");
        r6.setReadString("CCAAAGAA");
        r6.setAlignmentStart(10);
        r6.setBaseQualityString("BBBBBBBB");
        SAMRecord r7 = new SAMRecord(header);
        r7.setReadName("7");
        r7.setReadString("AACCCTTCCC");
        r7.setAlignmentStart(4);
        r7.setBaseQualityString("BBBBBBBBBB");
        reads.add(r1);
        reads.add(r2);
        reads.add(r3);
        reads.add(r4);
        reads.add(r5);
        reads.add(r6);
        reads.add(r7);
        clean(reads, reference, GenomeLocParser.createGenomeLoc(0,0));
    }

    private void testCleanWithDeletion() {
        String reference = "AAAAAACCCTTCCCAAAAAA";
        // the alternate reference is: "AAAAAACCCCCCAAAAAA";
        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        SAMFileHeader header = getToolkit().getEngine().getSAMHeader();
        SAMRecord r1 = new SAMRecord(header);
        r1.setReadName("1");
        r1.setReadString("ACCCTTCC");
        r1.setAlignmentStart(5);
        r1.setBaseQualityString("BBBBBBBB");
        SAMRecord r2 = new SAMRecord(header);
        r2.setReadName("2");
        r2.setReadString("AAAACCCC");
        r2.setAlignmentStart(2);
        r2.setBaseQualityString("BBBBBBBB");
        SAMRecord r3 = new SAMRecord(header);
        r3.setReadName("3");
        r3.setReadString("CCCC");
        r3.setAlignmentStart(6);
        r3.setBaseQualityString("BBBB");
        SAMRecord r4 = new SAMRecord(header);
        r4.setReadName("4");
        r4.setReadString("CCCCAA");
        r4.setAlignmentStart(10);
        r4.setBaseQualityString("BBBBBB");
        SAMRecord r5 = new SAMRecord(header);
        r5.setReadName("5");
        r5.setReadString("AAAGAACC");
        r5.setAlignmentStart(0);
        r5.setBaseQualityString("BBBBBBBB");
        SAMRecord r6 = new SAMRecord(header);
        r6.setReadName("6");
        r6.setReadString("CCAAAGAA");
        r6.setAlignmentStart(10);
        r6.setBaseQualityString("BBBBBBBB");
        SAMRecord r7 = new SAMRecord(header);
        r7.setReadName("7");
        r7.setReadString("AAAACCCG");
        r7.setAlignmentStart(2);
        r7.setBaseQualityString("BBBBBBBB");
        SAMRecord r8 = new SAMRecord(header);
        r8.setReadName("8");
        r8.setReadString("AACCCCCC");
        r8.setAlignmentStart(4);
        r8.setBaseQualityString("BBBBBBBB");
        reads.add(r1);
        reads.add(r2);
        reads.add(r3);
        reads.add(r4);
        reads.add(r5);
        reads.add(r6);
        reads.add(r7);
        reads.add(r8);
        clean(reads, reference, GenomeLocParser.createGenomeLoc(0,0));
    }
}
