
package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ComparableSAMRecord;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWindowWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.indels.*;

import net.sf.samtools.*;

import java.util.*;
import java.io.File;

@WalkerName("IntervalCleaner")
public class  IntervalCleanerWalker extends LocusWindowWalker<Integer, Integer> {
    @Argument(fullName="maxReadLength", shortName="maxRead", doc="max read length", required=false)
    public int maxReadLength = -1;
    @Argument(fullName="OutputCleaned", shortName="O", required=true, doc="Output file (sam or bam) for improved (realigned) reads")
    public String OUT;
    @Argument(fullName="OutputAllReads", shortName="all", doc="print out all reads (otherwise, just those within the intervals)", required=false)
    public boolean printAllReads = false;

    public static final int MAX_QUAL = 99;

    private SAMFileWriter writer;
    // we need to sort the reads ourselves because SAM headers get messed up and claim to be "unsorted" sometimes
    private TreeSet<ComparableSAMRecord> readsToWrite;

    public void initialize() {
        SAMFileHeader header = getToolkit().getSamReader().getFileHeader();
        writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, new File(OUT));
        readsToWrite = new TreeSet<ComparableSAMRecord>();
    }

    // do we care about reads that are not part of our intervals?
    public boolean actOnNonIntervalReads() {
        return printAllReads;
    }

    // What do we do with the reads that are not part of our intervals?
    public void nonIntervalReadAction(SAMRecord read) {
        writer.addAlignment(read);
    }

    public Integer map(RefMetaDataTracker tracker, String ref, LocusContext context) {
        List<SAMRecord> reads = context.getReads();
        ArrayList<SAMRecord> goodReads = new ArrayList<SAMRecord>();
        for ( SAMRecord read : reads ) {
            if ( read.getReadLength() <= maxReadLength &&
                 !read.getReadUnmappedFlag() &&
                 !read.getNotPrimaryAlignmentFlag() &&
                 !read.getDuplicateReadFlag() &&
                 read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START )
                goodReads.add(read);
            else
                readsToWrite.add(new ComparableSAMRecord(read));
        }
        clean(goodReads, ref, context.getLocation().getStart());
        //bruteForceClean(goodReads, ref, context.getLocation().getStart());
        //testCleanWithDeletion();
        //testCleanWithInsertion();

        Iterator<ComparableSAMRecord> iter = readsToWrite.iterator();
        while ( iter.hasNext() )
            writer.addAlignment(iter.next().getRecord());
        readsToWrite.clear();
        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        out.println("Saw " + result + " intervals");
        writer.close();
    }

    private static int mismatchQualitySum(AlignedRead aRead, String ref, int refIndex) {
        String read = aRead.getReadString();
        String quals = aRead.getBaseQualityString();

        int sum = 0;
        for ( int readIndex = 0 ; readIndex < read.length() ; readIndex++, refIndex++ ) {
            if ( refIndex >= ref.length() )
                sum += MAX_QUAL;
            else if ( Character.toUpperCase(read.charAt(readIndex)) != Character.toUpperCase(ref.charAt(refIndex)) )
                sum += (int)quals.charAt(readIndex) - 33;
        }
        return sum;
    }

    private void clean(List<SAMRecord> reads, String reference, long leftmostIndex) {

        ArrayList<SAMRecord> refReads = new ArrayList<SAMRecord>();
        ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();
        ArrayList<Boolean> altAlignmentsToTest = new ArrayList<Boolean>();
        int totalMismatchSum = 0;

        // decide which reads potentially need to be cleaned
        for ( SAMRecord read : reads ) {
            AlignedRead aRead = new AlignedRead(read);
            int mismatchScore = mismatchQualitySum(aRead, reference, read.getAlignmentStart()-(int)leftmostIndex);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( mismatchScore > 0 ) {
                altReads.add(aRead);
                altAlignmentsToTest.add(true);
                totalMismatchSum += mismatchScore;
                aRead.setMismatchScoreToReference(mismatchScore);
            }
            // otherwise, we can emit it as is
            else {
                refReads.add(read);
            }
        }

        Consensus bestConsensus = null;

        // for each alternative consensus to test, align it to the reference and create an alternative consensus
        for ( int index = 0; index < altAlignmentsToTest.size(); index++ ) {
            if ( altAlignmentsToTest.get(index) ) {

                // do a pairwise alignment against the reference
                AlignedRead aRead = altReads.get(index);
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getReadString());
                int refIdx = swConsensus.getAlignmentStart2wrt1();
                if ( refIdx < 0 )
                    continue;

                // create the new consensus
                StringBuffer sb = new StringBuffer();
                sb.append(reference.substring(0, refIdx));
                Cigar c = swConsensus.getCigar();
                //out.println("CIGAR = " + cigarToString(c));

                int indelCount = 0;
                int altIdx = 0;
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
                            sb.append(aRead.getReadString().substring(altIdx, altIdx+ce.getLength()));
                            altIdx += ce.getLength();
                            indelCount++;
                            break;
                    }
                }
                // make sure that there is at most only a single indel and it aligns appropriately!
                if ( !ok_flag || indelCount > 1 || reference.length() < refIdx )
                    continue;

                sb.append(reference.substring(refIdx));
                String altConsensus =  sb.toString();

                // for each imperfect match to the reference, score it against this alternative
                Consensus consensus = new Consensus(altConsensus, c, swConsensus.getAlignmentStart2wrt1());
                for ( int j = 0; j < altReads.size(); j++ ) {
                    AlignedRead toTest = altReads.get(j);
                    Pair<Integer, Integer> altAlignment = findBestOffset(altConsensus, toTest);

                    // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                    int myScore = altAlignment.getSecond();
                    if ( myScore >= toTest.getMismatchScoreToReference() )
                        myScore = toTest.getMismatchScoreToReference();
                    // keep track of reads that align better to the alternate consensus
                    else
                        consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.getFirst()));

                    logger.info(aRead.getReadString() +  " vs. " + toTest.getReadString() + " => " + myScore + " - " + altAlignment.getFirst());
                    consensus.mismatchSum += myScore;
                    if ( myScore == 0 )
                        // we already know that this is its consensus, so don't bother testing it later
                        altAlignmentsToTest.set(j, false);
                }
                logger.info(aRead.getReadString() +  " " + consensus.mismatchSum);
                if ( bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                    bestConsensus = consensus;
                    logger.info(aRead.getReadString() +  " " + consensus.mismatchSum);
                }
            }
        }

        // if the best alternate consensus has a smaller sum of quality score mismatches, then clean!
        if ( bestConsensus != null && bestConsensus.mismatchSum < totalMismatchSum ) {
            logger.info("CLEAN: " + bestConsensus.str);

            // clean the appropriate reads
            for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes )
                updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.getSecond(), altReads.get(indexPair.getFirst()), (int)leftmostIndex);
        }

        // write them out
        for ( SAMRecord rec : refReads )
            readsToWrite.add(new ComparableSAMRecord(rec));
        for ( AlignedRead aRec : altReads )
            readsToWrite.add(new ComparableSAMRecord(aRec.getRead()));
    }

    private Pair<Integer, Integer> findBestOffset(String ref, AlignedRead read) {
        int attempts = ref.length() - read.getReadLength() + 1;
        int bestScore = mismatchQualitySum(read, ref, 0);
        int bestIndex = 0;
        for ( int i = 1; i < attempts; i++ ) {
            // we can't get better than 0!
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
            int score = mismatchQualitySum(read, ref, i);
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
            aRead.getRead().setAlignmentStart(leftmostIndex + myPosOnAlt);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.getRead().setCigar(readCigar);
            return;
        }

        CigarElement altCE1 = altCigar.getCigarElement(0);
        CigarElement altCE2 = altCigar.getCigarElement(1);

        // the easiest thing to do is to take each case separately
        int endOfFirstBlock = altPosOnRef + altCE1.getLength();
        boolean sawAlignmentStart = false;

        // for reads starting before the indel
        if ( myPosOnAlt < endOfFirstBlock) {
            aRead.getRead().setAlignmentStart(leftmostIndex + myPosOnAlt);
            sawAlignmentStart = true;

            // for reads ending before the indel
            if ( myPosOnAlt + aRead.getReadLength() <= endOfFirstBlock) {
                readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
                aRead.getRead().setCigar(readCigar);
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
                aRead.getRead().setCigar(readCigar);
                return;
            }

            // for reads that start in an insertion
            if ( !sawAlignmentStart && myPosOnAlt < endOfFirstBlock + altCE2.getLength() ) {
                aRead.getRead().setAlignmentStart(leftmostIndex + endOfFirstBlock);
                readCigar.add(new CigarElement(myPosOnAlt - endOfFirstBlock, CigarOperator.I));
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
            aRead.getRead().setAlignmentStart(leftmostIndex + myPosOnAlt + indelOffsetOnRef - indelOffsetOnRead);
            readCigar.add(new CigarElement(aRead.getReadLength(), CigarOperator.M));
            aRead.getRead().setCigar(readCigar);
            return;
        }

        int readRemaining = aRead.getReadLength();
        for ( CigarElement ce : readCigar.getCigarElements() ) {
            if ( ce.getOperator() != CigarOperator.D )
                readRemaining -= ce.getLength();
        }
        if ( readRemaining > 0 )
            readCigar.add(new CigarElement(readRemaining, CigarOperator.M));
        aRead.getRead().setCigar(readCigar);
    }

    private class AlignedRead {
        SAMRecord read;
        int mismatchScoreToReference;

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
            return read.getCigar();
        }

        public void setCigar(Cigar cigar) {
            read.setCigar(cigar);
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
    }

    private void testCleanWithInsertion() {
        String reference = "AAAAAACCCCCCAAAAAA";
        // the alternate reference is: "AAAAAACCCTTCCCAAAAAA";
        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        SAMFileHeader header = getToolkit().getSamReader().getFileHeader();
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
        clean(reads, reference, 0);
    }

    private void testCleanWithDeletion() {
        String reference = "AAAAAACCCTTCCCAAAAAA";
        // the alternate reference is: "AAAAAACCCCCCAAAAAA";
        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        SAMFileHeader header = getToolkit().getSamReader().getFileHeader();
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
        clean(reads, reference, 0);
    }

    private void bruteForceClean(List<SAMRecord> reads, String reference, long leftmostIndex) {

        ArrayList<SAMRecord> refReads = new ArrayList<SAMRecord>();
        ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();
        int totalMismatchSum = 0;

        // decide which reads potentially need to be cleaned
        for ( SAMRecord read : reads ) {
            AlignedRead aRead = new AlignedRead(read);
            int mismatchScore = mismatchQualitySum(aRead, reference, read.getAlignmentStart()-(int)leftmostIndex);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( mismatchScore > 0 ) {
                altReads.add(aRead);
                totalMismatchSum += mismatchScore;
                aRead.setMismatchScoreToReference(mismatchScore);
            }
            // otherwise, we can emit it as is
            else {
                refReads.add(read);
            }
        }

        Consensus bestConsensus = null;

        // for each alternative consensus to test, align it to the reference and create an alternative consensus
        for ( int indelSize = 1; indelSize <= 5; indelSize++ ) {
            for ( int index = 1; index < reference.length(); index++ ) {
                for ( int inOrDel = 0; inOrDel < 2; inOrDel++ ) {

                    // create the new consensus
                    Cigar c = new Cigar();
                    c.add(new CigarElement(index, CigarOperator.M));
                    StringBuffer sb = new StringBuffer();
                    sb.append(reference.substring(0, index));
                    if ( inOrDel == 0 ) {
                        c.add(new CigarElement(indelSize, CigarOperator.D));
                        c.add(new CigarElement(reference.length()-index-indelSize, CigarOperator.M));
                        if ( reference.length() > index+indelSize )
                            sb.append(reference.substring(index+indelSize));
                    } else {
                        c.add(new CigarElement(indelSize, CigarOperator.I));
                        c.add(new CigarElement(reference.length()-index+indelSize, CigarOperator.M));
                        for ( int i = 0; i < indelSize; i++ )
                            sb.append("A");
                        sb.append(reference.substring(index));
                    }
                    String altConsensus =  sb.toString();

                    // for each imperfect match to the reference, score it against this alternative
                    Consensus consensus = new Consensus(altConsensus, c, 0);
                    for ( int j = 0; j < altReads.size(); j++ ) {
                        AlignedRead toTest = altReads.get(j);
                        Pair<Integer, Integer> altAlignment = findBestOffset(altConsensus, toTest);

                        // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                        int myScore = altAlignment.getSecond();
                        if ( myScore >= toTest.getMismatchScoreToReference() )
                            myScore = toTest.getMismatchScoreToReference();
                        // keep track of reads that align better to the alternate consensus
                        else
                            consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.getFirst()));

                        consensus.mismatchSum += myScore;
                    }
                    if ( bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                        bestConsensus = consensus;
                        logger.info(altConsensus +  " " + consensus.mismatchSum);
                    }
                }
            }
        }

        // if the best alternate consensus has a smaller sum of quality score mismatches, then clean!
        if ( bestConsensus != null && bestConsensus.mismatchSum < totalMismatchSum ) {
            logger.info("CLEAN: " + bestConsensus.str);

            // clean the appropriate reads
            for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes )
                updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.getSecond(), altReads.get(indexPair.getFirst()), (int)leftmostIndex);
        }

        // write them out
        for ( SAMRecord rec : refReads )
            readsToWrite.add(new ComparableSAMRecord(rec));
        for ( AlignedRead aRec : altReads )
            readsToWrite.add(new ComparableSAMRecord(aRec.getRead()));
    }

    public static String cigarToString(Cigar cig) {
        StringBuilder b = new StringBuilder();

        for ( int i = 0 ; i < cig.numCigarElements() ; i++ ) {
            char c='?';
            switch ( cig.getCigarElement(i).getOperator() ) {
            case M : c = 'M'; break;
            case D : c = 'D'; break;
            case I : c = 'I'; break;
            }
            b.append(cig.getCigarElement(i).getLength());
            b.append(c);
        }
        return b.toString();
    }
}