package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * testing of the LEGACY version of LocusIteratorByState
 */
public class LegacyLocusIteratorByStateUnitTest extends BaseTest {
    private static SAMFileHeader header;
    private LegacyLocusIteratorByState li;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    private LegacyLocusIteratorByState makeLTBS(List<SAMRecord> reads, ReadProperties readAttributes) {
        return new LegacyLocusIteratorByState(new FakeCloseableIterator<SAMRecord>(reads.iterator()), readAttributes, genomeLocParser, LegacyLocusIteratorByState.sampleListForSAMWithoutReadGroups());
    }

    @Test
    public void testXandEQOperators() {
        final byte[] bases1 = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] bases2 = new byte[] {'A','A','A','C','A','A','A','A','A','C'};

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();

        SAMRecord r1 = ArtificialSAMUtils.createArtificialRead(header,"r1",0,1,10);
        r1.setReadBases(bases1);
        r1.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r1.setCigarString("10M");

        SAMRecord r2 = ArtificialSAMUtils.createArtificialRead(header,"r2",0,1,10);
        r2.setReadBases(bases2);
        r2.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r2.setCigarString("3=1X5=1X");

        SAMRecord r3 = ArtificialSAMUtils.createArtificialRead(header,"r3",0,1,10);
        r3.setReadBases(bases2);
        r3.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r3.setCigarString("3=1X5M1X");

        SAMRecord r4  = ArtificialSAMUtils.createArtificialRead(header,"r4",0,1,10);
        r4.setReadBases(bases2);
        r4.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r4.setCigarString("10M");

        List<SAMRecord> reads = Arrays.asList(r1, r2, r3, r4);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads,readAttributes);

        while (li.hasNext()) {
            AlignmentContext context = li.next();
            ReadBackedPileup pileup = context.getBasePileup();
            Assert.assertEquals(pileup.depthOfCoverage(), 4);
        }
    }

    @Test
    public void testIndelsInRegularPileup() {
        final byte[] bases = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] indelBases = new byte[] {'A','A','A','A','C','T','A','A','A','A','A','A'};

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();

        SAMRecord before = ArtificialSAMUtils.createArtificialRead(header,"before",0,1,10);
        before.setReadBases(bases);
        before.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        before.setCigarString("10M");

        SAMRecord during = ArtificialSAMUtils.createArtificialRead(header,"during",0,2,10);
        during.setReadBases(indelBases);
        during.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        during.setCigarString("4M2I6M");

        SAMRecord after  = ArtificialSAMUtils.createArtificialRead(header,"after",0,3,10);
        after.setReadBases(bases);
        after.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        after.setCigarString("10M");

        List<SAMRecord> reads = Arrays.asList(before, during, after);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads,readAttributes);

        boolean foundIndel = false;
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            ReadBackedPileup pileup = context.getBasePileup().getBaseFilteredPileup(10);
            for (PileupElement p : pileup) {
                if (p.isBeforeInsertion()) {
                    foundIndel = true;
                    Assert.assertEquals(p.getEventLength(), 2, "Wrong event length");
                    Assert.assertEquals(p.getEventBases(), "CT", "Inserted bases are incorrect");
                    break;
               }
            }

         }

         Assert.assertTrue(foundIndel,"Indel in pileup not found");
    }

    @Test
    public void testWholeIndelReadInIsolation() {
        final int firstLocus = 44367789;

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();

        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header, "indelOnly", 0, firstLocus, 76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte)'A',76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte) '@', 76));
        indelOnlyRead.setCigarString("76I");

        List<SAMRecord> reads = Arrays.asList(indelOnlyRead);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, readAttributes);

        // Traditionally, reads that end with indels bleed into the pileup at the following locus.  Verify that the next pileup contains this read
        // and considers it to be an indel-containing read.
        Assert.assertTrue(li.hasNext(),"Should have found a whole-indel read in the normal base pileup without extended events enabled");
        AlignmentContext alignmentContext = li.next();
        Assert.assertEquals(alignmentContext.getLocation().getStart(), firstLocus, "Base pileup is at incorrect location.");
        ReadBackedPileup basePileup = alignmentContext.getBasePileup();
        Assert.assertEquals(basePileup.getReads().size(),1,"Pileup is of incorrect size");
        Assert.assertSame(basePileup.getReads().get(0), indelOnlyRead, "Read in pileup is incorrect");
    }

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) do
     * not negatively influence the ordering of the pileup.
     */
    @Test
    public void testWholeIndelRead() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        SAMRecord leadingRead = ArtificialSAMUtils.createArtificialRead(header,"leading",0,firstLocus,76);
        leadingRead.setReadBases(Utils.dupBytes((byte)'A',76));
        leadingRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        leadingRead.setCigarString("1M75I");

        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header,"indelOnly",0,secondLocus,76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte) 'A', 76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        indelOnlyRead.setCigarString("76I");

        SAMRecord fullMatchAfterIndel = ArtificialSAMUtils.createArtificialRead(header,"fullMatch",0,secondLocus,76);
        fullMatchAfterIndel.setReadBases(Utils.dupBytes((byte)'A',76));
        fullMatchAfterIndel.setBaseQualities(Utils.dupBytes((byte)'@',76));
        fullMatchAfterIndel.setCigarString("75I1M");

        List<SAMRecord> reads = Arrays.asList(leadingRead, indelOnlyRead, fullMatchAfterIndel);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, createTestReadProperties());
        int currentLocus = firstLocus;
        int numAlignmentContextsFound = 0;

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            Assert.assertEquals(alignmentContext.getLocation().getStart(),currentLocus,"Current locus returned by alignment context is incorrect");

            if(currentLocus == firstLocus) {
                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
                Assert.assertEquals(readsAtLocus.size(),1,"Wrong number of reads at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at locus " + currentLocus);
            }
            else if(currentLocus == secondLocus) {
                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
                Assert.assertEquals(readsAtLocus.size(),2,"Wrong number of reads at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(0),indelOnlyRead,"indelOnlyRead absent from pileup at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(1),fullMatchAfterIndel,"fullMatchAfterIndel absent from pileup at locus " + currentLocus);
            }

            currentLocus++;
            numAlignmentContextsFound++;
        }

        Assert.assertEquals(numAlignmentContextsFound, 2, "Found incorrect number of alignment contexts");
    }

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) are represented properly
     */
    @Test
    public void testWholeIndelReadRepresentedTest() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,secondLocus,1);
        read1.setReadBases(Utils.dupBytes((byte) 'A', 1));
        read1.setBaseQualities(Utils.dupBytes((byte) '@', 1));
        read1.setCigarString("1I");

        List<SAMRecord> reads = Arrays.asList(read1);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, createTestReadProperties());

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            PileupElement pe = p.iterator().next();
            Assert.assertTrue(pe.isBeforeInsertion());
            Assert.assertFalse(pe.isAfterInsertion());
            Assert.assertEquals(pe.getEventBases(), "A");
        }

        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,secondLocus,10);
        read2.setReadBases(Utils.dupBytes((byte) 'A', 10));
        read2.setBaseQualities(Utils.dupBytes((byte) '@', 10));
        read2.setCigarString("10I");

        reads = Arrays.asList(read2);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, createTestReadProperties());

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            PileupElement pe = p.iterator().next();
            Assert.assertTrue(pe.isBeforeInsertion());
            Assert.assertFalse(pe.isAfterInsertion());
            Assert.assertEquals(pe.getEventBases(), "AAAAAAAAAA");
        }
    }

    ////////////////////////////////////////////
    // comprehensive LIBS/PileupElement tests //
    ////////////////////////////////////////////

    private static class LIBSTest {


        final String cigar;
        final int readLength;

        private LIBSTest(final String cigar, final int readLength) {
            this.cigar = cigar;
            this.readLength = readLength;
        }
    }

    @DataProvider(name = "LIBSTest")
    public Object[][] createLIBSTestData() {

        //TODO -- when LIBS is fixed this should be replaced to provide all possible permutations of CIGAR strings

        return new Object[][]{
                {new LIBSTest("1I", 1)},
                {new LIBSTest("10I", 10)},
                {new LIBSTest("2M2I2M", 6)},
                {new LIBSTest("2M2I", 4)},
                //TODO -- uncomment these when LIBS is fixed
                //{new LIBSTest("2I2M", 4, Arrays.asList(2,3), Arrays.asList(IS_AFTER_INSERTION_FLAG,0))},
                //{new LIBSTest("1I1M1D1M", 3, Arrays.asList(0,1), Arrays.asList(IS_AFTER_INSERTION_FLAG | IS_BEFORE_DELETION_START_FLAG | IS_BEFORE_DELETED_BASE_FLAG,IS_AFTER_DELETED_BASE_FLAG | IS_AFTER_DELETION_END_FLAG))},
                //{new LIBSTest("1S1I1M", 3, Arrays.asList(2), Arrays.asList(IS_AFTER_INSERTION_FLAG))},
                //{new LIBSTest("1M2D2M", 3)},
                {new LIBSTest("1S1M", 2)},
                {new LIBSTest("1M1S", 2)},
                {new LIBSTest("1S1M1I", 3)}
        };
    }

    @Test(dataProvider = "LIBSTest")
    public void testLIBS(LIBSTest params) {
        final int locus = 44367788;

        SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, locus, params.readLength);
        read.setReadBases(Utils.dupBytes((byte) 'A', params.readLength));
        read.setBaseQualities(Utils.dupBytes((byte) '@', params.readLength));
        read.setCigarString(params.cigar);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(Arrays.asList(read), createTestReadProperties());
        final LIBS_position tester = new LIBS_position(read);

        while ( li.hasNext() ) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            PileupElement pe = p.iterator().next();

            tester.stepForwardOnGenome();

            Assert.assertEquals(pe.isBeforeDeletedBase(), tester.isBeforeDeletedBase);
            Assert.assertEquals(pe.isBeforeDeletionStart(), tester.isBeforeDeletionStart);
            Assert.assertEquals(pe.isAfterDeletedBase(), tester.isAfterDeletedBase);
            Assert.assertEquals(pe.isAfterDeletionEnd(), tester.isAfterDeletionEnd);
            Assert.assertEquals(pe.isBeforeInsertion(), tester.isBeforeInsertion);
            Assert.assertEquals(pe.isAfterInsertion(), tester.isAfterInsertion);
            Assert.assertEquals(pe.isNextToSoftClip(), tester.isNextToSoftClip);
            Assert.assertEquals(pe.getOffset(), tester.getCurrentReadOffset());
        }
    }

    ////////////////////////////////////////////////
    // End comprehensive LIBS/PileupElement tests //
    ////////////////////////////////////////////////

    private static ReadProperties createTestReadProperties() {
        return new ReadProperties(
                Collections.<SAMReaderID>emptyList(),
                new SAMFileHeader(),
                SAMFileHeader.SortOrder.coordinate,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                new ValidationExclusion(),
                Collections.<ReadFilter>emptyList(),
                Collections.<ReadTransformer>emptyList(),
                false,
                (byte) -1
        );
    }
}

class FakeCloseableIterator<T> implements CloseableIterator<T> {
    Iterator<T> iterator;

    public FakeCloseableIterator(Iterator<T> it) {
        iterator = it;
    }

    @Override
    public void close() {}

    @Override
    public boolean hasNext() {
        return iterator.hasNext();
    }

    @Override
    public T next() {
        return iterator.next();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Don't remove!");
    }
}


final class LIBS_position {

    SAMRecord read;

    final int numOperators;
    int currentOperatorIndex = 0;
    int currentPositionOnOperator = 0;
    int currentReadOffset = 0;

    boolean isBeforeDeletionStart = false;
    boolean isBeforeDeletedBase = false;
    boolean isAfterDeletionEnd = false;
    boolean isAfterDeletedBase = false;
    boolean isBeforeInsertion = false;
    boolean isAfterInsertion = false;
    boolean isNextToSoftClip = false;

    boolean sawMop = false;

    public LIBS_position(final SAMRecord read) {
        this.read = read;
        numOperators = read.getCigar().numCigarElements();
    }

    public int getCurrentReadOffset() {
        return Math.max(0, currentReadOffset - 1);
    }

    /**
     * Steps forward on the genome.  Returns false when done reading the read, true otherwise.
     */
    public boolean stepForwardOnGenome() {
        if ( currentOperatorIndex == numOperators )
            return false;

        CigarElement curElement = read.getCigar().getCigarElement(currentOperatorIndex);
        if ( currentPositionOnOperator >= curElement.getLength() ) {
            if ( ++currentOperatorIndex == numOperators )
                return false;

            curElement = read.getCigar().getCigarElement(currentOperatorIndex);
            currentPositionOnOperator = 0;
        }

        switch ( curElement.getOperator() ) {
            case I: // insertion w.r.t. the reference
                if ( !sawMop )
                    break;
            case S: // soft clip
                currentReadOffset += curElement.getLength();
            case H: // hard clip
            case P: // padding
                currentOperatorIndex++;
                return stepForwardOnGenome();

            case D: // deletion w.r.t. the reference
            case N: // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                currentPositionOnOperator++;
                break;

            case M:
            case EQ:
            case X:
                sawMop = true;
                currentReadOffset++;
                currentPositionOnOperator++;
                break;
            default:
                throw new IllegalStateException("No support for cigar op: " + curElement.getOperator());
        }

        final boolean isFirstOp = currentOperatorIndex == 0;
        final boolean isLastOp = currentOperatorIndex == numOperators - 1;
        final boolean isFirstBaseOfOp = currentPositionOnOperator == 1;
        final boolean isLastBaseOfOp = currentPositionOnOperator == curElement.getLength();

        isBeforeDeletionStart = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.D, isLastOp, isLastBaseOfOp);
        isBeforeDeletedBase = isBeforeDeletionStart || (!isLastBaseOfOp && curElement.getOperator() == CigarOperator.D);
        isAfterDeletionEnd = isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.D, isFirstOp, isFirstBaseOfOp);
        isAfterDeletedBase  = isAfterDeletionEnd || (!isFirstBaseOfOp && curElement.getOperator() == CigarOperator.D);
        isBeforeInsertion = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.I, isLastOp, isLastBaseOfOp)
                || (!sawMop && curElement.getOperator() == CigarOperator.I);
        isAfterInsertion = isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.I, isFirstOp, isFirstBaseOfOp);
        isNextToSoftClip = isBeforeOp(read.getCigar(), currentOperatorIndex, CigarOperator.S, isLastOp, isLastBaseOfOp)
                || isAfterOp(read.getCigar(), currentOperatorIndex, CigarOperator.S, isFirstOp, isFirstBaseOfOp);

        return true;
    }

    private static boolean isBeforeOp(final Cigar cigar,
                                      final int currentOperatorIndex,
                                      final CigarOperator op,
                                      final boolean isLastOp,
                                      final boolean isLastBaseOfOp) {
        return  !isLastOp && isLastBaseOfOp && cigar.getCigarElement(currentOperatorIndex+1).getOperator() == op;
    }

    private static boolean isAfterOp(final Cigar cigar,
                                     final int currentOperatorIndex,
                                     final CigarOperator op,
                                     final boolean isFirstOp,
                                     final boolean isFirstBaseOfOp) {
        return  !isFirstOp && isFirstBaseOfOp && cigar.getCigarElement(currentOperatorIndex-1).getOperator() == op;
    }
}
