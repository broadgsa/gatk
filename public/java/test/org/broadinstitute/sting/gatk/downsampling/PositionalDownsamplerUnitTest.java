package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.*;

// TODO: generalize these tests so that all possible arrangements of 1-4 stacks can be tested
public class PositionalDownsamplerUnitTest extends BaseTest {

    /**
     * -------
     * -------
     *   -------
     *   -------
     *     -------
     *     -------
     */
    @Test
    public void testThreeOverlappingIdenticalStacks() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        PositionalDownsampler<SAMRecord> downsampler = new PositionalDownsampler<SAMRecord>(1000);

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 1, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 25, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 50, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampler.consumeDownsampledItems());

        System.out.println("testThreeOverlappingIdenticalStacks: Downsampled Stack sizes: " + downsampledStackSizes);

        Assert.assertEquals(downsampledStackSizes.size(), 3);
        Assert.assertTrue(downsampledStackSizes.get(0) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(2) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(0) + downsampledStackSizes.get(1) + downsampledStackSizes.get(2) <= 1000);
    }

    /**
     * -------
     * -------
     *         -------
     *         -------
     *                 -------
     *                 -------
     */
    @Test
    public void testThreeNonOverlappingIdenticalStacks() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        PositionalDownsampler<SAMRecord> downsampler = new PositionalDownsampler<SAMRecord>(1000);

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 1, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 201, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 301, 100));
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampler.consumeDownsampledItems());

        System.out.println("testThreeNonOverlappingIdenticalStacks: Downsampled Stack sizes: " + downsampledStackSizes);

        Assert.assertEquals(downsampledStackSizes.size(), 3);
        Assert.assertTrue(downsampledStackSizes.get(0) == 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) == 1000);
        Assert.assertTrue(downsampledStackSizes.get(2) == 1000);
    }

    /**
     * ---
     * ---
     *   -------
     *   -------
     *     -------
     *     -------
     */
    @Test
    public void testThreeStacksWithShortStackAtBeginning() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        PositionalDownsampler<SAMRecord> downsampler = new PositionalDownsampler<SAMRecord>(1000);

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 1, 25));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 20, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 50, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampler.consumeDownsampledItems());

        System.out.println("testThreeStacksWithShortStackAtBeginning: Downsampled Stack sizes: " + downsampledStackSizes);

        Assert.assertEquals(downsampledStackSizes.size(), 3);
        Assert.assertTrue(downsampledStackSizes.get(0) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(2) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(0) + downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) + downsampledStackSizes.get(2) <= 1000);
    }

    /**
     * -------
     * -------
     *   ---
     *   ---
     *      -------
     *      -------
     */
    @Test
    public void testThreeStacksWithShortStackInMiddle() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        PositionalDownsampler<SAMRecord> downsampler = new PositionalDownsampler<SAMRecord>(1000);

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 1, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 25, 25));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 75, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampler.consumeDownsampledItems());

        System.out.println("testThreeStacksWithShortStackInMiddle: Downsampled Stack sizes: " + downsampledStackSizes);

        Assert.assertEquals(downsampledStackSizes.size(), 3);
        Assert.assertTrue(downsampledStackSizes.get(0) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(2) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(0) + downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(0) + downsampledStackSizes.get(2) <= 1000);
    }

    /**
     * ------
     * ------
     *   -------
     *   -------
     *        ---
     *        ---
     */
    @Test
    public void testThreeStacksWithShortStackAtEnd() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        PositionalDownsampler<SAMRecord> downsampler = new PositionalDownsampler<SAMRecord>(1000);

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 1, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 50, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(1500, header, "foo", 0, 135, 25));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampler.consumeDownsampledItems());

        System.out.println("testThreeStacksWithShortStackAtEnd: Downsampled Stack sizes: " + downsampledStackSizes);

        Assert.assertEquals(downsampledStackSizes.size(), 3);
        Assert.assertTrue(downsampledStackSizes.get(0) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(2) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(0) + downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) + downsampledStackSizes.get(2) <= 1000);
    }

    /**
     * -------
     * ----
     *      -------
     *      ----
     *           -------
     *           -------
     */
    @Test
    public void testThreePartiallyOverlappingStacks() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        PositionalDownsampler<SAMRecord> downsampler = new PositionalDownsampler<SAMRecord>(1000);

        downsampler.submit(createStackOfVaryingReads(2000, header, "foo", 0, 1, 100, 50));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfVaryingReads(2000, header, "foo", 0, 75, 100, 50));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(2000, header, "foo", 0, 150, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampler.consumeDownsampledItems());

        System.out.println("testThreePartiallyOverlappingStacks: Downsampled Stack sizes: " + downsampledStackSizes);

        Assert.assertEquals(downsampledStackSizes.size(), 3);
        Assert.assertTrue(downsampledStackSizes.get(0) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(1) <= 1000);
        Assert.assertTrue(downsampledStackSizes.get(2) <= 1000);

        // TODO: need to examine per-base coverage here
    }

    @Test
    public void testNoDownsamplingRequired() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        PositionalDownsampler<SAMRecord> downsampler = new PositionalDownsampler<SAMRecord>(1000);

        downsampler.submit(createStackOfIdenticalReads(300, header, "foo", 0, 1, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(300, header, "foo", 0, 25, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.submit(createStackOfIdenticalReads(300, header, "foo", 0, 50, 100));
        Assert.assertFalse(downsampler.hasDownsampledItems());
        Assert.assertTrue(downsampler.hasPendingItems());

        downsampler.signalEndOfInput();
        Assert.assertTrue(downsampler.hasDownsampledItems());
        Assert.assertFalse(downsampler.hasPendingItems());

        List<Integer> downsampledStackSizes = getDownsampledStackSizesAndVerifySortedness(downsampler.consumeDownsampledItems());

        System.out.println("testNoDownsamplingRequired: Downsampled Stack sizes: " + downsampledStackSizes);

        Assert.assertEquals(downsampledStackSizes.size(), 3);
        Assert.assertTrue(downsampledStackSizes.get(0) == 300);
        Assert.assertTrue(downsampledStackSizes.get(1) == 300);
        Assert.assertTrue(downsampledStackSizes.get(2) == 300);
    }

    @Test
    public void testGATKSAMRecordSupport() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);
        PositionalDownsampler<GATKSAMRecord> downsampler = new PositionalDownsampler<GATKSAMRecord>(1000);

        List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        for ( int i = 0; i < 10; i++ ) {
            reads.add(ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 10, 20 * i + 10));
        }

        downsampler.submit(reads);
        downsampler.signalEndOfInput();
        List<GATKSAMRecord> downsampledReads = downsampler.consumeDownsampledItems();

        Assert.assertTrue(downsampledReads.size() == 10);
    }

    private ArrayList<SAMRecord> createStackOfIdenticalReads( int stackSize, SAMFileHeader header, String name, int refIndex, int alignmentStart, int length ) {
        ArrayList<SAMRecord> stack = new ArrayList<SAMRecord>(stackSize);
        for ( int i = 1; i <= stackSize; i++ ) {
            stack.add(ArtificialSAMUtils.createArtificialRead(header, name, refIndex, alignmentStart, length));
        }
        return stack;
    }

    private ArrayList<SAMRecord> createStackOfVaryingReads( int stackSize, SAMFileHeader header, String name, int refIndex, int alignmentStart, int firstLength, int secondLength ) {
        ArrayList<SAMRecord> stack = createStackOfIdenticalReads(stackSize / 2, header, name, refIndex, alignmentStart, firstLength);
        stack.addAll(createStackOfIdenticalReads(stackSize / 2, header, name, refIndex, alignmentStart, secondLength));
        return stack;
    }

    private List<Integer> getDownsampledStackSizesAndVerifySortedness( List<SAMRecord> downsampledReads ) {
        List<Integer> stackSizes = new ArrayList<Integer>();
        Iterator<SAMRecord> iter = downsampledReads.iterator();
        Assert.assertTrue(iter.hasNext());

        SAMRecord previousRead = iter.next();
        int currentStackSize = 1;

        while ( iter.hasNext() ) {
            SAMRecord currentRead = iter.next();

            if ( ! currentRead.getReferenceIndex().equals(previousRead.getReferenceIndex()) || currentRead.getAlignmentStart() > previousRead.getAlignmentStart() ) {
                stackSizes.add(currentStackSize);
                currentStackSize = 1;
            }
            else if ( currentRead.getAlignmentStart() < previousRead.getAlignmentStart() ) {
                Assert.fail(String.format("Reads are out of order: %s %s", previousRead, currentRead));
            }
            else {
                currentStackSize++;
            }

            previousRead = currentRead;
        }

        stackSizes.add(currentStackSize);
        return stackSizes;
    }
}

