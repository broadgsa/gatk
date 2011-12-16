package org.broadinstitute.sting.utils.clipreads;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;


/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 11/27/11
 * Time: 5:17 AM
 * To change this template use File | Settings | File Templates.
 */
public class ClippingOpUnitTest extends BaseTest {

    ClippingOp clippingOp;
    GATKSAMRecord read;

    @BeforeTest
    public void init() {
        read = ClipReadsTestUtils.makeRead();
    }

    @Test
    private void testHardClip() {
//        List<TestParameter> testList = new LinkedList<TestParameter>();
//        testList.add(new TestParameter(0, 0, 1, 4, "1H3M"));//clip 1 base at start
//        testList.add(new TestParameter(3, 3, 0, 3, "3M1H"));//clip 1 base at end
//        testList.add(new TestParameter(0, 1, 2, 4, "2H2M"));//clip 2 bases at start
//        testList.add(new TestParameter(2, 3, 0, 2, "2M2H"));//clip 2 bases at end
//        testList.add(new TestParameter(0, 2, 3, 4, "3H1M"));//clip 3 bases at start
//        testList.add(new TestParameter(1, 3, 0, 1, "1M3H"));//clip 3 bases at end
//
//        for (TestParameter p : testList) {
//            init();
//            clippingOp = new ClippingOp(p.inputStart, p.inputStop);
//            logger.warn("Testing Parameters: " + p.inputStart + "," + p.inputStop + "," + p.substringStart + "," + p.substringStop + "," + p.cigar);
//            ClipReadsTestUtils.testBaseQualCigar(clippingOp.apply(ClippingRepresentation.HARDCLIP_BASES, read),
//                    ClipReadsTestUtils.BASES.substring(p.substringStart, p.substringStop).getBytes(),
//                    ClipReadsTestUtils.QUALS.substring(p.substringStart, p.substringStop).getBytes(),
//                    p.cigar);
//        }

    }

    @Test
    private void testSoftClip() {
//        List<TestParameter> testList = new LinkedList<TestParameter>();
//        testList.add(new TestParameter(0, 0, -1, -1, "1S3M"));//clip 1 base at start
//        testList.add(new TestParameter(3, 3, -1, -1, "3M1S"));//clip 1 base at end
//        testList.add(new TestParameter(0, 1, -1, -1, "2S2M"));//clip 2 bases at start
//        testList.add(new TestParameter(2, 3, -1, -1, "2M2S"));//clip 2 bases at end
//        testList.add(new TestParameter(0, 2, -1, -1, "3S1M"));//clip 3 bases at start
//        testList.add(new TestParameter(1, 3, -1, -1, "1M3S"));//clip 3 bases at end
//
//        for (TestParameter p : testList) {
//            init();
//            clippingOp = new ClippingOp(p.inputStart, p.inputStop);
//            logger.warn("Testing Parameters: " + p.inputStart + "," + p.inputStop + "," + p.cigar);
//            ClipReadsTestUtils.testBaseQualCigar(clippingOp.apply(ClippingRepresentation.SOFTCLIP_BASES, read),
//                    ClipReadsTestUtils.BASES.getBytes(), ClipReadsTestUtils.QUALS.getBytes(), p.cigar);
//        }

    }
}
