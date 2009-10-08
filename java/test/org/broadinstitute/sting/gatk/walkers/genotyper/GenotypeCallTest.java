package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class SSGenotypeCallTest
 *         <p/>
 *         test the SS Genotype call class
 */
public class GenotypeCallTest extends BaseTest {

    // we need a fake GenotypeLikelihoods class
    public class GenotypeLikelihoodsImpl extends GenotypeLikelihoods {
        public boolean cacheIsEnabled() { return false; }

        protected GenotypeLikelihoods getSetCache( char observedBase, byte qualityScore, int ploidy,
                                                 SAMRecord read, int offset, GenotypeLikelihoods val ) {
            return null;
        }

        /**
         * Must be overridden by concrete subclasses
         *
         * @param observedBase
         * @param chromBase
         * @param read
         * @param offset
         *
         * @return
         */
        @Override
        protected double log10PofTrueBaseGivenMiscall(char observedBase, char chromBase, SAMRecord read, int offset) {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public GenotypeLikelihoodsImpl() {
            this.posteriors = new double[10];
            for (int x = 0; x < 10; x++) {
                posteriors[x] = -(x);
            }
            this.likelihoods = new double[10];
            for (int x = 0; x < 10; x++) {
                likelihoods[x] = -(x);
            }
        }
    }


    // make a fake read pile-up
    public Pair<ReadBackedPileup, GenomeLoc> makePileup() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        List<Integer> offset = new ArrayList<Integer>();
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10000);
        for (int x = 0; x < 10; x++) {
            reads.add(ArtificialSAMUtils.createArtificialRead(header, "test", 0, 1, 100));
            offset.add(10 - x);
        }
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
        GenomeLoc loc = GenomeLocParser.createGenomeLoc("chr1", 1, 10);
        return new Pair(new ReadBackedPileup(loc, 'A', reads, offset), loc);
    }


    @Test
    public void testBestVrsRefSame() {
        Pair<ReadBackedPileup, GenomeLoc> myPair = makePileup();
        GenotypeCall call = new GenotypeCall("TESTSAMPLE",myPair.second, 'A', new GenotypeLikelihoodsImpl(), myPair.first);
        Assert.assertEquals(0, call.toVariation().getNegLog10PError(), 0.0000001);
    }

    @Test
    public void testBestVrsRef2() {
        Pair<ReadBackedPileup, GenomeLoc> myPair = makePileup();
        GenotypeCall call2 = new GenotypeCall("TESTSAMPLE",myPair.second, 'T', new GenotypeLikelihoodsImpl(), myPair.first);
        Assert.assertEquals(9, call2.toVariation().getNegLog10PError(), 0.0000001);
    }

    @Test
    public void testBestVrsRef3() {
        Pair<ReadBackedPileup, GenomeLoc> myPair = makePileup();
        GenotypeCall call3 = new GenotypeCall("TESTSAMPLE",myPair.second, 'C', new GenotypeLikelihoodsImpl(), myPair.first);
        Assert.assertEquals(4, call3.toVariation().getNegLog10PError(), 0.0000001);
    }


    @Test
    public void testBestVrsNextSame() {
        Pair<ReadBackedPileup, GenomeLoc> myPair = makePileup();
        GenotypeCall call = new GenotypeCall("TESTSAMPLE",myPair.second, 'A', new GenotypeLikelihoodsImpl(), myPair.first);
        Assert.assertEquals(1, call.getNegLog10PError(), 0.0000001);
    }

    @Test
    public void testBestVrsNext2() {
        Pair<ReadBackedPileup, GenomeLoc> myPair = makePileup();
        GenotypeCall call2 = new GenotypeCall("TESTSAMPLE",myPair.second, 'A', new GenotypeLikelihoodsImpl(), myPair.first);
        Assert.assertEquals(1, call2.getNegLog10PError(), 0.0000001);
    }

    @Test
    public void testBestVrsNext3() {
        Pair<ReadBackedPileup, GenomeLoc> myPair = makePileup();
        GenotypeCall call3 = new GenotypeCall("TESTSAMPLE",myPair.second, 'C', new GenotypeLikelihoodsImpl(), myPair.first);
        Assert.assertEquals(1, call3.getNegLog10PError(), 0.0000001);
    }
}
