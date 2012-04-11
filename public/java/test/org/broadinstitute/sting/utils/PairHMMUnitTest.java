/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public class PairHMMUnitTest extends BaseTest {
    final static boolean EXTENSIVE_TESTING = true;
    PairHMM hmm = new PairHMM( false ); // reference implementation
    PairHMM bandedHMM = new PairHMM( true ); // algorithm with banding

    // --------------------------------------------------------------------------------
    //
    // Provider
    //
    // --------------------------------------------------------------------------------

    private class BasicLikelihoodTestProvider extends TestDataProvider {
        final String ref, read;
        final byte[] refBasesWithContext, readBasesWithContext;
        final int baseQual, insQual, delQual, gcp;
        final int expectedQual;
        final static String CONTEXT = "ACGTAATGACGATTGCA";
        final static String LEFT_FLANK = "GATTTATCATCGAGTCTGC";
        final static String RIGHT_FLANK = "CATGGATCGTTATCAGCTATCTCGAGGGATTCACTTAACAGTTTTA";

        public BasicLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp) {
            this(ref, read, baseQual, insQual, delQual, expectedQual, gcp, false, false);
        }

        public BasicLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp, final boolean left, final boolean right) {
            super(BasicLikelihoodTestProvider.class, String.format("ref=%s read=%s b/i/d/c quals = %d/%d/%d/%d l/r flank = %b/%b e[qual]=%d", ref, read, baseQual, insQual, delQual, gcp, left, right, expectedQual));
            this.baseQual = baseQual;
            this.delQual = delQual;
            this.insQual = insQual;
            this.gcp = gcp;
            this.read = read;
            this.ref = ref;
            this.expectedQual = expectedQual;

            refBasesWithContext = asBytes(ref, left, right);
            readBasesWithContext = asBytes(read, false, false);
        }

        public double expectedLogL() {
            return expectedQual / -10.0;
        }

        public double tolerance() {
            return 0.1; // TODO FIXME arbitrary
        }

        public double calcLogL() {

            double logL = hmm.computeReadLikelihoodGivenHaplotype(
                    refBasesWithContext, readBasesWithContext,
                    qualAsBytes(baseQual, false), qualAsBytes(insQual, true), qualAsBytes(delQual, true),
                    qualAsBytes(gcp, false));

            return logL;
        }

        private final byte[] asBytes(final String bases, final boolean left, final boolean right) {
            return ( (left ? LEFT_FLANK : "") + CONTEXT + bases + CONTEXT + (right ? RIGHT_FLANK : "")).getBytes();
        }

        private byte[] qualAsBytes(final int phredQual, final boolean doGOP) {
            final byte phredQuals[] = new byte[readBasesWithContext.length];
            // initialize everything to MASSIVE_QUAL so it cannot be moved by HMM
            Arrays.fill(phredQuals, (byte)100);

            // update just the bases corresponding to the provided micro read with the quality scores
            if( doGOP ) {
                phredQuals[0 + CONTEXT.length()] = (byte)phredQual;
            } else {
                for ( int i = 0; i < read.length(); i++)
                    phredQuals[i + CONTEXT.length()] = (byte)phredQual;
            }

            return phredQuals;
        }
    }

    final Random random = new Random(87865573);
    private class BandedLikelihoodTestProvider extends TestDataProvider {
        final String ref, read;
        final byte[] refBasesWithContext, readBasesWithContext;
        final int baseQual, insQual, delQual, gcp;
        final int expectedQual;
        final static String LEFT_CONTEXT = "ACGTAATGACGCTACATGTCGCCAACCGTC";
        final static String RIGHT_CONTEXT = "TACGGCTTCATATAGGGCAATGTGTGTGGCAAAA";
        final static String LEFT_FLANK = "GATTTATCATCGAGTCTGTT";
        final static String RIGHT_FLANK = "CATGGATCGTTATCAGCTATCTCGAGGGATTCACTTAACAGTTTCCGTA";
        final byte[] baseQuals, insQuals, delQuals, gcps;

        public BandedLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp) {
            this(ref, read, baseQual, insQual, delQual, expectedQual, gcp, false, false);
        }

        public BandedLikelihoodTestProvider(final String ref, final String read, final int baseQual, final int insQual, final int delQual, final int expectedQual, final int gcp, final boolean left, final boolean right) {
            super(BandedLikelihoodTestProvider.class, String.format("BANDED: ref=%s read=%s b/i/d/c quals = %d/%d/%d/%d l/r flank = %b/%b e[qual]=%d", ref, read, baseQual, insQual, delQual, gcp, left, right, expectedQual));
            this.baseQual = baseQual;
            this.delQual = delQual;
            this.insQual = insQual;
            this.gcp = gcp;
            this.read = read;
            this.ref = ref;
            this.expectedQual = expectedQual;

            refBasesWithContext = asBytes(ref, left, right);
            readBasesWithContext = asBytes(read, false, false);
            baseQuals = qualAsBytes(baseQual);
            insQuals = qualAsBytes(insQual);
            delQuals = qualAsBytes(delQual);
            gcps = qualAsBytes(gcp, false);
        }

        public double expectedLogL() {
            double logL = hmm.computeReadLikelihoodGivenHaplotype(
                    refBasesWithContext, readBasesWithContext,
                    baseQuals, insQuals, delQuals, gcps);

            return logL;
        }

        public double tolerance() {
            return 0.2; // TODO FIXME arbitrary
        }

        public double calcLogL() {

            double logL = bandedHMM.computeReadLikelihoodGivenHaplotype(
                    refBasesWithContext, readBasesWithContext,
                    baseQuals, insQuals, delQuals, gcps);

            return logL;
        }

        private final byte[] asBytes(final String bases, final boolean left, final boolean right) {
            return ( (left ? LEFT_FLANK : "") + LEFT_CONTEXT + bases + RIGHT_CONTEXT + (right ? RIGHT_FLANK : "")).getBytes();
        }

        private byte[] qualAsBytes(final int phredQual) {
            return qualAsBytes(phredQual, true);
        }

        private byte[] qualAsBytes(final int phredQual, final boolean addRandom) {
            final byte phredQuals[] = new byte[readBasesWithContext.length];
            Arrays.fill(phredQuals, (byte)phredQual);
            if(addRandom) {
                for( int iii = 0; iii < phredQuals.length; iii++) {
                    phredQuals[iii] = (byte) ((int) phredQuals[iii] + (random.nextInt(7) - 3));
                }
            }
            return phredQuals;
        }
    }

    @DataProvider(name = "BasicLikelihoodTestProvider")
    public Object[][] makeBasicLikelihoodTests() {
        // context on either side is ACGTTGCA REF ACGTTGCA
        // test all combinations
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30, 40, 50) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30, 40, 50) : Arrays.asList(40);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(10, 20, 30) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(2,3,4,5,6,7,8,9,10,20) : Arrays.asList(2);

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {

                    // test substitutions
                    for ( final byte refBase : BaseUtils.BASES ) {
                        for ( final byte readBase : BaseUtils.BASES ) {
                            final String ref  = new String(new byte[]{refBase});
                            final String read = new String(new byte[]{readBase});
                            final int expected = refBase == readBase ? 0 : baseQual;
                            new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                        }
                    }

                    // test insertions and deletions
                    for ( final int size : sizes ) {
                        for ( final byte base : BaseUtils.BASES ) {
                            final int expected = indelQual + (size - 2) * gcp;

                            for ( boolean insertionP : Arrays.asList(true, false)) {
                                final String small = Utils.dupString((char)base, 1);
                                final String big = Utils.dupString((char)base, size);

                                final String ref = insertionP ? small : big;
                                final String read = insertionP ? big : small;

                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, true, false);
                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, false, true);
                                new BasicLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, true, true);
                            }
                        }
                    }
                }
            }
        }

        return BasicLikelihoodTestProvider.getTests(BasicLikelihoodTestProvider.class);
    }

    @Test(dataProvider = "BasicLikelihoodTestProvider", enabled = true)
    public void testBasicLikelihoods(BasicLikelihoodTestProvider cfg) {
        double calculatedLogL = cfg.calcLogL();
        double expectedLogL = cfg.expectedLogL();
        logger.warn(String.format("Test: logL calc=%.2f expected=%.2f for %s", calculatedLogL, expectedLogL, cfg.toString()));
        Assert.assertEquals(calculatedLogL, expectedLogL, cfg.tolerance());
    }

    @DataProvider(name = "BandedLikelihoodTestProvider")
    public Object[][] makeBandedLikelihoodTests() {
        // context on either side is ACGTTGCA REF ACGTTGCA
        // test all combinations
        final List<Integer> baseQuals = EXTENSIVE_TESTING ? Arrays.asList(25, 30, 40, 50) : Arrays.asList(30);
        final List<Integer> indelQuals = EXTENSIVE_TESTING ? Arrays.asList(30, 40, 50) : Arrays.asList(40);
        final List<Integer> gcps = EXTENSIVE_TESTING ? Arrays.asList(10, 12) : Arrays.asList(10);
        final List<Integer> sizes = EXTENSIVE_TESTING ? Arrays.asList(2,3,4,5,6,7,8,9,10,20) : Arrays.asList(2);

        for ( final int baseQual : baseQuals ) {
            for ( final int indelQual : indelQuals ) {
                for ( final int gcp : gcps ) {

                    // test substitutions
                    for ( final byte refBase : BaseUtils.BASES ) {
                        for ( final byte readBase : BaseUtils.BASES ) {
                            final String ref  = new String(new byte[]{refBase});
                            final String read = new String(new byte[]{readBase});
                            final int expected = refBase == readBase ? 0 : baseQual;
                            new BandedLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                        }
                    }

                    // test insertions and deletions
                    for ( final int size : sizes ) {
                        for ( final byte base : BaseUtils.BASES ) {
                            final int expected = indelQual + (size - 2) * gcp;

                            for ( boolean insertionP : Arrays.asList(true, false)) {
                                final String small = Utils.dupString((char)base, 1);
                                final String big = Utils.dupString((char)base, size);

                                final String ref = insertionP ? small : big;
                                final String read = insertionP ? big : small;

                                new BandedLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp);
                                new BandedLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, true, false);
                                new BandedLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, false, true);
                                new BandedLikelihoodTestProvider(ref, read, baseQual, indelQual, indelQual, expected, gcp, true, true);
                            }
                        }
                    }
                }
            }
        }

        return BandedLikelihoodTestProvider.getTests(BandedLikelihoodTestProvider.class);
    }

    @Test(dataProvider = "BandedLikelihoodTestProvider", enabled = true)
    public void testBandedLikelihoods(BandedLikelihoodTestProvider cfg) {
        double calculatedLogL = cfg.calcLogL();
        double expectedLogL = cfg.expectedLogL();
        logger.warn(String.format("Test: logL calc=%.2f expected=%.2f for %s", calculatedLogL, expectedLogL, cfg.toString()));
        Assert.assertEquals(calculatedLogL, expectedLogL, cfg.tolerance());
    }
}