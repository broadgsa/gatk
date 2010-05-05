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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ExtendedPileupElement;

import java.util.*;
import net.sf.samtools.SAMRecord;

// todo -- rename to haplotype penalty
public class HaplotypeScore implements InfoFieldAnnotation, WorkInProgressAnnotation {
    private final static boolean DEBUG = false;
    private final static int MIN_CONTEXT_WING_SIZE = 10;

    // if true, we compute a second haplotype from the reads, instead of constraining ourselves to the reference
    // as a haplotype itself
    private final static boolean USE_NON_REFERENCE_SECOND_HAPLOTYPE = true;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( !vc.isBiallelic() || !vc.isSNP() || stratifiedContexts.size() == 0 ) // size 0 means that call was made by someone else and we have no data here
            return null;

        AlignmentContext context = StratifiedAlignmentContext.joinContexts(stratifiedContexts.values(), StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

        int contextWingSize = Math.min(((int)ref.getWindow().size() - 1)/2, MIN_CONTEXT_WING_SIZE);
        int contextSize = contextWingSize * 2 + 1;

        // calculate
        Haplotype refHaplotype = calcRefHaplotype(vc, ref, context, contextSize);
        Haplotype altHaplotype = new Haplotype(getPileupOfAllele(vc.getAlternateAllele(0), context.getBasePileup()), contextSize);

        //System.exit(1);
        // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
        double score = scoreReadsAgainstHaplotypes(Arrays.asList(refHaplotype, altHaplotype), context.getBasePileup(), contextSize);

        // return the score
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%.2f", score));
        return map;
    }

    // todo -- note that the refPileup.size() won't deal correctly with the situation where the current site is hom-var
    // todo -- but there's nearby het size.  In order to really handle this we need to group reads into two clusters,
    // todo -- but if we are going to do this we might as well just assemble the whole region
    public Haplotype calcRefHaplotype(VariantContext vc, ReferenceContext ref, AlignmentContext context, int contextSize) {
        ReadBackedPileup refPileup = getPileupOfAllele(vc.getReference(), context.getBasePileup());
        if ( USE_NON_REFERENCE_SECOND_HAPLOTYPE && refPileup.size() > 0 ) {
            // we are calculating the reference haplotype from the reads itself -- effectively allows us to
            // have het haplotypes that are hom-var in the surrounding context, indicating that the individual
            // as two alt haplotypes
            return new Haplotype(refPileup, contextSize);
        } else {
            // we are constraining the reference haplotype to really be the reference itself
            int contextWingSize = (contextSize - 1) / 2;
            int refMiddle = (int)(ref.getWindow().size() - 1) / 2;
            int refStart = refMiddle - contextWingSize;
            int refStop = refMiddle + contextWingSize + 1;
            String refString = new String(ref.getBases()).substring(refStart, refStop);
            return new Haplotype(refString.getBytes(), 60);
        }
    }

    // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
    private double scoreReadsAgainstHaplotypes(List<Haplotype> haplotypes, ReadBackedPileup pileup, int contextSize) {
        if ( DEBUG ) System.out.printf("HAP1: %s%n", haplotypes.get(0));
        if ( DEBUG ) System.out.printf("HAP1: %s%n", haplotypes.get(1));

        double[][] haplotypeScores = new double[pileup.size()][haplotypes.size()];
        for ( ExtendedPileupElement p : pileup.extendedForeachIterator() ) {
            SAMRecord read = p.getRead();
            int readOffsetFromPileup = p.getOffset();

            if ( DEBUG ) System.out.printf("--------------------------------------------- Read %s%n", read.getReadName());
            double m = 10000000;
            for ( int i = 0; i < haplotypes.size(); i++ ) {
                Haplotype haplotype = haplotypes.get(i);
                int start = readOffsetFromPileup - (contextSize - 1)/2;
                double score = scoreReadAgainstHaplotype(read, start, contextSize, haplotype);
                haplotypeScores[p.getPileupOffset()][i] = score;
                if ( DEBUG ) System.out.printf("  vs. haplotype %d = %f%n", i, score);
                m = Math.min(score, m);
            }
            if ( DEBUG ) System.out.printf("  => best score was %f%n", m);
        }

        double overallScore = 0.0;
        for ( double[] readHaplotypeScores : haplotypeScores ) {
            overallScore += MathUtils.arrayMin(readHaplotypeScores);
        }

        return overallScore;
    }

    private double scoreReadAgainstHaplotype(SAMRecord read, int baseOffsetStart, int contextSize, Haplotype haplotype ) {
        double expected = 0.0;
        double mismatches = 0.0;

        // What's the expected mismatch rate under the model that this read is actually sampled from
        // this haplotype?  Let's assume the consensus base c is a random choice one of A, C, G, or T, and that
        // the observed base is actually from a c with an error rate e.  Since e is the rate at which we'd
        // see a miscalled c, the expected mismatch rate is really e.  So the expected number of mismatches
        // is just sum_i e_i for i from 1..n for n sites
        //
        // Now, what's the probabilistic sum of mismatches?  Suppose that the base b is equal to c.  Well, it could
        // actually be a miscall in a matching direction, which would happen at a e / 3 rate.  If b != c, then
        // the chance that it is actually a mismatch is 1 - e, since any of the other 3 options would be a mismatch.
        // so the probability-weighted mismatch rate is sum_i ( matched ? e_i / 3 : 1 - e_i ) for i = 1 ... n
        for ( int i = 0; i < contextSize; i++ ) {
            int baseOffset = i + baseOffsetStart;
            if ( baseOffset < 0 )
                continue;
            if ( baseOffset >= read.getReadLength() )
                break;

            byte haplotypeBase = haplotype.bases[i];
            byte readBase = read.getReadBases()[baseOffset];

            boolean matched = BaseUtils.basesAreEqual(readBase, haplotypeBase );
            double e = QualityUtils.qualToErrorProb(read.getBaseQualities()[baseOffset]);
            expected += e;
            mismatches += matched ? e : 1 - e / 3;

            // a more sophisticated calculation would include the reference quality, but it's nice to actually penalize
            // the mismatching of poorly determined regions of the consensus

            if ( DEBUG ) System.out.printf("Read %s: scoring %c vs. %c => e = %f Q%d esum %f vs. msum %f%n",
                    read.getReadName(), (char)haplotypeBase, (char)readBase, e, read.getBaseQualities()[baseOffset], expected, mismatches);
        }

        return mismatches - expected;
    }


    private static final double[] FLAT_BASE_PRIORS = new double[BaseUtils.Base.values().length];
    static {
        for ( int i = 0; i < BaseUtils.Base.values().length; i++ )
            FLAT_BASE_PRIORS[i] = Math.log10(1.0 / BaseUtils.Base.values().length);
    }

    private class Haplotype {
        byte[] bases = null;
        byte[] quals = null;

        /**
         * Create a simple consensus sequence with provided bases and a uniform quality over all bases of qual
         *
         * @param bases
         * @param qual
         */
        Haplotype(byte[] bases, int qual) {
            this.bases = bases;
            quals = new byte[bases.length];
            Arrays.fill(quals, (byte)qual);
        }

        Haplotype(ReadBackedPileup pileup, int contextSize ) {
            this.bases = new byte[contextSize];
            this.quals = new byte[contextSize];
            calculateConsensusOverWindow(pileup, contextSize, (contextSize - 1) / 2);
        }

        private void calculateConsensusOverWindow(ReadBackedPileup pileup, int contextSize, int pileupOffset) {
            // for each context position
            for ( int i = 0; i < contextSize; i++ ) {
                int offsetFromPileup = i - pileupOffset;
                ReadBackedPileup offsetPileup = pileupAtOffset(pileup, offsetFromPileup);
                if ( DEBUG ) System.out.printf("pileup is %s%n", offsetPileup);
                BaseQual bq = calcConsensusAtLocus(offsetPileup, FLAT_BASE_PRIORS);
                this.bases[i] = bq.getBase();
                this.quals[i] = bq.getQual();
                if ( DEBUG ) System.out.printf("  At %d: offset %d bq = %c / %d%n", i, offsetFromPileup, (char)bq.getBase(), bq.getQual());
            }
        }

        private BaseQual calcConsensusAtLocus( ReadBackedPileup pileup, double[] log10priors ) {
            double[] log10BaseLikelihoods = new double[BaseUtils.Base.values().length];

            // loop over a, c, g, t and determine the most likely hypothesis
            for ( BaseUtils.Base base : BaseUtils.Base.values() ) {
                double log10L = log10priors[base.getIndex()];

                for ( PileupElement p : pileup ) {
                    byte qual = p.getQual();
                    if ( qual > 5 ) {
                        double baseP = QualityUtils.qualToProb(qual);
                        double L = base.sameBase(p.getBase()) ? baseP : 1 - baseP;
                        if ( Double.isInfinite(Math.log10(L)) )
                            throw new StingException("BUG -- base likelihood is infinity!");
                        log10L += Math.log10(L);
                    }
                }

                log10BaseLikelihoods[base.getIndex()] = log10L;
            }

            double[] posteriors = MathUtils.normalizeFromLog10(log10BaseLikelihoods, false);
            int mostLikelyIndex = MathUtils.maxElementIndex(posteriors);
            byte mostLikelyBase = BaseUtils.Base.values()[mostLikelyIndex].getBase();                // get the most likely option
            double MAX_CONSENSUS_QUALITY = 0.000001;
            byte qual = QualityUtils.probToQual(posteriors[mostLikelyIndex],MAX_CONSENSUS_QUALITY);  // call posterior calculator here over L over bases
            return new BaseQual(mostLikelyBase, qual);
        }

        public String toString() { return new String(this.bases); }
    }

    private class BaseQual extends Pair<Byte, Byte> {
        public BaseQual(byte base, byte qual) {
            super(base, qual);
        }

        public byte getBase() { return getFirst(); }
        public byte getQual() { return getSecond(); }
    }

    private static ReadBackedPileup pileupAtOffset(ReadBackedPileup pileup, int offsetFromPileup) {
        ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
        ArrayList<Integer> offsets = new ArrayList<Integer>();

        // go through the pileup, read by read, collecting up the context at offset
        for ( PileupElement p : pileup ) {
            // not safe -- doesn't work for indel-containing reads!

            // whole algorithm should be restructured to handle other reads in the window or use LocusIteratorByState
            SAMRecord read = p.getRead();
            int readOffsetInPileup = p.getOffset();
            int neededReadOffset = readOffsetInPileup + offsetFromPileup;
            if ( neededReadOffset >= 0 && neededReadOffset < read.getReadLength() ) {
                reads.add(p.getRead());
                offsets.add(neededReadOffset);
            }
        }

        return new ReadBackedPileup(pileup.getLocation(), reads, offsets);
    }

    private static ReadBackedPileup getPileupOfAllele( Allele allele, ReadBackedPileup pileup ) {
        ArrayList<PileupElement> filteredPileup = new ArrayList<PileupElement>();
        byte alleleBase = allele.getBases()[0]; // assumes SNP 

        for ( PileupElement p : pileup ) {
            if ( p.getBase() == alleleBase ) {
                filteredPileup.add(p);
            }
        }

        return new ReadBackedPileup(pileup.getLocation(), filteredPileup);
    }

    public String getKeyName() { return "HaplotypeScore"; }
    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("HaplotypeScore", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Consistency of the site with two (and only two) segregating haplotypes"); }
}