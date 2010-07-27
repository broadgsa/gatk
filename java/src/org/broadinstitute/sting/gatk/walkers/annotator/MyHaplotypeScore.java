package org.broadinstitute.sting.gatk.walkers.annotator;

import net.sf.samtools.SAMRecord;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ExtendedPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;

import java.util.*;
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


import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.*;

import java.util.*;
import java.util.regex.Pattern;

import net.sf.samtools.SAMRecord;

// todo -- rename to haplotype penalty
public class MyHaplotypeScore implements InfoFieldAnnotation, StandardAnnotation {
    private final static boolean DEBUG = false;
    private final static int MIN_CONTEXT_WING_SIZE = 10;
    private final static String REGEXP_WILDCARD = ".";

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( !vc.isBiallelic() || !vc.isSNP() || stratifiedContexts.size() == 0 ) // size 0 means that call was made by someone else and we have no data here
            return null;

        AlignmentContext context = StratifiedAlignmentContext.joinContexts(stratifiedContexts.values());

        int contextWingSize = Math.min(((int)ref.getWindow().size() - 1)/2, MIN_CONTEXT_WING_SIZE);
        int contextSize = contextWingSize * 2 + 1;

        // Compute all haplotypes consistent with the current read pileup
        List<HaplotypePair> haplotypePairs = computeHaplotypes(context.getBasePileup(), contextSize);

        // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
        double score = scoreReadsAgainstHaplotypes(haplotypePairs, context.getBasePileup(), contextSize);

        // return the score
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", score));
        return map;
    }

    private List<HaplotypePair> computeHaplotypes(ReadBackedPileup pileup, int contextSize) {
        // Compute all possible haplotypes consistent with current pileup
        ArrayList<HaplotypePair> haplotypeList = new ArrayList<HaplotypePair>();
        ArrayList<Integer> readsPerHaplotype = new ArrayList<Integer>();

        for ( ExtendedPileupElement p : pileup.extendedForeachIterator() ) {
            SAMRecord read = p.getRead();
            int readOffsetFromPileup = p.getOffset();
            int baseOffsetStart = readOffsetFromPileup - (contextSize - 1)/2;

            byte[] haplotypeBases = new byte[contextSize];

            for(int i=0; i < contextSize; i++) {
                haplotypeBases[i] = REGEXP_WILDCARD.getBytes()[0];
            }

            int numUsedLocations = 0;

            for (int i = 0; i < contextSize; i++ ) {
                int baseOffset = i + baseOffsetStart;
                if ( baseOffset < 0 )
                    continue;
                if ( baseOffset >= read.getReadLength() )
                    break;

                haplotypeBases[i] = read.getReadBases()[baseOffset];
                numUsedLocations++;
            }

 //           if (DEBUG)
 //               System.out.println(new String(haplotypeBases));

            boolean foundHaplotypeMatch = false;
            for (int pos = 0; pos < haplotypeList.size(); pos++ ) {
                HaplotypePair elem = haplotypeList.get(pos);

                if (patternMatches(elem.getHaplotype().toString(), new String(haplotypeBases)))  {
                    // this haplotype is consistent with element in list: store whichever has more bases.
                    if (elem.getWeight() < numUsedLocations) {
                        // current read has more bases: remove old and keep new element
                        haplotypeList.set(pos,new HaplotypePair(new Haplotype(haplotypeBases,contextSize), numUsedLocations));
                    }
                    // increment by one the count of reads for current haplotype
                    Integer currReads = readsPerHaplotype.get(pos);
                    currReads++;
                    readsPerHaplotype.set(pos,currReads);
                    foundHaplotypeMatch = true;
                    break;
                }
            }

            if (!foundHaplotypeMatch) {
                haplotypeList.add(new HaplotypePair(new Haplotype(haplotypeBases,contextSize), numUsedLocations));
                readsPerHaplotype.add(1);
            }
        }

        // Now retrieve two most popular haplotypes
        // TODO - quick and dirty solution, could use better data structures to do this automatically
        int bestIdx=0, secondBestIdx=0, bestIdxVal=-1, secondBestIdxVal = -1;

        for (int k=0; k < haplotypeList.size(); k++) {
            int readCount = readsPerHaplotype.get(k);
            if (readCount >= bestIdxVal) {
                secondBestIdx = bestIdx;
                secondBestIdxVal = bestIdxVal;
                bestIdx = k;
                bestIdxVal = readCount;
            }
            else if (readCount >= secondBestIdxVal) {
                // check if current is second best
                secondBestIdx = k;
                secondBestIdxVal = readCount;
            }
        }
        return Arrays.asList(haplotypeList.get(bestIdx),haplotypeList.get(secondBestIdx));
    }

    private boolean patternMatches(String a, String b) {
        // fast poor man's version of Pattern.matches
        if (a.length() != b.length())
            throw new StingException("String a and b must be of same length");

        char chA, chB;
        char wc = REGEXP_WILDCARD.charAt(0);
        
        for (int i=0; i < a.length(); i++) {
            chA = a.charAt(i);
            chB = b.charAt(i);

            if (chA == wc)
                continue;

            if (chB == wc)
                continue;

            if (chA != chB)
                return false;

        }

        return true;
    }
    // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
    private double scoreReadsAgainstHaplotypes(List<HaplotypePair> haplotypePairs, ReadBackedPileup pileup, int contextSize) {
//        if ( DEBUG ) System.out.printf("HAP1: %s%n", haplotypes.get(0));
//        if ( DEBUG ) System.out.printf("HAP1: %s%n", haplotypes.get(1));

        double[][] haplotypeScores = new double[pileup.size()][haplotypePairs.size()];
        for ( ExtendedPileupElement p : pileup.extendedForeachIterator() ) {
            SAMRecord read = p.getRead();
            int readOffsetFromPileup = p.getOffset();

            if ( DEBUG ) System.out.printf("--------------------------------------------- Read %s%n", read.getReadName());
            double m = 10000000;
            for ( int i = 0; i < haplotypePairs.size(); i++ ) {
                Haplotype haplotype = haplotypePairs.get(i).getHaplotype();
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

        public String toString() { return new String(this.bases); }
    }

    private static class HaplotypePair extends Pair<Haplotype,Integer> {
        public HaplotypePair(Haplotype h, Integer i) {
            super(h, i);
        }

        public Haplotype getHaplotype() { return getFirst(); }
        public Integer getWeight() { return getSecond(); }

    }

    public List<String> getKeyNames() { return Arrays.asList("MyHaplotypeScore"); }
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("MyHaplotypeScore", 1, VCFHeaderLineType.Float, "Consistency of the site with two (and only two) segregating haplotypes")); }
}