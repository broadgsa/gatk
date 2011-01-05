/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.genotype.Haplotype;
import org.broadinstitute.sting.utils.pileup.*;

import java.util.*;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

public class HaplotypeScore implements InfoFieldAnnotation, StandardAnnotation {
    private final static boolean DEBUG = false;
    private final static int MIN_CONTEXT_WING_SIZE = 20;
    private final static int MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER = 50;
    private final static char REGEXP_WILDCARD = '.';

    public boolean useRead(PileupElement p) {
	    //return true; // Use all reads
        //return !ReadUtils.is454Read(p.getRead()); // Use all non-454 reads
        return p.getOffset() == -1 || ((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset()); // Use all reads from the filtered context
    }

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( !vc.isBiallelic() || !vc.isSNP() || stratifiedContexts.size() == 0 ) // size 0 means that call was made by someone else and we have no data here
            return null;
        
        final AlignmentContext context = StratifiedAlignmentContext.joinContexts(stratifiedContexts.values());

        final int contextWingSize = Math.min(((int)ref.getWindow().size() - 1)/2, MIN_CONTEXT_WING_SIZE);
        final int contextSize = contextWingSize * 2 + 1;

        final int locus = ref.getLocus().getStart() + (ref.getLocus().getStop() - ref.getLocus().getStart()) / 2;

        // Compute all haplotypes consistent with the current read pileup
        final List<Haplotype> haplotypes = computeHaplotypes(context.getBasePileup(), contextSize, locus);

	    final MathUtils.RunningAverage scoreRA = new MathUtils.RunningAverage();
        if (haplotypes != null) {
            final Set<Map.Entry<String, Genotype>> genotypes = vc.getGenotypes().entrySet();
            for ( final Map.Entry<String, Genotype> genotype : genotypes ) {
                final StratifiedAlignmentContext thisContext = stratifiedContexts.get(genotype.getKey());
                if ( thisContext != null ) {
                    double thisScore = scoreReadsAgainstHaplotypes(haplotypes, thisContext.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup(), contextSize, locus);
                    scoreRA.add(thisScore); // Taking the simple average of all sample's score since the score can be negative and the RMS doesn't make sense
                }
            }
        }

        // annotate the score in the info field
        final Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.4f", scoreRA.mean()));
        return map;
    }

    private class HaplotypeComparator implements Comparator<Haplotype>{

        public int compare(Haplotype a, Haplotype b) {
            if (a.getQualitySum() < b.getQualitySum())
                return 1;
            if (a.getQualitySum() > b.getQualitySum()){
                return -1;
            }
            return 0;
        }
    }

    private List<Haplotype> computeHaplotypes(final ReadBackedPileup pileup, final int contextSize, final int locus) {
        // Compute all possible haplotypes consistent with current pileup
        final PriorityQueue<Haplotype> candidateHaplotypeQueue = new PriorityQueue<Haplotype>(100, new HaplotypeComparator());
        final PriorityQueue<Haplotype> consensusHaplotypeQueue = new PriorityQueue<Haplotype>(MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER, new HaplotypeComparator());

        for ( final PileupElement p : pileup ) {
            if ( useRead(p) ) {
                final Haplotype haplotypeFromRead = getHaplotypeFromRead(p, contextSize, locus);
                candidateHaplotypeQueue.add(haplotypeFromRead);
            }
        }

        // Now that priority queue has been built with all reads at context, we need to merge and find possible segregating haplotypes
        Haplotype elem;
        while ((elem = candidateHaplotypeQueue.poll()) != null)  {
            boolean foundHaplotypeMatch = false;
            Haplotype lastCheckedHaplotype = null;
            for ( final Haplotype haplotypeFromList : consensusHaplotypeQueue ) {
                final Haplotype consensusHaplotype = getConsensusHaplotype(elem, haplotypeFromList);
                if (consensusHaplotype != null)  {
                    foundHaplotypeMatch = true;
                    if (consensusHaplotype.getQualitySum() > haplotypeFromList.getQualitySum()) {
                        consensusHaplotypeQueue.remove(haplotypeFromList);
                        consensusHaplotypeQueue.add(consensusHaplotype);
                    }
                    break;
                }
                else {
                    lastCheckedHaplotype = haplotypeFromList;
                }
            }

            if (!foundHaplotypeMatch && consensusHaplotypeQueue.size() < MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER) {
                consensusHaplotypeQueue.add(elem);
            } else if (!foundHaplotypeMatch && lastCheckedHaplotype != null && elem.getQualitySum() > lastCheckedHaplotype.getQualitySum() ) {
                consensusHaplotypeQueue.remove(lastCheckedHaplotype);
                consensusHaplotypeQueue.add(elem);
            }
        }

        // Now retrieve the two most popular haplotypes
        if (consensusHaplotypeQueue.size() > 0) {
            // Since the consensus haplotypes are in a quality-ordered priority queue, the two best haplotypes are just the first two in the queue
            final Haplotype haplotype1 = consensusHaplotypeQueue.poll();
            Haplotype haplotype2 = consensusHaplotypeQueue.poll();
            if(haplotype2 == null ) { haplotype2 = haplotype1; } // Sometimes only the reference haplotype can be found
            return Arrays.asList(new Haplotype(haplotype1.getBasesAsBytes(), 60), new Haplotype(haplotype2.getBasesAsBytes(), 20)); // These qual values aren't used for anything
        } else {
            return null;
        }
    }

    private Haplotype getHaplotypeFromRead(final PileupElement p, final int contextSize, final int locus) {
        final SAMRecord read = p.getRead();
        int readOffsetFromPileup = p.getOffset();
        final byte[] haplotypeBases = new byte[contextSize];

        for(int i=0; i < contextSize; i++) {
            haplotypeBases[i] = (byte)REGEXP_WILDCARD;
        }

        final double[] baseQualities = new double[contextSize];
        Arrays.fill(baseQualities, 0.0);

        byte[] readBases = read.getReadBases();
        readBases = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getBaseQualities();
        readQuals = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readQuals); // Shift the location of the qual scores based on the Cigar string

        readOffsetFromPileup = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), readOffsetFromPileup, p.getRead().getAlignmentStart(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1)/2;

        for (int i = 0; i < contextSize; i++ ) {
            final int baseOffset = i + baseOffsetStart;
            if ( baseOffset < 0 ) {
                continue;
            }
            if ( baseOffset >= readBases.length ) {
                break;
            }
            if( readQuals[baseOffset] == PileupElement.DELETION_BASE) { readQuals[baseOffset] = PileupElement.DELETION_QUAL; }
            if( readBases[baseOffset] == BaseUtils.N ) { readBases[baseOffset] = (byte)REGEXP_WILDCARD; readQuals[baseOffset] = (byte) 4; } // N's shouldn't be treated as distinct bases
            if( ((double)readQuals[baseOffset]) < 4.0 ) { readQuals[baseOffset] = (byte) 4; } // quals less than 4 are used as codes and don't have actual probabilistic meaning behind them
            haplotypeBases[i] = readBases[baseOffset];
            baseQualities[i] = (double)readQuals[baseOffset];
        }

        return new Haplotype(haplotypeBases, baseQualities);
    }

    private Haplotype getConsensusHaplotype(final Haplotype haplotypeA, final Haplotype haplotypeB) {
        final byte[] a = haplotypeA.getBasesAsBytes();
        final byte[] b = haplotypeB.getBasesAsBytes();

        if (a.length != b.length) {
            throw new ReviewedStingException("Haplotypes a and b must be of same length");
        }

        byte chA, chB;
        final byte wc = (byte)REGEXP_WILDCARD;

        final int length = a.length;
        final byte[] consensusChars = new byte[length];
        final double[] consensusQuals = new double[length];

        final double[] qualsA = haplotypeA.getQuals();
        final double[] qualsB = haplotypeB.getQuals();

        for (int i=0; i < length; i++) {
            chA = a[i];
            chB = b[i];

            if ((chA != chB) && (chA != wc) && (chB != wc))
                return null;

            if ((chA == wc) && (chB == wc)) {
                consensusChars[i] = wc;
                consensusQuals[i] = 0.0;
            }
            else if ((chA == wc)) {
                consensusChars[i] = chB;
                consensusQuals[i] = qualsB[i];
            }
            else if ((chB == wc)){
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i];
            } else {
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i]+qualsB[i];
            }
        }

        return new Haplotype(consensusChars, consensusQuals);
    }

    // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
    private double scoreReadsAgainstHaplotypes(final List<Haplotype> haplotypes, final ReadBackedPileup pileup, final int contextSize, final int locus) {
//        if ( DEBUG ) System.out.printf("HAP1: %s%n", haplotypes.get(0));
//        if ( DEBUG ) System.out.printf("HAP2: %s%n", haplotypes.get(1));

        final ArrayList<double[]> haplotypeScores = new ArrayList<double[]>();
        for ( final PileupElement p : pileup ) {
            // Score all the reads in the pileup, even the filtered ones
            final double[] scores = new double[haplotypes.size()];
            for ( int i = 0; i < haplotypes.size(); i++ ) {
                final Haplotype haplotype = haplotypes.get(i);
                final double score = scoreReadAgainstHaplotype(p, contextSize, haplotype, locus);
                scores[i] = score;
                if ( DEBUG ) { System.out.printf("  vs. haplotype %d = %f%n", i, score); }
            }
            haplotypeScores.add(scores);
        }

        double overallScore = 0.0;
        for ( final double[] readHaplotypeScores : haplotypeScores ) {
            overallScore += MathUtils.arrayMin(readHaplotypeScores);
        }

        return overallScore;
    }

    private double scoreReadAgainstHaplotype(final PileupElement p, final int contextSize, final Haplotype haplotype, final int locus ) {
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
        final byte[] haplotypeBases = haplotype.getBasesAsBytes();
        final SAMRecord read = p.getRead();
        byte[] readBases = read.getReadBases();

        readBases = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getBaseQualities();
        readQuals = AlignmentUtils.readToAlignmentByteArray(p.getRead().getCigar(), readQuals); // Shift the location of the qual scores based on the Cigar string
        int readOffsetFromPileup = p.getOffset();
        readOffsetFromPileup = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), readOffsetFromPileup, p.getRead().getAlignmentStart(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1)/2;

        for ( int i = 0; i < contextSize; i++ ) {
            final int baseOffset = i + baseOffsetStart;
            if ( baseOffset < 0 ) {
                continue;
            }
            if ( baseOffset >= readBases.length ) {
                break;
            }

            final byte haplotypeBase = haplotypeBases[i];
            final byte readBase = readBases[baseOffset];


            final boolean matched = (readBase == haplotypeBase); // Purposefully counting wildcards in the chosen haplotype as a mismatch
            byte qual = readQuals[baseOffset];
            if( qual == PileupElement.DELETION_BASE ) { qual = PileupElement.DELETION_QUAL; }
            if( ((double) qual) < 4.0 ) { qual = (byte) 4; } // quals less than 4 are used as codes and don't have actual probabilistic meaning behind them
            final double e = QualityUtils.qualToErrorProb(qual);
            expected += e;
            mismatches += matched ? e : 1.0 - e / 3.0;

            // a more sophisticated calculation would include the reference quality, but it's nice to actually penalize
            // the mismatching of poorly determined regions of the consensus

//            if ( DEBUG ) {
//                System.out.printf("Read %s: scoring %c vs. %c => e = %f Q%d esum %f vs. msum %f%n",
//                    read.getReadName(), (char)haplotypeBase, (char)readBase, e, readQuals[baseOffset], expected, mismatches);
//            }
        }

        return mismatches - expected;
    }

    public List<String> getKeyNames() { return Arrays.asList("HaplotypeScore"); }
    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("HaplotypeScore", 1, VCFHeaderLineType.Float, "Consistency of the site with two (and only two) segregating haplotypes")); }
}
