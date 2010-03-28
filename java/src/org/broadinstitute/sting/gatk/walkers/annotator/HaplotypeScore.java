package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ExtendedPileupElement;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.*;

import net.sf.samtools.SAMRecord;


public class HaplotypeScore implements InfoFieldAnnotation, WorkInProgressAnnotation {
    private final static boolean DEBUG = false;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( !vc.isBiallelic() || !vc.isSNP() || stratifiedContexts.size() == 0 ) // size 0 means that call was made by someone else and we have no data here
            return null;

        AlignmentContext context = StratifiedAlignmentContext.joinContexts(stratifiedContexts.values(), StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

        int contextWingSize = Math.min(((int)ref.getWindow().size() - 1)/2, 10);
        int contextSize = contextWingSize * 2 + 1;

        int refMiddle = (int)(ref.getWindow().size() - 1) / 2;
        int refStart = refMiddle - contextWingSize;
        int refStop = refMiddle + contextWingSize + 1;
        String refString = new String(ref.getBases()).substring(refStart, refStop);
        Consensus refConsensus = new Consensus(refString.getBytes(), 60);

        ReadBackedPileup altPileup = getPileupOfAllele(vc.getAlternateAllele(0), context.getBasePileup());
        Consensus altConsensus = new Consensus(altPileup, contextSize);

        //System.exit(1);
        // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
        double score = scoreReadsAgainstHaplotypes(Arrays.asList(refConsensus, altConsensus), context.getBasePileup(), contextSize);

        // return the score
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%.2f", score));
        return map;
    }

    // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
    private double scoreReadsAgainstHaplotypes(List<Consensus> haplotypes, ReadBackedPileup pileup, int contextSize) {
        if ( DEBUG ) System.out.printf("REF: %s%n", haplotypes.get(0));
        if ( DEBUG ) System.out.printf("ALT: %s%n", haplotypes.get(1));

        double[][] haplotypeScores = new double[pileup.size()][haplotypes.size()];
        for ( ExtendedPileupElement p : pileup.extendedForeachIterator() ) {
            SAMRecord read = p.getRead();
            int readOffsetFromPileup = p.getOffset();

            if ( DEBUG ) System.out.printf("--------------------------------------------- Read %s%n", read.getReadName());
            double m = 10000000;
            for ( int i = 0; i < haplotypes.size(); i++ ) {
                Consensus haplotype = haplotypes.get(i);
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

    private double scoreReadAgainstHaplotype(SAMRecord read, int baseOffsetStart, int contextSize, Consensus haplotype ) {
        double log10LExpected = 0.0;
        double log10LMismatches = 0.0;

        for ( int i = 0; i < contextSize; i++ ) {
            int baseOffset = i + baseOffsetStart;
            if ( baseOffset < 0 )
                continue;
            if ( baseOffset >= read.getReadLength() )
                break;

            byte haplotypeBase = haplotype.bases[i];
            byte readBase = read.getReadBases()[baseOffset];
            if ( ! BaseUtils.basesAreEqual(readBase, haplotypeBase ) ) {
                log10LMismatches += read.getBaseQualities()[baseOffset];
            }

            //System.out.printf("Read %s: scoring %c vs. %c => %f%n", read.getReadName(), (char)haplotypeBase, (char)readBase, log10LMismatches);
        }

        return log10LMismatches; //  - log10LExpected;
    }


    private static final double[] FLAT_BASE_PRIORS = new double[BaseUtils.Base.values().length];
    static {
        for ( int i = 0; i < BaseUtils.Base.values().length; i++ )
            FLAT_BASE_PRIORS[i] = Math.log10(1.0 / BaseUtils.Base.values().length);
    }

    private class Consensus {
        byte[] bases = null;
        byte[] quals = null;

        /**
         * Create a simple consensus sequence with provided bases and a uniform quality over all bases of qual
         *
         * @param bases
         * @param qual
         */
        Consensus(byte[] bases, int qual) {
            this.bases = bases;
            quals = new byte[bases.length];
            Arrays.fill(quals, (byte)qual);
        }

        Consensus(ReadBackedPileup pileup, int contextSize ) {
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
                    double baseP = QualityUtils.qualToProb(p.getQual());
                    double L = base.sameBase(p.getBase()) ? baseP : 1 - baseP;
                    log10L += Math.log10(L);
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