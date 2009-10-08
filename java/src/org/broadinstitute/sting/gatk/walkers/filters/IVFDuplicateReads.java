package org.broadinstitute.sting.gatk.walkers.filters;

import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.HashSet;
import java.util.ArrayList;

import org.broadinstitute.sting.gatk.walkers.genotyper.SingleSampleGenotyper;
import org.broadinstitute.sting.gatk.walkers.genotyper.SSGenotypeCall;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

public class IVFDuplicateReads implements IndependentVariantFeature {
    private SingleSampleGenotyper SSG;
    private double[] likelihoods = new double[10];

    public void initialize(String arguments) {
        SSG = new SingleSampleGenotyper();
        SSG.initialize();
    }

    public void compute(VariantContextWindow contextWindow) {
        List<SAMRecord> reads = contextWindow.getContext().getAlignmentContext(useZeroQualityReads()).getReads();
        List<Integer> offsets = contextWindow.getContext().getAlignmentContext(useZeroQualityReads()).getOffsets();

        HashSet<Integer> alignmentStarts = new HashSet<Integer>();
        ArrayList<SAMRecord> dupeRemovedReads = new ArrayList<SAMRecord>();
        ArrayList<Integer> dupeRemovedOffsets = new ArrayList<Integer>();

        for (int readIndex = 0; readIndex < reads.size(); readIndex++) {
            SAMRecord read = reads.get(readIndex);
            Integer offset = offsets.get(readIndex);
            Integer start = read.getAlignmentStart();

            if (!alignmentStarts.contains(start)) {
                dupeRemovedReads.add(read);
                dupeRemovedOffsets.add(offset);
            }
            
            alignmentStarts.add(start);
        }

        AlignmentContext ac = new AlignmentContext(contextWindow.getContext().getAlignmentContext(useZeroQualityReads()).getLocation(), dupeRemovedReads, dupeRemovedOffsets);
        SSGenotypeCall ssgcall = SSG.map(contextWindow.getContext().getTracker(), contextWindow.getContext().getReferenceContext(), ac);

        double[] newlikelihoods = ssgcall.getLikelihoods();
        for (int likelihoodIndex = 0; likelihoodIndex < newlikelihoods.length; likelihoodIndex++) {
            likelihoods[likelihoodIndex] = newlikelihoods[likelihoodIndex] - contextWindow.getContext().getVariant().getGenotypeLikelihoods()[likelihoodIndex];
        }
    }

    public double[] getLikelihoods() {
        return likelihoods;
    }

    public String getStudyHeader() {
        return "";
    }

    public String getStudyInfo() {
        return "";
    }

    public boolean useZeroQualityReads() { return false; }
}
