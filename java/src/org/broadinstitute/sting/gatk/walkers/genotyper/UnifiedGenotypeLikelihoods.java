package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.ReadBackedPileup;

import net.sf.samtools.SAMRecord;


/**
 A wrapper around the GenotypeLikelihoods class for the UnifiedGenotyper.
 This class incorporates both strand-based likelihoods and a combined likelihoods over both strands.
*/

public class UnifiedGenotypeLikelihoods {

    private GenotypeLikelihoods forwardStrandGL;
    private GenotypeLikelihoods reverseStrandGL;
    private GenotypeLikelihoods combinedGL;

    public UnifiedGenotypeLikelihoods(BaseMismatchModel baseModel, DiploidGenotypePriors priors, EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform defaultPlatform, boolean VERBOSE) {
        forwardStrandGL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
        forwardStrandGL.setVerbose(VERBOSE);
        reverseStrandGL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
        reverseStrandGL.setVerbose(VERBOSE);
        combinedGL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
        combinedGL.setVerbose(VERBOSE);
    }

    public GenotypeLikelihoods getForwardStrandGenotypeLikelihoods() { return forwardStrandGL; }
    public GenotypeLikelihoods getReverseStrandGenotypeLikelihoods() { return reverseStrandGL; }
    public GenotypeLikelihoods getGenotypeLikelihoods() { return combinedGL; }

    public void add(ReadBackedPileup pileup) {
        for (int i = 0; i < pileup.getReads().size(); i++) {
            int offset = pileup.getOffsets().get(i);
            // ignore deletions
            if ( offset == -1 )
                continue;

            SAMRecord read = pileup.getReads().get(i);
            char base = read.getReadString().charAt(offset);
            add(base, read, offset);
        }
    }

    public void add(char base, SAMRecord read, int offset) {
        if ( !read.getReadNegativeStrandFlag() )
            forwardStrandGL.add(base, read.getBaseQualities()[offset], read, offset);
        else
            reverseStrandGL.add(base, read.getBaseQualities()[offset], read, offset);
        combinedGL.add(base, read.getBaseQualities()[offset], read, offset);
    }

    public void setPriors(DiploidGenotypePriors priors) {
        forwardStrandGL.setPriors(priors);
        reverseStrandGL.setPriors(priors);
        combinedGL.setPriors(priors);
    }
}
