package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.io.File;
import java.io.IOException;

import net.sf.picard.reference.ReferenceSequence;


/**
 * Created by IntelliJ IDEA.
 * User: mmelgar
 * Date: Jul 28, 2009
 * Time: 3:59:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class VECHomopolymer implements VariantExclusionCriterion {

    private int extent = 5;
    private float frac = 0.8f;
    public IndexedFastaSequenceFile refSeq;
    private String contig;
    public ReferenceSequence contigSeq;

    private boolean exclude = false;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            String[] argPieces = arguments.split(",");
            extent = Integer.valueOf(argPieces[0]);
            frac = Float.valueOf(argPieces[1]);
        }

        File refFile = new File ("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta");

        try {
            refSeq = new IndexedFastaSequenceFile(refFile);
        } catch (IOException e) {
            refSeq = null;
        }
    }

    public boolean useZeroQualityReads() { return false; }

    public String getStudyHeader() {
        return "";
    }

    public String getStudyInfo() {
        return "";
    }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        AlignmentContext alignmentContext = context.getAlignmentContext();
        char[] geno = context.getVariant().getBestGenotype().toCharArray();
        int[] genotype = {-1,-1};
        for (int h = 0; h < 2; h ++) {
            genotype[h] = BaseUtils.simpleBaseToBaseIndex(geno[h]);
        }

        if (contig == null || !alignmentContext.getContig().equals(contig)) {
            contig = alignmentContext.getContig();
            contigSeq = refSeq.getSequence(contig);
        }

        String bases = new String(contigSeq.getBases());
        String prevbases = bases.substring((int) (alignmentContext.getPosition() - extent - 1), (int) (alignmentContext.getPosition() - 1));
        String nextbases = bases.substring((int) (alignmentContext.getPosition() + 1), (int) (alignmentContext.getPosition() + extent + 1));

        int[] prevCounts = {0,0,0,0};
        int[] nextCounts = {0,0,0,0};
        for (int i = 0; i < extent; i ++) {
            int prevBaseIndex = BaseUtils.simpleBaseToBaseIndex(prevbases.toCharArray()[i]);
            int nextBaseIndex = BaseUtils.simpleBaseToBaseIndex(nextbases.toCharArray()[i]);
            if (prevBaseIndex != -1) {
                prevCounts[prevBaseIndex] ++;
            }
            if (nextBaseIndex != -1) {
                nextCounts[nextBaseIndex] ++;
            }
        }

        int prevMaxBase = 0;
        int nextMaxBase = 0;
        int prevMax = prevCounts[0];
        int nextMax = nextCounts[0];
        for (int j = 1; j < 4; j ++) {
            if (prevCounts[j] > prevMax) {
                prevMax = prevCounts[j];
                prevMaxBase = j;
            }
            if (nextCounts[j] > nextMax) {
                nextMax = nextCounts[j];
                nextMaxBase = j;
            }
        }

        int[] homopolymer = {-1,-1};
        if ((float)(prevMax)/(float)(extent) > frac) {homopolymer[0] = prevMaxBase;}
        if ((float)(nextMax)/(float)(extent) > frac) {homopolymer[1] = nextMaxBase;}

        char ref = context.getReferenceContext().getBase();
        boolean checkOne = homopolymer[0] == genotype[0] && homopolymer[0] != BaseUtils.simpleBaseToBaseIndex(ref);
        boolean checkTwo = homopolymer[0] == genotype[1] && homopolymer[0] != BaseUtils.simpleBaseToBaseIndex(ref);
        boolean checkThree = homopolymer[1] == genotype[0] && homopolymer[1] != BaseUtils.simpleBaseToBaseIndex(ref);
        boolean checkFour = homopolymer[1] == genotype[1] && homopolymer[1] != BaseUtils.simpleBaseToBaseIndex(ref);

        exclude =  checkOne || checkTwo || checkThree || checkFour;
    }

    public double inclusionProbability() {
        return exclude ? 0.0 : 1.0;
    }
}