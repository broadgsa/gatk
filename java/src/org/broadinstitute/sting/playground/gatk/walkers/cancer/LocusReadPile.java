package org.broadinstitute.sting.playground.gatk.walkers.cancer;

import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;

/**
 * Created by IntelliJ IDEA.
 * User: kcibul
 * Date: Jun 27, 2009
 * Time: 3:41:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class LocusReadPile {
    public List<SAMRecord> reads = new ArrayList<SAMRecord>();
    public List<Integer> offsets = new ArrayList<Integer>();
    public char refBase;
    public int minQualityScore;
    public GenotypeLikelihoods likelihoods = new GenotypeLikelihoods();
    public QualitySums qualitySums = new QualitySums();

    public LocusReadPile(char refBase) {
        this(refBase, 0);
    }

    public LocusReadPile(char refBase, int minQualityScore) {
        this.refBase = refBase;
        this.minQualityScore = minQualityScore;
    }

    public void add(SAMRecord read, int offset) {
        add(read, offset, false);
    }

    public void add(SAMRecord read, int offset, boolean allowMapq0ForQualSum) {
        char base = read.getReadString().charAt(offset);
        byte qual = read.getBaseQualities()[offset];

        if (base == 'N' || base == 'n') { return; }

        if (read.getMappingQuality() > 0) {
            reads.add(read);
            offsets.add(offset);

            likelihoods.add(refBase, base, qual);
        }

        if (read.getMappingQuality() == 0 && !allowMapq0ForQualSum) { return; }

        if (qual > this.minQualityScore) qualitySums.incrementSum(base, qual);
    }

    public String getLocusBases() {
        return getLocusBases(0);
    }

    public String getLocusBases(int locusOffset) {
        StringBuilder sb = new StringBuilder();
        for(int i=0; i<reads.size(); i++) {
            SAMRecord read = reads.get(i);
            int readOffset = offsets.get(i);

            int offset = readOffset + locusOffset;
            if (offset >= 0 && offset < read.getReadString().length()) {
                char base = read.getReadString().charAt(offset);
                sb.append(base);
            }
        }
        return sb.toString();
    }

    public double getAltVsRef(char altAllele) {
        double[] tumorRefAlt = extractRefAlt(this.likelihoods, this.refBase, altAllele);
        double tumorLod = Math.log10(tumorRefAlt[1] + tumorRefAlt[2]) - tumorRefAlt[0];
        return tumorLod;
    }

    public double getRefVsAlt(char altAllele) {
        double[] refAlt = extractRefAlt(this.likelihoods, this.refBase, altAllele);
        double normalLod = refAlt[0] - Math.log10(refAlt[1] + refAlt[2]);
        return normalLod;
    }

    /**
     * Extract the LOD comparing ref:ref to ref:alt and alt:alt
     */
    private double[] extractRefAlt(GenotypeLikelihoods gl, char ref, char altAllele) {
        double refRef = 0;
        double altRef = 0;
        double altAlt = 0;
        for(int j=0; j<10; j++) {
            String gt = gl.genotypes[j];
            double likelihood = gl.likelihoods[j];

            // the ref:mutant theory
            if ( (gt.charAt(0) == ref && gt.charAt(1) == altAllele) ||
                 (gt.charAt(0) == altAllele && gt.charAt(1) == ref) ) {
                altRef += Math.pow(10, likelihood);
            }

            if ( gt.charAt(0) == altAllele && gt.charAt(1) == altAllele) {
                altAlt += Math.pow(10, likelihood);
            }

            if ( gt.charAt(0) == ref && gt.charAt(1) == ref) {
                refRef = likelihood;
            }

        }
        return new double[]{refRef, altRef, altAlt};
    }

    private GenotypeLikelihoods getLikelihood(int locusOffset) {
        GenotypeLikelihoods likelihoods = new GenotypeLikelihoods();


        for(int i=0; i<reads.size(); i++) {
            SAMRecord read = reads.get(i);
            int readOffset = offsets.get(i);

            int offset = readOffset + locusOffset;
            if (offset >= 0 && offset < read.getReadString().length()) {

                char base = read.getReadString().charAt(offset);
                byte qual = read.getBaseQualities()[offset];

                likelihoods.add(refBase, base, qual);
            }

        }
        return likelihoods;
    }

    public double[] getNormalizedProbs(int locusOffset, int qThreshold) {
        GenotypeLikelihoods likelihoods = new GenotypeLikelihoods();


        // FIXME: gaddy suggested this "correction" which evidently has a name...
        //likelihoods.add(refBase, refBase, (byte) 20);

        for(int i=0; i<reads.size(); i++) {
            SAMRecord read = reads.get(i);
            int readOffset = offsets.get(i);

            int offset = readOffset + locusOffset;
            if (offset >= 0 && offset < read.getReadString().length()) {

                char base = read.getReadString().charAt(offset);
                byte qual = read.getBaseQualities()[offset];

                if (qual >= qThreshold) {
                    likelihoods.add(refBase, base, qual);
//                        if (locusOffset == -43) System.out.println("\tUSING " + base + " " + qual);
                } else {
//                        if (locusOffset == -43) System.out.println("\tDropping " + base + " " + qual);
                }
            }

        }


        double[] logLikelihood = likelihoods.likelihoods;
        double[] nonLogLikelihood = new double[10];
        double sum = 0;
        for(int i=0; i<10; i++) {
            nonLogLikelihood[i] = Math.pow(10, logLikelihood[i]);
            sum += nonLogLikelihood[i];
        }

        double[] normalizedProbs = new double[10];
        for(int i=0; i<10; i++) {
            normalizedProbs[i] = nonLogLikelihood[i] / sum;
        }


        //quick sanity check
//            sum=0;
//            for(int i=0; i<10; i++) {
//                sum += normalizedProbs[i];
//            }
//            System.out.println("normalized probs = " + sum);


        return normalizedProbs;
    }


}
