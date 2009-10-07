package org.broadinstitute.sting.utils.genotype;

import edu.mit.broad.picard.genotype.DiploidGenotype;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.StingException;

import java.util.HashMap;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class LikelyhoodObject
 *         <p/>
 *         An object used to store likelyhood information for genotypes.  Genotype
 *         likelihoods are assumed to be infinite (negitive log likelihood), unless set.
 *         This allows the consumer to make an empty LikelihoodObject, and just set
 *         those values which have associated likelihood values.
 */
public class LikelihoodObject {


    // our possible genotypes, in order according to GLFv3
    public enum GENOTYPE {
        AA, AC, AG, AT, CC, CG, CT, GG, GT, TT
    }

    // our pileup of bases
    //final private String basePileup;

    // possible types of likihoods to store

    public enum LIKELIHOOD_TYPE {
        NEGATIVE_LOG, LOG, RAW;
    }

    // our liklihood storage type
    protected LIKELIHOOD_TYPE mLikelihoodType = LIKELIHOOD_TYPE.NEGATIVE_LOG;

    // default the bestGenotype likelihoods to the allele AA
    protected GENOTYPE bestGenotype = GENOTYPE.AA;

    // how many genotypes we're storing
    public static final int genoTypeCount = GENOTYPE.values().length;

    // the associated negitive log likelihoods for each genotype
    protected final HashMap<GENOTYPE, Double> likelihoods = new HashMap<GENOTYPE, Double>();

    /** create a blank likelihood object */
    public LikelihoodObject() {
        for (GENOTYPE type : GENOTYPE.values()) {
            likelihoods.put(type, Double.MAX_VALUE);
        }
    }

    /**
     * create a likelihood object, given a picard style GenotypeLikelihoods object.  The
     * GenotypeLikelihoods stores likelihoods in log likelihood format, and we want them in
     * negitive log likelihood
     *
     * @param lk the likelihood object
     */
    public LikelihoodObject(GenotypeLikelihoods lk) {
        mLikelihoodType = LIKELIHOOD_TYPE.LOG;
        Double minValue = Double.MAX_VALUE;
        for (GENOTYPE type : GENOTYPE.values()) {
            byte[] bases = new byte[2];
            bases[0] = (byte) type.toString().charAt(0);
            bases[1] = (byte) type.toString().charAt(1);
            double val = -1.0d * lk.getLikelihood(DiploidGenotype.fromBases(bases));
            likelihoods.put(type, val);
            if (val < minValue) {
                bestGenotype = type;
            }
        }
    }

    /**
     * create a likelyhood object, given an array of genotype scores in GLFv3 ordering
     *
     * @param values an array of int's from 0 to 255, representing the negitive log likelihoods.
     * @param type   the likelihood storage type
     */
    public LikelihoodObject(double[] values, LIKELIHOOD_TYPE type) {
        mLikelihoodType = type;
        if (values.length != GENOTYPE.values().length) {
            throw new IllegalArgumentException("invalid array passed to LikelihoodObject, should be size " + GENOTYPE.values().length);
        }
        findBestLikelihood(values);
    }

    /**
     * find the best likelihood
     * @param values
     */
    private void findBestLikelihood(double[] values) {
        int index = 0;
        double lowestScore = Double.MAX_VALUE;
        for (GENOTYPE t : GENOTYPE.values()) {
            likelihoods.put(t, values[index]);
            if (values[index] < lowestScore) {
                lowestScore = values[index];
                bestGenotype = t;
            }
            ++index;
        }
    }

    /**
     * set the likelihood, given it's probability and the genotype
     *
     * @param type the genotype
     * @param lh   the likelihood as a double
     */
    public void setLikelihood(GENOTYPE type, double lh) {
        likelihoods.put(type, lh);
        if (lh < likelihoods.get(this.bestGenotype)) {
            this.bestGenotype = type;
        }
    }

    /**
     * find the minimum likelihood value stored in the set.  This represents the most likely genotype,
     * since genotypes are represented as negitive log likeihoods
     *
     * @return the min value
     */
    public double getBestLikelihood() {
        return likelihoods.get(this.bestGenotype);
    }

    /**
     * return a byte array representation of the likelihood object, in GLFv3 specified order.
     * The return type is short[] instead of byte[], since signed bytes only store -127 to 127,
     * not the 255 range we need.
     *
     * @return a byte array of the genotype values
     */
    public short[] toByteArray() {
        short ret[] = new short[GENOTYPE.values().length];
        int index = 0;
        for (GENOTYPE type : GENOTYPE.values()) {
            ret[index] = (likelihoods.get(type).intValue() > 254) ? 255 : (short) likelihoods.get(type).intValue();
            ++index;
        }
        return ret;
    }

    /**
     * create a float array of our genotype values, in order specified in the GENOTYPE enum (currently the GLF and
     * geli ordering).
     *
     * @return a float array containing our genotype likelihoods, as negitive log likelihoods
     */
    public double[] toDoubleArray() {
        // make an array of floats
        double[] ft = new double[10];
        int index = 0;
        for (GENOTYPE T : GENOTYPE.values()) {
            ft[index] = this.likelihoods.get(T).doubleValue();
            index++;
        }
        return ft;
    }

    /**
     * convert this object, with aditional information, to a GenotypeLikelihoods object.  This involves determining
     * what our underlying storage type is, and coverting our values to the appropriate (log likelihood) format.
     *
     * @return a GenotypeLikelihoods object representing our data
     */
    public GenotypeLikelihoods convertToGenotypeLikelihoods(SAMFileHeader samHeader, int seqIndex, int seqPosition, byte refBase) {
        double[] ft = toDoubleArray();
        float[] db = new float[ft.length];
        int index = 0;
        if (this.mLikelihoodType == LIKELIHOOD_TYPE.NEGATIVE_LOG) {
            for (; index < ft.length; index++) {
                db[index] = ((float) ft[index] * -1.0f);
            }
        } else if (this.mLikelihoodType == LIKELIHOOD_TYPE.RAW) {
            for (; index < ft.length; index++) {
                db[index] = (float) Math.log(ft[index]);
            }
        } else {
            for (int x = 0; x < ft.length; x++)
                db[x] = (float)ft[x];
        }
        return new GenotypeLikelihoods(samHeader, seqIndex, seqPosition, refBase, db);
    }

    /**
     * getter for the likelihood type
     *
     * @return our likelihood storage type
     */
    public LIKELIHOOD_TYPE getLikelihoodType() {
        return mLikelihoodType;
    }


    /**
     * validate a genotype score
     *
     * @param score the score to validate
     */
    public void validateScore(double score) {
        int x = 0;
        switch (mLikelihoodType) {
            case NEGATIVE_LOG:
                if (score < 0)
                    throw new StingException("Likelikhood score of " + score + " is invalid, for NEGATIVE_LOG it must be greater than or equal to 0");
                break;
            case LOG:
                if (score > 0)
                    throw new StingException("Likelikhood score of " + score + " is invalid, for LOG it must be less than or equal to 0");
                break;
            case RAW:
                if (score < 0 || score > 1)
                    throw new StingException("Likelikhood score of " + score + " is invalid, for RAW it must be [0,1]");
                break;
        }
    }


    /**
     * set our likelihood storage type, and adjust our current likelihood values to reflect
     * the new setting.
     *
     * @param likelihood the type to set the values to.
     */
    public void setLikelihoodType(LIKELIHOOD_TYPE likelihood) {
        if (likelihood == mLikelihoodType)
            return;
        if (mLikelihoodType == LIKELIHOOD_TYPE.RAW) {
            double mult = 1.0;
            if (likelihood == LIKELIHOOD_TYPE.NEGATIVE_LOG) {
                mult = -1.0;
            }
            // one of us in log, the other negitive log, it doesn't matter which
            for (GENOTYPE g : likelihoods.keySet()) {
                likelihoods.put(g, -1.0 * Math.log(likelihoods.get(g)));
            }
        } else if (likelihood == LIKELIHOOD_TYPE.RAW) {
            double mult = 1.0;
            if (mLikelihoodType == LIKELIHOOD_TYPE.NEGATIVE_LOG) {
                mult = -1.0;
            }
            // one of us in log, the other negitive log, it doesn't matter which
            for (GENOTYPE g : likelihoods.keySet()) {
                likelihoods.put(g, Math.pow(likelihoods.get(g) * mult, 10));
            }
        } else {
            // one of us in log, the other negitive log, it doesn't matter which
            for (GENOTYPE g : likelihoods.keySet()) {
                likelihoods.put(g, -1.0 * likelihoods.get(g));
            }
        }
        this.mLikelihoodType = likelihood;
    }
}


