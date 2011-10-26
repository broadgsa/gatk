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

package org.broadinstitute.sting.utils.variantcontext;

import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.gatk.io.DirectOutputTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.jgrapht.util.MathUtil;

import java.util.EnumMap;
import java.util.Map;

public class GenotypeLikelihoods {
    public static final boolean CAP_PLS = false;
    public static final int PL_CAP = 255;

    //
    // There are two objects here because we are lazy in creating both representations
    // for this object: a vector of log10 Probs and the PL phred-scaled string.  Supports
    // having one set during initializating, and dynamic creation of the other, if needed
    //
    private double[] log10Likelihoods = null;
    private String likelihoodsAsString_PLs = null;

    public final static GenotypeLikelihoods fromPLField(String PLs) {
        return new GenotypeLikelihoods(PLs);
    }

    public final static GenotypeLikelihoods fromGLField(String GLs) {
        return new GenotypeLikelihoods(parseDeprecatedGLString(GLs));
    }

    public final static GenotypeLikelihoods fromLog10Likelihoods(double[] log10Likelihoods) {
        return new GenotypeLikelihoods(log10Likelihoods);
    }

    //
    // You must use the factory methods now
    //
    protected GenotypeLikelihoods(String asString) {
        likelihoodsAsString_PLs = asString;
    }

    protected GenotypeLikelihoods(double[] asVector) {
        log10Likelihoods = asVector;
    }

    /**
     * Returns the genotypes likelihoods in negative log10 vector format.  pr{AA} = x, this
     * vector returns math.log10(x) for each of the genotypes.  Can return null if the
     * genotype likelihoods are "missing".
     *
     * @return
     */
    public double[] getAsVector() {
        // assumes one of the likelihoods vector or the string isn't null
        if ( log10Likelihoods == null ) {
            // make sure we create the GL string if it doesn't already exist
            log10Likelihoods = parsePLsIntoLikelihoods(likelihoodsAsString_PLs);
        }

        return log10Likelihoods;
    }

    public String toString() {
        return getAsString();
    }

    public String getAsString() {
        if ( likelihoodsAsString_PLs == null ) {
            // todo -- should we accept null log10Likelihoods and set PLs as MISSING?
            if ( log10Likelihoods == null )
                throw new TribbleException("BUG: Attempted to get likelihoods as strings and neither the vector nor the string is set!");
            likelihoodsAsString_PLs = convertLikelihoodsToPLString(log10Likelihoods);
        }

        return likelihoodsAsString_PLs;
    }

    //Return genotype likelihoods as an EnumMap with Genotypes as keys and likelihoods as values
    public EnumMap<Genotype.Type,Double> getAsMap(boolean normalizeFromLog10){
        //Make sure that the log10likelihoods are set
        double[] likelihoods = normalizeFromLog10 ? MathUtils.normalizeFromLog10(getAsVector()) : getAsVector();
        EnumMap<Genotype.Type,Double> likelihoodsMap = new EnumMap<Genotype.Type, Double>(Genotype.Type.class);
        likelihoodsMap.put(Genotype.Type.HOM_REF,likelihoods[Genotype.Type.HOM_REF.ordinal()-1]);
        likelihoodsMap.put(Genotype.Type.HET,likelihoods[Genotype.Type.HET.ordinal()-1]);
        likelihoodsMap.put(Genotype.Type.HOM_VAR, likelihoods[Genotype.Type.HOM_VAR.ordinal() - 1]);
        return likelihoodsMap;
    }

    //Return the neg log10 Genotype Quality (GQ) for the given genotype
    public double getNegLog10GQ(Genotype.Type genotype){

        double qual = Double.NEGATIVE_INFINITY;
        EnumMap<Genotype.Type,Double> likelihoods = getAsMap(false);
        for(Map.Entry<Genotype.Type,Double> likelihood : likelihoods.entrySet()){
            if(likelihood.getKey() == genotype)
                continue;
            if(likelihood.getValue() > qual)
                qual = likelihood.getValue();

        }

        //Quality of the most likely genotype = likelihood(most likely) - likelihood (2nd best)
        qual = likelihoods.get(genotype) - qual;

        //Quality of other genotypes 1-P(G)
        if (qual < 0) {
            double[] normalized = MathUtils.normalizeFromLog10(getAsVector());
            double chosenGenotype = normalized[genotype.ordinal()-1];
            qual = -1.0 * Math.log10(1.0 - chosenGenotype);
        }
        return qual;
    }

    private final static double[] parsePLsIntoLikelihoods(String likelihoodsAsString_PLs) {
        if ( !likelihoodsAsString_PLs.equals(VCFConstants.MISSING_VALUE_v4) ) {
            String[] strings = likelihoodsAsString_PLs.split(",");
            double[] likelihoodsAsVector = new double[strings.length];
            for ( int i = 0; i < strings.length; i++ ) {
                likelihoodsAsVector[i] = Integer.parseInt(strings[i]) / -10.0;
            }
            return likelihoodsAsVector;
        } else
            return null;
    }

    /**
     * Back-compatibility function to read old style GL formatted genotype likelihoods in VCF format
     * @param GLString
     * @return
     */
    private final static double[] parseDeprecatedGLString(String GLString) {
        if ( !GLString.equals(VCFConstants.MISSING_VALUE_v4) ) {
            String[] strings = GLString.split(",");
            double[] likelihoodsAsVector = new double[strings.length];
            for ( int i = 0; i < strings.length; i++ ) {
                likelihoodsAsVector[i] = Double.parseDouble(strings[i]);
            }
            return likelihoodsAsVector;
        }

        return null;
    }

    private final static String convertLikelihoodsToPLString(double[] GLs) {
        if ( GLs == null )
            return VCFConstants.MISSING_VALUE_v4;

        StringBuilder s = new StringBuilder();

        double adjust = Double.NEGATIVE_INFINITY;
        for ( double l : GLs ) adjust = Math.max(adjust, l);

        boolean first = true;
        for ( double l : GLs ) {
            if ( ! first )
                s.append(",");
            else
                first = false;

            long PL = Math.round(-10 * (l - adjust));
            if ( CAP_PLS )
                PL = Math.min(PL, PL_CAP);
            s.append(Long.toString(PL));
        }

        return s.toString();
    }
}
