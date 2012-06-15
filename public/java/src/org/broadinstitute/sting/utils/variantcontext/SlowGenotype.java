/*
 * Copyright (c) 2012, The Broad Institute
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


import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
@Deprecated
public class SlowGenotype extends Genotype {
    protected CommonInfo commonInfo;
    public final static double NO_LOG10_PERROR = CommonInfo.NO_LOG10_PERROR;
    protected List<Allele> alleles = null;
    protected boolean isPhased = false;

    protected SlowGenotype(String sampleName, List<Allele> alleles, double log10PError, Set<String> filters, Map<String, Object> attributes, boolean isPhased, double[] log10Likelihoods) {
        super(sampleName);

        if ( alleles == null || alleles.isEmpty() )
            this.alleles = Collections.emptyList();
        else
            this.alleles = Collections.unmodifiableList(alleles);
        commonInfo = new CommonInfo(sampleName, log10PError, filters, attributes);
        if ( log10Likelihoods != null )
            commonInfo.putAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods));
        this.isPhased = isPhased;
        validate();
    }

    @Override public List<Allele> getAlleles() {
        return alleles;
    }

    @Override public Allele getAllele(int i) {
        if ( getType() == GenotypeType.UNAVAILABLE )
            throw new ReviewedStingException("Requesting alleles for an UNAVAILABLE genotype");
        return alleles.get(i);
    }

    @Override public boolean isPhased() { return isPhased; }

    //
    // Useful methods for getting genotype likelihoods for a genotype object, if present
    //
    @Override public boolean hasLikelihoods() {
        return (commonInfo.hasAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY) && !commonInfo.getAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY).equals(VCFConstants.MISSING_VALUE_v4)) ||
                (commonInfo.hasAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY) && !commonInfo.getAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY).equals(VCFConstants.MISSING_VALUE_v4));
    }

    @Override public GenotypeLikelihoods getLikelihoods() {
        GenotypeLikelihoods x = getLikelihoods(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, true);
        if ( x != null )
            return x;
        else {
            x = getLikelihoods(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, false);
            return x;
        }
    }

    private GenotypeLikelihoods getLikelihoods(String key, boolean asPL) {
        Object x = commonInfo.getAttribute(key);
        if ( x instanceof String ) {
            if ( asPL )
                return GenotypeLikelihoods.fromPLField((String)x);
            else
                return GenotypeLikelihoods.fromGLField((String)x);
        }
        else if ( x instanceof GenotypeLikelihoods ) return (GenotypeLikelihoods)x;
        else return null;
    }

    private final void validate() {
        if ( alleles.size() == 0) return;

        for ( Allele allele : alleles ) {
            if ( allele == null )
                throw new IllegalArgumentException("BUG: allele cannot be null in Genotype");
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // 
    // get routines to access context info fields
    //
    // ---------------------------------------------------------------------------------------------------------
    @Override public List<String> getFilters()    { return new ArrayList<String>(commonInfo.getFilters()); }
    @Override public boolean filtersWereApplied() { return commonInfo.filtersWereApplied(); }
    @Override public boolean hasLog10PError()     { return commonInfo.hasLog10PError(); }
    @Override public double getLog10PError()      { return commonInfo.getLog10PError(); }

    @Override
    public boolean hasExtendedAttribute(String key)     { return commonInfo.hasAttribute(key); }

    @Override
    public Object getExtendedAttribute(String key)      { return commonInfo.getAttribute(key); }

    @Override
    public Object getExtendedAttribute(String key, Object defaultValue) {
        return commonInfo.getAttribute(key, defaultValue); 
    }

//    public String getAttributeAsString(String key, String defaultValue)   { return commonInfo.getAttributeAsString(key, defaultValue); }
//    public int getAttributeAsInt(String key, int defaultValue)            { return commonInfo.getAttributeAsInt(key, defaultValue); }
//    public double getAttributeAsDouble(String key, double  defaultValue)  { return commonInfo.getAttributeAsDouble(key, defaultValue); }
//    public boolean getAttributeAsBoolean(String key, boolean  defaultValue)  { return commonInfo.getAttributeAsBoolean(key, defaultValue); }

    @Override
    public int[] getPL() {
        return hasPL() ? getLikelihoods().getAsPLs() : null;
    }

    @Override
    public boolean hasPL() {
        return hasLikelihoods();
    }

    @Override
    public int getDP() {
        return commonInfo.getAttributeAsInt(VCFConstants.DEPTH_KEY, -1);
    }

    @Override
    public boolean hasDP() {
        return commonInfo.hasAttribute(VCFConstants.DEPTH_KEY);
    }

    @Override
    public int[] getAD() {
        if ( hasAD() ) {
            return (int[])commonInfo.getAttribute(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
        } else
            return null;
    }

    @Override
    public boolean hasAD() {
        return commonInfo.hasAttribute(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
    }

    @Override
    public int getGQ() {
        if ( commonInfo.hasLog10PError() )
            return (int)Math.round(commonInfo.getPhredScaledQual());
        else
            return -1;
    }

    @Override
    public boolean hasGQ() {
        return hasLog10PError();
    }

    @Override
    public Map<String, Object> getExtendedAttributes() {
        final Map<String, Object> ea = new LinkedHashMap<String, Object>(commonInfo.getAttributes());
        for ( final String primary : FastGenotype.PRIMARY_KEYS )
            ea.remove(primary);
        return ea;
    }
}