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


import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * A genotype has several key fields
 *
 * -- a sample name, must be a non-null string
 *
 * -- an ordered list of alleles, intrepreted as the genotype of the sample,
 *    each allele for each chromosome given in order.  If alleles = [a*, t]
 *    then the sample is a/t, with a (the reference from the *) the first
 *    chromosome and t on the second chromosome
 *
 * -- a isPhased marker indicting where the alleles are phased with respect to some global
 *    coordinate system.  See VCF4.1 spec for a detailed discussion
 *
 * -- Inline, optimized ints and int[] values for:
 *      -- GQ: the phred-scaled genotype quality, of -1 if it's missing
 *
 *      -- DP: the count of reads at this locus for this sample, of -1 if missing
 *
 *      -- AD: an array of counts of reads at this locus, one for each Allele at the site.
 *             that is, for each allele in the surrounding VariantContext.  Null if missing.
 *
 *      -- PL: phred-scaled genotype likelihoods in standard VCF4.1 order for
 *             all combinations of the alleles in the surrounding VariantContext, given
 *             the ploidy of the sample (from the alleles vector).  Null if missing.
 *
 * -- A general map from String keys to -> Object values for all other attributes in
 *    this genotype.  Note that this map should not contain duplicate values for the
 *    standard bindings for GQ, DP, AD, and PL.  Genotype filters can be put into
 *    this genotype, but it isn't respected by the GATK in analyses
 *
 * The only way to build a Genotype object is with a GenotypeBuilder, which permits values
 * to be set in any order, which means that GenotypeBuilder may at some in the chain of
 * sets pass through invalid states that are not permitted in a fully formed immutable
 * Genotype.
 *
 * Note this is a simplified, refactored Genotype object based on the original
 * generic (and slow) implementation from the original VariantContext + Genotype
 * codebase.
 *
 * @author Mark DePristo
 * @since 05/12
 */
public final class FastGenotype implements Comparable<FastGenotype> {
    private final String sampleName;
    private final List<Allele> alleles;
    private final boolean isPhased;
    private final int GQ;
    private final int DP;
    private final int[] AD;
    private final int[] PL;
    private final Map<String, Object> extendedAttributes;

    private Type type = null;

    /**
     * The only way to make one of these, for use by GenotypeBuilder only
     *
     * @param sampleName
     * @param alleles
     * @param isPhased
     * @param GQ
     * @param DP
     * @param AD
     * @param PL
     * @param extendedAttributes
     */
    @Requires({
            "sampleName != null",
            "alleles != null",
            "GQ >= -1",
            "DP >= -1",
            "validADorPLField(AD)",
            "validADorPLField(PL)",
            "extendedAttributes != null",
            "! hasForbiddenKey(extendedAttributes)"})
    protected FastGenotype(final String sampleName,
                           final List<Allele> alleles,
                           final boolean isPhased,
                           final int GQ,
                           final int DP,
                           final int[] AD,
                           final int[] PL,
                           final Map<String, Object> extendedAttributes) {
        this.sampleName = sampleName;
        this.alleles = alleles;
        this.isPhased = isPhased;
        this.GQ = GQ;
        this.DP = DP;
        this.AD = AD;
        this.PL = PL;
        this.extendedAttributes = extendedAttributes;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Basic getters
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the alleles for this genotype.  Cannot be null.  May be empty
     */
    @Ensures("result != null")
    public List<Allele> getAlleles() {
        return alleles;
    }

    /**
     * Returns how many times allele appears in this genotype object?
     *
     * @param allele
     * @return a value >= 0 indicating how many times the allele occurred in this sample's genotype
     */
    @Requires("allele != null")
    @Ensures("result >= 0")
    public int countAllele(final Allele allele) {
        int c = 0;
        for ( final Allele a : alleles )
            if ( a.equals(allele) )
                c++;

        return c;
    }

    /**
     * Get the ith allele in this genotype
     *
     * @param i the ith allele, must be < the ploidy, starting with 0
     * @return the allele at position i, which cannot be null
     */
    @Requires("i >=0 && i < getPloidy()")
    @Ensures("result != null")
    public Allele getAllele(int i) {
        if ( getType() == Type.UNAVAILABLE )
            throw new ReviewedStingException("Requesting alleles for an UNAVAILABLE genotype");
        return alleles.get(i);
    }

    /**
     * Are the alleles phased w.r.t. the global phasing system?
     *
     * @return true if yes
     */
    public boolean isPhased() {
        return isPhased;
    }

    /**
     * What is the ploidy of this sample?
     *
     * @return the ploidy of this genotype.  0 if the site is no-called.
     */
    @Ensures("result >= 0")
    public int getPloidy() {
        return alleles.size();
    }

    /**
     * @return the sequencing depth of this sample, or -1 if this value is missing
     */
    @Ensures("result >= -1")
    public int getDP() {
        return DP;
    }

    /**
     * @return the count of reads, one for each allele in the surrounding Variant context,
     *      matching the corresponding allele, or null if this value is missing.  MUST
     *      NOT BE MODIFIED!
     */
    public int[] getAD() {
        return AD;
    }

    @Ensures("result != null")
    public String getSampleName() {
        return sampleName;
    }

    /**
     * Returns a phred-scaled quality score, or -1 if none is available
     * @return
     */
    @Ensures("result >= -1")
    public int getGQ()  {
        return GQ;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // The type of this genotype
    //
    // ---------------------------------------------------------------------------------------------------------

    public enum Type {
        /** The sample is no-called (all alleles are NO_CALL */
        NO_CALL,
        /** The sample is homozygous reference */
        HOM_REF,
        /** The sample is heterozygous, with at least one ref and at least one one alt in any order */
        HET,
        /** All alleles are non-reference */
        HOM_VAR,
        /** There is no allele data availble for this sample (alleles.isEmpty) */
        UNAVAILABLE,
        /** Some chromosomes are NO_CALL and others are called */
        MIXED  // no-call and call in the same genotype
    }

    /**
     * @return the high-level type of this sample's genotype
     */
    @Ensures({"type != null", "result != null"})
    public Type getType() {
        if ( type == null ) {
            type = determineType();
        }
        return type;
    }

    /**
     * Internal code to determine the type of the genotype from the alleles vector
     * @return the type
     */
    protected Type determineType() {
        // TODO -- this code is slow and could be optimized for the diploid case
        if ( alleles.isEmpty() )
            return Type.UNAVAILABLE;

        boolean sawNoCall = false, sawMultipleAlleles = false;
        Allele observedAllele = null;

        for ( Allele allele : alleles ) {
            if ( allele.isNoCall() )
                sawNoCall = true;
            else if ( observedAllele == null )
                observedAllele = allele;
            else if ( !allele.equals(observedAllele) )
                sawMultipleAlleles = true;
        }

        if ( sawNoCall ) {
            if ( observedAllele == null )
                return Type.NO_CALL;
            return Type.MIXED;
        }

        if ( observedAllele == null )
            throw new ReviewedStingException("BUG: there are no alleles present in this genotype but the alleles list is not null");

        return sawMultipleAlleles ? Type.HET : observedAllele.isReference() ? Type.HOM_REF : Type.HOM_VAR;
    }

    /**
     * @return true if all observed alleles are the same (regardless of whether they are ref or alt); if any alleles are no-calls, this method will return false.
     */
    public boolean isHom()    { return isHomRef() || isHomVar(); }

    /**
     * @return true if all observed alleles are ref; if any alleles are no-calls, this method will return false.
     */
    public boolean isHomRef() { return getType() == Type.HOM_REF; }

    /**
     * @return true if all observed alleles are alt; if any alleles are no-calls, this method will return false.
     */
    public boolean isHomVar() { return getType() == Type.HOM_VAR; }
    
    /**
     * @return true if we're het (observed alleles differ); if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
     */
    public boolean isHet() { return getType() == Type.HET; }

    /**
     * @return true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF); if any alleles are not no-calls (even if some are), this method will return false.
     */
    public boolean isNoCall() { return getType() == Type.NO_CALL; }

    /**
     * @return true if this genotype is comprised of any alleles that are not no-calls (even if some are).
     */
    public boolean isCalled() { return getType() != Type.NO_CALL && getType() != Type.UNAVAILABLE; }

    /**
     * @return true if this genotype is comprised of both calls and no-calls.
     */
    public boolean isMixed() { return getType() == Type.MIXED; }

    /**
     * @return true if the type of this genotype is set.
     */
    public boolean isAvailable() { return getType() != Type.UNAVAILABLE; }

    // ------------------------------------------------------------------------------
    //
    // methods for getting genotype likelihoods for a genotype object, if present
    //
    // ------------------------------------------------------------------------------

    /**
     * @return Returns true if this Genotype has PL field values
     */
    public boolean hasLikelihoods() {
        return PL != null;
    }

    /**
     * Convenience function that returns a string representation of the PL field of this
     * genotype, or . if none is available.
     *
     * @return a non-null String representation for the PL of this sample
     */
    @Ensures("result != null")
    public String getLikelihoodsString() {
        return hasLikelihoods() ? getLikelihoods().toString() : VCFConstants.MISSING_VALUE_v4;
    }

    /**
     * Returns the GenotypesLikelihoods data associated with this Genotype, or null if missing
     * @return null or a GenotypesLikelihood object for this sample's PL field
     */
    @Ensures({"hasLikelihoods() && result != null", "! hasLikelihoods() && result == null"})
    public GenotypeLikelihoods getLikelihoods() {
        return hasLikelihoods() ? GenotypeLikelihoods.fromPLs(getPL()) : null;
    }

    /**
     * Unsafe low-level accessor the PL field itself, may be null.
     *
     * @return a pointer to the underlying PL data.  MUST NOT BE MODIFIED!
     */
    public int[] getPL() {
        return PL;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Many different string representations
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Return a VCF-like string representation for the alleles of this genotype.
     *
     * Does not append the reference * marker on the alleles.
     *
     * @return a string representing the genotypes, or null if the type is unavailable.
     */
    @Ensures("result != null || ! isAvailable()")
    public String getGenotypeString() {
        return getGenotypeString(true);
    }

    private final static String PHASED_ALLELE_SEPARATOR = "|";
    private final static String UNPHASED_ALLELE_SEPARATOR = "/";

    /**
     * Return a VCF-like string representation for the alleles of this genotype.
     *
     * If ignoreRefState is true, will not append the reference * marker on the alleles.
     *
     * @return a string representing the genotypes, or null if the type is unavailable.
     */
    @Ensures("result != null || ! isAvailable()")
    public String getGenotypeString(boolean ignoreRefState) {
        if ( alleles.size() == 0 )
            return null;

        // Notes:
        // 1. Make sure to use the appropriate separator depending on whether the genotype is phased
        // 2. If ignoreRefState is true, then we want just the bases of the Alleles (ignoring the '*' indicating a ref Allele)
        // 3. So that everything is deterministic with regards to integration tests, we sort Alleles (when the genotype isn't phased, of course)
        return ParsingUtils.join(isPhased() ? PHASED_ALLELE_SEPARATOR : UNPHASED_ALLELE_SEPARATOR,
                ignoreRefState ? getAlleleStrings() : (isPhased() ? getAlleles() : ParsingUtils.sortList(getAlleles())));
    }

    /**
     * Utility that returns a list of allele strings corresponding to the alleles in this sample
     * @return
     */
    private List<String> getAlleleStrings() {
        List<String> al = new ArrayList<String>();
        for ( Allele a : alleles )
            al.add(a.getBaseString());

        return al;
    }

    public String toString() {
        return String.format("[%s %s%s%s%s%s%s]",
                getSampleName(),
                getGenotypeString(false),
                toStringIfExists(VCFConstants.GENOTYPE_QUALITY_KEY, GQ),
                toStringIfExists(VCFConstants.DEPTH_KEY, DP),
                toStringIfExists(VCFConstants.GENOTYPE_ALLELE_DEPTHS, AD),
                toStringIfExists(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, PL),
                sortedString(getExtendedAttributes()));
    }

    public String toBriefString() {
        return String.format("%s:Q%d", getGenotypeString(false), getGQ());
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Comparison operations
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * comparable genotypes -> compareTo on the sample names
     * @param genotype
     * @return
     */
    @Override
    public int compareTo(final FastGenotype genotype) {
        return getSampleName().compareTo(genotype.getSampleName());
    }

    public boolean sameGenotype(final FastGenotype other) {
        return sameGenotype(other, true);
    }

    public boolean sameGenotype(final FastGenotype other, boolean ignorePhase) {
        if (getPloidy() != other.getPloidy())
            return false; // gotta have the same number of allele to be equal

        // By default, compare the elements in the lists of alleles, element-by-element
        Collection<Allele> thisAlleles = this.getAlleles();
        Collection<Allele> otherAlleles = other.getAlleles();

        if (ignorePhase) { // do not care about order, only identity of Alleles
            thisAlleles = new TreeSet<Allele>(thisAlleles);   //implemented Allele.compareTo()
            otherAlleles = new TreeSet<Allele>(otherAlleles);
        }

        return thisAlleles.equals(otherAlleles);
    }

    // ---------------------------------------------------------------------------------------------------------
    // 
    // get routines for extended attributes
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Returns the extended attributes for this object
     * @return is never null, but is often isEmpty()
     */
    @Ensures("result != null")
    public Map<String, Object> getExtendedAttributes() {
        return extendedAttributes;
    }

    /**
     * Is key associated with a value (even a null one) in the extended attributes?
     *
     * Note this will not return true for the inline attributes DP, GQ, AD, or PL
     *
     * @param key a non-null string key to check for an association
     * @return true if key has a value in the extendedAttributes
     */
    @Requires("key != null")
    public boolean hasAttribute(final String key) {
        return extendedAttributes.containsKey(key);
    }

    /**
     * Get the extended attribute value associated with key, if possible
     *
     * @param key a non-null string key to fetch a value for
     * @param defaultValue the value to return if key isn't in the extended attributes
     * @return a value (potentially) null associated with key, or defaultValue if no association exists
     */
    @Requires("key != null")
    @Ensures("hasAttribute(key) || result == defaultValue")
    public Object getAttribute(final String key, final Object defaultValue) {
        return hasAttribute(key) ? extendedAttributes.get(key) : defaultValue;
    }

    // TODO -- add getAttributesAsX interface here

    // ------------------------------------------------------------------------------
    //
    // private utilities
    //
    // ------------------------------------------------------------------------------

    /**
     * a utility method for generating sorted strings from a map key set.
     * @param c the map
     * @param <T> the key type
     * @param <V> the value type
     * @return a sting, enclosed in {}, with comma seperated key value pairs in order of the keys
     */
    @Requires("c != null")
    private static <T extends Comparable<T>, V> String sortedString(Map<T, V> c) {

        // NOTE -- THIS IS COPIED FROM GATK UTILS TO ALLOW US TO KEEP A SEPARATION BETWEEN THE GATK AND VCF CODECS
        final List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);

        final List<String> pairs = new ArrayList<String>();
        for (final T k : t) {
            pairs.add(k + "=" + c.get(k));
        }

        return "{" + ParsingUtils.join(", ", pairs.toArray(new String[pairs.size()])) + "}";
    }

    /**
     * Returns a display name for field name with value v if this isn't -1.  Otherwise returns ""
     * @param name of the field ("AD")
     * @param v the value of the field, or -1 if missing
     * @return a non-null string for display if the field is not missing
     */
    @Requires("name != null")
    @Ensures("result != null")
    private final static String toStringIfExists(final String name, final int v) {
        return v == -1 ? "" : " " + name + " " + v;
    }

    /**
     * Returns a display name for field name with values vs if this isn't null.  Otherwise returns ""
     * @param name of the field ("AD")
     * @param vs the value of the field, or null if missing
     * @return a non-null string for display if the field is not missing
     */
    @Requires("name != null")
    @Ensures("result != null")
    private final static String toStringIfExists(final String name, final int[] vs) {
        if ( vs == null )
            return "";
        else {
            StringBuilder b = new StringBuilder();
            b.append(" ").append(name).append(" ");
            for ( int i = 0; i < vs.length; i++ ) {
                if ( i != 0 ) b.append(",");
                b.append(vs[i]);
            }
            return b.toString();
        }
    }

    /**
     * A list of genotype field keys corresponding to values we
     * manage inline in the Genotype object.  They must not appear in the
     * extended attributes map
     */
    private final static Collection<String> FORBIDDEN_KEYS = Arrays.asList(
            VCFConstants.GENOTYPE_KEY,
            VCFConstants.GENOTYPE_QUALITY_KEY,
            VCFConstants.DEPTH_KEY,
            VCFConstants.GENOTYPE_ALLELE_DEPTHS,
            VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY);

    /**
     * Does the attribute map have a mapping involving a forbidden key (i.e.,
     * one that's managed inline by this Genotypes object?
     *
     * @param attributes the extended attributes key
     * @return
     */
    private final static boolean hasForbiddenKey(final Map<String, Object> attributes) {
        for ( final String forbidden : FORBIDDEN_KEYS )
            if ( attributes.containsKey(forbidden) )
                return true;
        return false;
    }

    /**
     * Is values a valid AD or PL field
     * @param values
     * @return
     */
    private final static boolean validADorPLField(final int[] values) {
        if ( values != null )
            for ( int v : values )
                if ( v < 0 )
                    return false;
        return true;
    }
}