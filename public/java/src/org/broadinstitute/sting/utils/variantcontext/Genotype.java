package org.broadinstitute.sting.utils.variantcontext;


import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
public class Genotype {

    public final static String PHASED_ALLELE_SEPARATOR = "|";
    public final static String UNPHASED_ALLELE_SEPARATOR = "/";

    protected InferredGeneticContext commonInfo;
    public final static double NO_NEG_LOG_10PERROR = InferredGeneticContext.NO_NEG_LOG_10PERROR;
    protected List<Allele> alleles = null; // new ArrayList<Allele>();
    protected Type type = null;

    protected boolean isPhased = false;
    protected boolean filtersWereAppliedToContext;

    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased) {
        this(sampleName, alleles, negLog10PError, filters, attributes, isPhased, null);
    }

    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean isPhased, double[] log10Likelihoods) {
        if ( alleles != null )
            this.alleles = Collections.unmodifiableList(alleles);
        commonInfo = new InferredGeneticContext(sampleName, negLog10PError, filters, attributes);
        if ( log10Likelihoods != null )
            commonInfo.putAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, GenotypeLikelihoods.fromLog10Likelihoods(log10Likelihoods));
        filtersWereAppliedToContext = filters != null;
        this.isPhased = isPhased;
        validate();
    }

    /**
     * Creates a new Genotype for sampleName with genotype according to alleles.
     * @param sampleName
     * @param alleles
     * @param negLog10PError the confidence in these alleles
     * @param log10Likelihoods a log10 likelihoods for each of the genotype combinations possible for alleles, in the standard VCF ordering, or null if not known
     */
    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError, double[] log10Likelihoods) {
        this(sampleName, alleles, negLog10PError, null, null, false, log10Likelihoods);
    }

    public Genotype(String sampleName, List<Allele> alleles, double negLog10PError) {
        this(sampleName, alleles, negLog10PError, null, null, false);
    }

    public Genotype(String sampleName, List<Allele> alleles) {
        this(sampleName, alleles, NO_NEG_LOG_10PERROR, null, null, false);
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Partial-cloning routines (because Genotype is immutable).
    //
    // ---------------------------------------------------------------------------------------------------------

    public static Genotype modifyName(Genotype g, String name) {
        return new Genotype(name, g.getAlleles(), g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, g.getAttributes(), g.isPhased());
    }

    public static Genotype modifyAttributes(Genotype g, Map<String, Object> attributes) {
        return new Genotype(g.getSampleName(), g.getAlleles(), g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, attributes, g.isPhased());
    }

    public static Genotype modifyAlleles(Genotype g, List<Allele> alleles) {
        return new Genotype(g.getSampleName(), alleles, g.getNegLog10PError(), g.filtersWereApplied() ? g.getFilters() : null, g.getAttributes(), g.isPhased());
    }

    /**
     * @return the alleles for this genotype
     */
    public List<Allele> getAlleles() {
        return alleles;
    }

    public List<Allele> getAlleles(Allele allele) {
        if ( getType() == Type.UNAVAILABLE )
            throw new ReviewedStingException("Requesting alleles for an UNAVAILABLE genotype");

        List<Allele> al = new ArrayList<Allele>();
        for ( Allele a : alleles )
            if ( a.equals(allele) )
                al.add(a);

        return Collections.unmodifiableList(al);
    }

    public Allele getAllele(int i) {
        if ( getType() == Type.UNAVAILABLE )
            throw new ReviewedStingException("Requesting alleles for an UNAVAILABLE genotype");
        return alleles.get(i);
    }

    public boolean isPhased() { return isPhased; }

    /**
     * @return the ploidy of this genotype
     */
    public int getPloidy() {
        if ( alleles == null )
            throw new ReviewedStingException("Requesting ploidy for an UNAVAILABLE genotype");
        return alleles.size();
    }

    public enum Type {
        NO_CALL,
        HOM_REF,
        HET,
        HOM_VAR,
        UNAVAILABLE,
        MIXED  // no-call and call in the same genotype
    }

    public Type getType() {
        if ( type == null ) {
            type = determineType();
        }
        return type;
    }

    protected Type determineType() {
        if ( alleles == null )
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

    //
    // Useful methods for getting genotype likelihoods for a genotype object, if present
    //
    public boolean hasLikelihoods() {
        return (hasAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY) && !getAttribute(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY).equals(VCFConstants.MISSING_VALUE_v4)) ||
                (hasAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY) && !getAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY).equals(VCFConstants.MISSING_VALUE_v4));
    }
    
    public GenotypeLikelihoods getLikelihoods() {
        GenotypeLikelihoods x = getLikelihoods(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, true);
        if ( x != null )
            return x;
        else {
            x = getLikelihoods(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, false);
            if ( x != null ) return x;
            else
                throw new IllegalStateException("BUG: genotype likelihood field in " + this.getSampleName() + " sample are not either a string or a genotype likelihood class!");
        }
    }

    private GenotypeLikelihoods getLikelihoods(String key, boolean asPL) {
        Object x = getAttribute(key);
        if ( x instanceof String ) {
            if ( asPL )
                return GenotypeLikelihoods.fromPLField((String)x);
            else
                return GenotypeLikelihoods.fromGLField((String)x);
        }
        else if ( x instanceof GenotypeLikelihoods ) return (GenotypeLikelihoods)x;
        else return null;
    }

    public void validate() {
        if ( alleles == null ) return;
        if ( alleles.size() == 0) throw new IllegalArgumentException("BUG: alleles cannot be of size 0");

        // int nNoCalls = 0;
        for ( Allele allele : alleles ) {
            if ( allele == null )
                throw new IllegalArgumentException("BUG: allele cannot be null in Genotype");
            // nNoCalls += allele.isNoCall() ? 1 : 0;
        }

        // Technically, the spec does allow for the below case so this is not an illegal state
        //if ( nNoCalls > 0 && nNoCalls != alleles.size() )
        //    throw new IllegalArgumentException("BUG: alleles include some No Calls and some Calls, an illegal state " + this);
    }

    public String getGenotypeString() {
        return getGenotypeString(true);
    }

    public String getGenotypeString(boolean ignoreRefState) {
        if ( alleles == null )
            return null;

        // Notes:
        // 1. Make sure to use the appropriate separator depending on whether the genotype is phased
        // 2. If ignoreRefState is true, then we want just the bases of the Alleles (ignoring the '*' indicating a ref Allele)
        // 3. So that everything is deterministic with regards to integration tests, we sort Alleles (when the genotype isn't phased, of course)
        return ParsingUtils.join(isPhased() ? PHASED_ALLELE_SEPARATOR : UNPHASED_ALLELE_SEPARATOR,
                ignoreRefState ? getAlleleStrings() : (isPhased() ? getAlleles() : ParsingUtils.sortList(getAlleles())));
    }

    private List<String> getAlleleStrings() {
        List<String> al = new ArrayList<String>();
        for ( Allele a : alleles )
            al.add(a.getBaseString());

        return al;
    }

    public String toString() {
        int Q = (int)Math.round(getPhredScaledQual());
        return String.format("[%s %s Q%s %s]", getSampleName(), getGenotypeString(false),
                Q == -10 ? "." : String.format("%2d",Q), sortedString(getAttributes()));
    }

    public String toBriefString() {
        return String.format("%s:Q%.2f", getGenotypeString(false), getPhredScaledQual());
    }

    public boolean sameGenotype(Genotype other) {
        return sameGenotype(other, true);
    }

    public boolean sameGenotype(Genotype other, boolean ignorePhase) {
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

    /**
     * a utility method for generating sorted strings from a map key set.
     * @param c the map
     * @param <T> the key type
     * @param <V> the value type
     * @return a sting, enclosed in {}, with comma seperated key value pairs in order of the keys
     */
    private static <T extends Comparable<T>, V> String sortedString(Map<T, V> c) {
        // NOTE -- THIS IS COPIED FROM GATK UTILS TO ALLOW US TO KEEP A SEPARATION BETWEEN THE GATK AND VCF CODECS
        List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);

        List<String> pairs = new ArrayList<String>();
        for (T k : t) {
            pairs.add(k + "=" + c.get(k));
        }

        return "{" + ParsingUtils.join(", ", pairs.toArray(new String[pairs.size()])) + "}";
    }


    // ---------------------------------------------------------------------------------------------------------
    // 
    // get routines to access context info fields
    //
    // ---------------------------------------------------------------------------------------------------------
    public String getSampleName()       { return commonInfo.getName(); }
    public Set<String> getFilters()     { return commonInfo.getFilters(); }
    public boolean isFiltered()         { return commonInfo.isFiltered(); }
    public boolean isNotFiltered()      { return commonInfo.isNotFiltered(); }
    public boolean filtersWereApplied() { return filtersWereAppliedToContext; }
    public boolean hasNegLog10PError()  { return commonInfo.hasNegLog10PError(); }
    public double getNegLog10PError()   { return commonInfo.getNegLog10PError(); }
    public double getPhredScaledQual()  { return commonInfo.getPhredScaledQual(); }

    public Map<String, Object> getAttributes()  { return commonInfo.getAttributes(); }
    public boolean hasAttribute(String key)     { return commonInfo.hasAttribute(key); }
    public Object getAttribute(String key)      { return commonInfo.getAttribute(key); }
    
    public Object getAttribute(String key, Object defaultValue) {
        return commonInfo.getAttribute(key, defaultValue); 
    }

    public String getAttributeAsString(String key, String defaultValue)   { return commonInfo.getAttributeAsString(key, defaultValue); }
    public int getAttributeAsInt(String key, int defaultValue)            { return commonInfo.getAttributeAsInt(key, defaultValue); }
    public double getAttributeAsDouble(String key, double  defaultValue)  { return commonInfo.getAttributeAsDouble(key, defaultValue); }
    public boolean getAttributeAsBoolean(String key, boolean  defaultValue)  { return commonInfo.getAttributeAsBoolean(key, defaultValue); }
}