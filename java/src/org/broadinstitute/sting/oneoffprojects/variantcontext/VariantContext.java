package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.refdata.*;

import java.util.*;
import org.apache.commons.jexl.*;


/**
 * @author ebanks
 *         <p/>
 *         Class VariantContext
 *         <p/>
 *         This class represents a context that unifies one or more variants
 */
public class VariantContext {
    private GenomeLoc loc;

    private Set<Allele> alleles = new HashSet<Allele>();

    private Set<Genotype> genotypes = new HashSet<Genotype>();

    private HashMap<Object, Object> attributes = new HashMap<Object, Object>();

    Type type = null;

    private double negLog10PError = 0.0;    // todo - fixme

    /** Have we checked this VariantContext already? */
    private boolean validatedP = false;

//    public VariantContext(VariationRod rod) {
//
//        // TODO -- VariationRod should eventually require classes to implement toVariationContext()
//        // TODO --  (instead of using a temporary adapter class)
//
//        loc = rod.getLocation();
//        reference = new Allele(Allele.AlleleType.REFERENCE, rod.getReference());
//
//        // TODO -- populate the alleles and genotypes through an adapter
//        alleles = new HashSet<Allele>();
//        genotypes = new HashSet<Genotype>();
//
//        attributes = new HashMap<Object, Object>();
//    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // ---------------------------------------------------------------------------------------------------------

    public VariantContext(GenomeLoc loc) {
        if ( loc == null ) { throw new StingException("GenomeLoc cannot be null"); }

        this.loc = loc;
    }

    protected VariantContext(VariantContext parent, Set<Genotype> genotypes, HashMap<Object, Object> attributes) {
        this(parent.getLocation(), parent.getAlleles(), genotypes, attributes);
    }

    public VariantContext(GenomeLoc loc, Set<Allele> alleles, Set<Genotype> genotypes, HashMap<Object, Object> attributes) {
        this(loc);

        // todo -- add extensive testing here

        // todo -- check that exactly one allele is tagged as reference

        this.alleles = new HashSet<Allele>(alleles);
        this.genotypes = new HashSet<Genotype>(genotypes);
        this.attributes = new HashMap<Object, Object>(attributes);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // type operations
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * see: http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=handbook&part=ch5&rendertype=table&id=ch5.ch5_t3
     *
     * Format:
     * dbSNP variation class
     * Rules for assigning allele classes
     * Sample allele definition
     *
     * Single Nucleotide Polymorphisms (SNPs)a
     *   Strictly defined as single base substitutions involving A, T, C, or G.
     *   A/T
     *
     * Deletion/Insertion Polymorphisms (DIPs)
     *   Designated using the full sequence of the insertion as one allele, and either a fully
     *   defined string for the variant allele or a '-' character to specify the deleted allele.
     *   This class will be assigned to a variation if the variation alleles are of different lengths or
     *   if one of the alleles is deleted ('-').
     *   T/-/CCTA/G
     *
     * No-variation
     *   Reports may be submitted for segments of sequence that are assayed and determined to be invariant
     *   in the sample.
     *   (NoVariation)
     *
     * Mixed
     *   Mix of other classes
     *
     *
     * Not currently supported:
     *
     * Heterozygous sequencea
     * The term heterozygous is used to specify a region detected by certain methods that do not
     * resolve the polymorphism into a specific sequence motif. In these cases, a unique flanking
     * sequence must be provided to define a sequence context for the variation.
     * (heterozygous)
     *
     * Microsatellite or short tandem repeat (STR)
     * Alleles are designated by providing the repeat motif and the copy number for each allele.
     * Expansion of the allele repeat motif designated in dbSNP into full-length sequence will
     * be only an approximation of the true genomic sequence because many microsatellite markers are
     * not fully sequenced and are resolved as size variants only.
     * (CAC)8/9/10/11
     *
     * Named variant
     * Applies to insertion/deletion polymorphisms of longer sequence features, such as retroposon
     * dimorphism for Alu or line elements. These variations frequently include a deletion '-' indicator
     * for the absent allele.
     * (alu) / -
     *
     * Multi-Nucleotide Polymorphism (MNP)
     *   Assigned to variations that are multi-base variations of a single, common length
     *   GGA/AGT
     */

    public enum Type {
        NO_VARIATION,
        SNP,
        INDEL,
        MIXED
    }

    /**
     * convenience method for switching over the allele type
     *
     * @return the AlleleType of this allele
     **/
    public Type getType() {
        if ( type == null )
            determineType();

        return type;
    }

    /**
     * convenience method for SNPs
     *
     * @return true if this is a SNP, false otherwise
     */
    public boolean isSNP() { return getType() == Type.SNP; }

    /**
     * convenience method for variants
     *
     * @return true if this is a variant allele, false if it's reference
     */
    public boolean isVariant() { return getType() != Type.NO_VARIATION; }


    /**
     * convenience method for indels
     *
     * @return true if this is an indel, false otherwise
     */
    public boolean isIndel() { return getType() == Type.INDEL; }


    /**
     * convenience method for indels
     *
     * @return true if this is an mixed variation, false otherwise
     */
    public boolean isMixed() { return getType() == Type.MIXED; }    

    // ---------------------------------------------------------------------------------------------------------
    //
    // Generic accessors
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the location of this context
     */
    public GenomeLoc getLocation() { return loc; }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with alleles
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the reference allele for this context
     */
    public Allele getReference() {
        for ( Allele allele : getAlleles() )
            if ( allele.isReference() )
                return allele;

        throw new StingException("BUG: no reference allele found at " + this);
    }

    /**
     * @return true if the context is strictly bi-allelic
     */
    public boolean isBiallelic() {
        //return getAlternateAlleles().size() == 1;
        return getAlleles().size() == 2;
    }

    /**
     * Gets the alleles.  This method should return all of the alleles present at the location,
     * including the reference allele.  There are no constraints imposed on the ordering of alleles
     * in the set. If the reference is not an allele in this context it will not be included.
     *
     * @return the set of alleles
     */
    public Set<Allele> getAlleles() { return alleles; }

    /**
     * Gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference allele.  There are no constraints imposed on the ordering of alleles
     * in the set.
     *
     * @return the set of alternate alleles
     */
    public Set<Allele> getAlternateAlleles() {
        HashSet<Allele> altAlleles = new HashSet<Allele>();
        for ( Allele allele : alleles ) {
            if ( allele.isNonReference() )
                altAlleles.add(allele);
        }

        return altAlleles;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with genotypes
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return true if the context represents variants with associated genotypes
     */
    public boolean hasGenotypes() { return genotypes.size() > 0; }

    /**
     * @return set of all Genotypes associated with this context
     */

    // todo -- genotypes should really be stored as map, not set
    public Set<Genotype> getGenotypes() { return genotypes; }

    public Map<String, Genotype> getGenotypeMap() {
        HashMap<String, Genotype> map = new HashMap<String, Genotype>();
        for ( Genotype g : genotypes )
            map.put(g.getSample(), g);
        return map;
    }


    /**
     * @return the set of all sample names in this context
     */
    public Set<String> getSampleNames() {
        return getGenotypeMap().keySet();
    }

    /**
     * @param sample  the sample name
     *
     * @return the Genotype associated with the given sample in this context or null if the sample is not in this context
     */
    public Genotype getGenotype(String sample) {
        return getGenotypeMap().get(sample);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with attributes
    //
    // ---------------------------------------------------------------------------------------------------------

    // todo -- refactor into AttributedObject and have VariantContext and Genotype inherit from them
    
    // todo -- define common attributes as enum

    /**
     * Sets the given attribute
     *
     * @param key    the attribute key
     * @param value  the attribute value
     */
    public void putAttribute(Object key, Object value) {
        attributes.put(key, value);
    }

    public void putAttributes(Map<? extends Object, Object> map) {
        attributes.putAll(map);
    }

    public boolean hasAttribute(Object key) {
        return attributes.containsKey(key);
    }

    public int getNumAttributes() {
        return attributes.size();
    }

    /**
     * @param key    the attribute key
     *
     * @return the attribute value for the given key (or null if not set)
     */
    public Object getAttribute(Object key) {
        return attributes.get(key);
    }

    public Object getAttribute(Object key, Object defaultValue) {
        if ( hasAttribute(key) )
            return attributes.get(key);
        else
            return defaultValue;
    }

    public String getAttributeAsString(Object key)      { return (String)getAttribute(key); }
    public int getAttributeAsInt(Object key)            { return (Integer)getAttribute(key); }
    public double getAttributeAsDouble(Object key)      { return (Double)getAttribute(key); }

    public String getAttributeAsString(Object key, String defaultValue)   { return (String)getAttribute(key, defaultValue); }
    public int getAttributeAsInt(Object key, int defaultValue)            { return (Integer)getAttribute(key, defaultValue); }
    public double getAttributeAsDouble(Object key, double defaultValue)   { return (Double)getAttribute(key, defaultValue); }

    /**
     * @return the attribute map
     */
    public Map<Object, Object> getAttributes() {
        return attributes;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // validation
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * To be called by any modifying routines
     */
    private void invalidate() { validatedP = false; }

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(boolean throwException) {
        if ( ! validatedP ) {
            boolean valid = false;
            // todo -- add extensive validation checking here
            if ( valid ) {
                validatedP = valid;
            } else if ( throwException ) {
                throw new StingException(this + " failed validation");
            }

            return valid;
        } else {
            return validatedP;
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // utility routines
    //
    // ---------------------------------------------------------------------------------------------------------

    private void determineType() {
        if ( type == null ) {
            // todo -- figure out the variation type
        }
    }

    // todo -- toString() method

    /**
     * @return true if the context represents point alleles only (i.e. no indels or structural variants)
     */
//    public boolean isPointAllele() {
//        for ( Allele allele : alleles ) {
//            if ( allele.isVariant() && !allele.isSNP() )
//                return false;
//        }
//        return true;
//    }
//

//    /**
//     * @return set of all subclasses within this context
//     */
//    public Set<Object> getSubclasses() {
//        Set<Object> subclasses = new HashSet<Object>();
//        for ( Genotype g : genotypes )
//            subclasses.addAll(g.getAttributes().keySet());
//        return subclasses;
//    }

    /**
     * @param allele  the allele to be queried
     *
     * @return the frequency of the given allele in this context
     */
    public double getAlleleFrequency(Allele allele) {
        int alleleCount = 0;
        int totalCount = 0;

        for ( Genotype g : genotypes ) {
            for ( Allele a : g.getAlleles() ) {
                totalCount++;
                if ( allele.equals(a) )
                    alleleCount++;
            }
        }

        return totalCount == 0 ? 0.0 : (double)alleleCount / (double)totalCount;
    }
}