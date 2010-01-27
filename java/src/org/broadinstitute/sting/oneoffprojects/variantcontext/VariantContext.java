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
public class VariantContext extends AttributedObject {
    private GenomeLoc loc;

    private Set<Allele> alleles = new HashSet<Allele>();

    private Map<String, Genotype> genotypes = new HashMap<String, Genotype>();

    Type type = null;

    // todo -- add QUAL and FILTER

    //private double negLog10PError = 0.0;    // todo - fixme

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

    // todo -- add more convenience methods
    public VariantContext(GenomeLoc loc, Set<Allele> alleles)   { this(loc, alleles, null); }
    public VariantContext(GenomeLoc loc, List<Allele> alleles ) { this(loc, alleles, null); }

    public VariantContext(GenomeLoc loc, List<Allele> alleles, Map<Object, Object> attributes) {
        this(loc);

        HashSet<Allele> alleleSet = new HashSet<Allele>();
        for ( Allele a : alleles ) {
            if ( alleleSet.contains(a) )
                throw new IllegalArgumentException("List<Alleles> contains duplicate elements " + loc + " " + alleles );
            alleleSet.add(a);
        }

        setAlleles(alleleSet);
        setAttributes(attributes);
        validate();
    }

    public VariantContext(GenomeLoc loc, Set<Allele> alleles, Map<Object, Object> attributes) {
        this(loc);
        setAlleles(alleles);
        setAttributes(attributes);
        validate();
    }                     

    public VariantContext(GenomeLoc loc, Set<Allele> alleles, Set<Genotype> genotypes, Map<Object, Object> attributes) {
        this(loc);
        setAlleles(alleles);
        setGenotypes(genotypes);
        setAttributes(attributes);
        validate();
    }


    public VariantContext(GenomeLoc loc, Set<Allele> alleles, Map<String, Genotype> genotypes, Map<Object, Object> attributes) {
        this(loc);
        setAlleles(alleles);
        setGenotypes(genotypes);
        setAttributes(attributes);
        validate();
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

    // todo -- implement, looking at reference allele
    //public boolean isInsertion() { return getType() == Type.INDEL; }
    //public boolean isDeletion() { return getType() == Type.INDEL; }

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
        Allele ref = getReferenceWithoutError();
        if ( ref == null )
            throw new StingException("BUG: no reference allele found at " + this);
        return ref;
    }

    private Allele getReferenceWithoutError() {
        for ( Allele allele : getAlleles() )
            if ( allele.isReference() )
                return allele;
        return null;
    }


    /**
     * @return true if the context is strictly bi-allelic
     */
    public boolean isBiallelic() {
        return getNAlleles() == 2;
    }

    public int getNAlleles() {
        return alleles.size();
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

    public Allele getAlternateAllele(int count) {
        int n = 0;

        for ( Allele allele : alleles ) {
            if ( allele.isNonReference() && n++ == count )
                return allele;
        }

        throw new IllegalArgumentException("Requested " + count + " alternative allele but there are only " + n + " alternative alleles " + this);
    }


    public void setAlleles(Set<Allele> alleles) {
        this.alleles.clear();
        for ( Allele a : alleles )
            addAllele(a);
    }

    public void addAllele(Allele allele) {
        addAllele(allele, false);
    }

    public void addAllele(Allele allele, boolean allowDuplicates) {
        for ( Allele a : alleles ) {
            if ( a.basesMatch(allele) && ! allowDuplicates )
                throw new IllegalArgumentException("Duplicate allele added to VariantContext" + this);
        }

        // we are a novel allele
        alleles.add(allele);
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
    public Map<String, Genotype> getGenotypes() { return genotypes; }

    /**
     * @return the set of all sample names in this context
     */
    public Set<String> getSampleNames() {
        return getGenotypes().keySet();
    }

    /**
     * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS
     *
     * @return
     */
    public int getChromosomeCount() {
        // todo -- return the number of ! no_call alleles
        return 0;
    }

    /**
     * Returns the number of chromosomes carrying allele A in the genotypes
     *
     * @param a
     * @return
     */
    public int getChromosomeCount(Allele a) {
        // todo -- walk through genotypes and count genotypes with allele
        return 0;
    }

    /**
     * These are genotype-specific functions
     *
     * @return
     */
    public boolean isMonomorphic() {
        return ! isVariant() || getChromosomeCount(getReference()) == getChromosomeCount();
    }

    public boolean isPolymorphic() {
        return ! isMonomorphic();
    }


    /**
     * @param sample  the sample name
     *
     * @return the Genotype associated with the given sample in this context or null if the sample is not in this context
     */
    public Genotype getGenotype(String sample) {
        return getGenotypes().get(sample);
    }

    public boolean hasGenotype(String sample) {
        return getGenotypes().containsKey(sample);
    }

    public void setGenotypes(Genotype genotype) {
        this.genotypes.clear();
        addGenotype(genotype);
    }

    public void setGenotypes(Collection<Genotype> genotypes) {
        this.genotypes.clear();

        for ( Genotype g : genotypes ) {
            addGenotype(g.getSample(), g);
        }
    }

    public void setGenotypes(Map<String, Genotype> genotypes) {
        this.genotypes.clear();

        for ( Map.Entry<String, Genotype> elt : genotypes.entrySet() ) {
            addGenotype(elt.getKey(), elt.getValue());
        }
    }

    public void addGenotype(Genotype genotype) {
        addGenotype(genotype.getSample(), genotype, false);
    }


    public void addGenotype(String sampleName, Genotype genotype) {
        addGenotype(sampleName, genotype, false);
    }

    public void addGenotype(String sampleName, Genotype genotype, boolean allowOverwrites) {
        if ( hasGenotype(sampleName) && ! allowOverwrites )
            throw new StingException("Attempting to overwrite sample->genotype binding: " + sampleName + " this=" + this);

        if ( ! sampleName.equals(genotype.getSample()) )
            throw new StingException("Sample name doesn't equal genotype.getSample(): " + sampleName + " genotype=" + genotype);

        this.genotypes.put(sampleName, genotype);
    }

    public void removeGenotype(String sampleName) {
        this.genotypes.remove(sampleName);
    }

    public void removeGenotype(Genotype genotype) {
        removeGenotype(genotype.getSample());
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with attributes
    //
    // ---------------------------------------------------------------------------------------------------------

    // todo -- define common attributes as enum


    // ---------------------------------------------------------------------------------------------------------
    //
    // validation
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * To be called by any modifying routines
     */
    //private void invalidate() { validatedP = false; }

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(boolean throwException) {
        // todo -- add extensive testing here
        // todo -- check that exactly one allele is tagged as reference
        // todo -- check that there's only one null allele

        try {
            // check alleles
            boolean alreadySeenRef = false, alreadySeenNull = false;
            for ( Allele allele : alleles ) {
                if ( allele.isReference() ) {
                    if ( alreadySeenRef ) throw new IllegalArgumentException("BUG: Received two reference tagged alleles in VariantContext " + alleles + " this=" + this);
                    alreadySeenRef = true;
                }

                if ( allele.isNullAllele() ) {
                    if ( alreadySeenNull ) throw new IllegalArgumentException("BUG: Received two null alleles in VariantContext " + alleles + " this=" + this);
                    alreadySeenNull = true;
                }
            }

            if ( ! alreadySeenRef )
                throw new IllegalArgumentException("No reference allele found in VariantContext");
        } catch ( IllegalArgumentException e ) {
            if ( throwException )
                throw e;
            else
                return false;
        }

        return true;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // utility routines
    //
    // ---------------------------------------------------------------------------------------------------------

    private void determineType() {
        if ( type == null ) {
            if ( alleles.size() == 0 ) {
                throw new StingException("Unexpected requested type of VariantContext with no alleles!" + this);
            } else if ( alleles.size() == 1 ) {
                type = Type.NO_VARIATION;
                // note that this doesn't require a reference allele.  You can be monomorphic independent of having a
                // reference allele
            } else if ( isSNPAllele(alleles) ) {
                type = Type.SNP;
            } else if ( isDIPAllele(alleles) ) {
                type = Type.INDEL;
            } else {
                type = Type.MIXED;
            }
        }
    }

    private static boolean isSNPAllele(Set<Allele> alleles) {
        if ( alleles.size() < 2 )
            return false;

        for ( Allele allele : alleles ) {
            if ( allele.length() != 1 )
                return false;
        }

        return true;
    }

    private static boolean isDIPAllele(Set<Allele> alleles) {
        if ( alleles.size() != 2 )
            return false;

        Iterator<Allele> it = alleles.iterator();
        Allele a1 = it.next();
        Allele a2 = it.next();
        return a1.length() != a2.length();
    }
    

    public String toString() {
        return String.format("[VC @ %s of type=%s alleles=%s attr=%s GT=%s",
                getLocation(), this.getType(), this.getAlleles(), this.getAttributes(), this.getGenotypes());
    }

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

    // todo -- move to utils
    /**
     * @param allele  the allele to be queried
     *
     * @return the frequency of the given allele in this context
     */
    public double getAlleleFrequency(Allele allele) {
        int alleleCount = 0;
        int totalCount = 0;

        for ( Genotype g : getGenotypes().values() ) {
            for ( Allele a : g.getAlleles() ) {
                totalCount++;
                if ( allele.equals(a) )
                    alleleCount++;
            }
        }

        return totalCount == 0 ? 0.0 : (double)alleleCount / (double)totalCount;
    }
}