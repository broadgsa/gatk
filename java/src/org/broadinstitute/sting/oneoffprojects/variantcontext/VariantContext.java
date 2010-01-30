package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.*;


/**
 * @author depristo
 *         <p/>
 *         Class VariantContext
 *         <p/>
 *
 * 
 */
public class VariantContext extends AttributedObject {
    private GenomeLoc loc;
    private Type type = Type.UNDETERMINED;
    private Set<Allele> alleles = new HashSet<Allele>();
    private Map<String, Genotype> genotypes = new HashMap<String, Genotype>();
    private Set<Object> filters = new HashSet<Object>();

    // ---------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // ---------------------------------------------------------------------------------------------------------

    public VariantContext(GenomeLoc loc) {
        super();

        if ( loc == null ) { throw new StingException("GenomeLoc cannot be null"); }
        this.loc = loc;
    }

    public VariantContext(GenomeLoc loc, Collection<Allele> alleles) {
        this(loc, alleles, (Map<Object, Object>)null);
    }

    public VariantContext(GenomeLoc loc, Collection<Allele> alleles, Map<Object, Object> attributes) {
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

    public VariantContext(GenomeLoc loc, Collection<Allele> alleles, Collection<Genotype> genotypes) {
        this(loc);
        setAlleles(alleles);
        setGenotypes(genotypes);
        validate();
    }

    public VariantContext(GenomeLoc loc, Collection<Allele> alleles,
                          Collection<Genotype> genotypes, Map<Object, Object> attributes,
                          double negLog10PError, Collection<Object> filters) {
        this(loc);
        setAlleles(alleles);
        setGenotypes(genotypes);
        setAttributes(attributes);
        setNegLog10PError(negLog10PError);
        setFilters(filters);
        validate();
    }

    public VariantContext(GenomeLoc loc, Collection<Allele> alleles,
                          Map<String, Genotype> genotypes, Map<Object, Object> attributes,
                          double negLog10PError, Collection<Object> filters) {
        this(loc);
        setAlleles(alleles);
        setGenotypes(genotypes);
        setAttributes(attributes);
        setNegLog10PError(negLog10PError);
        setFilters(filters);
        validate();
    }

    // todo -- add clone method

    // ---------------------------------------------------------------------------------------------------------
    //
    // Selectors
    //
    // ---------------------------------------------------------------------------------------------------------

    public VariantContext subContextFromGenotypes(Genotype genotype) {
        return subContextFromGenotypes(Arrays.asList(genotype));
    }

    public VariantContext subContextFromGenotypes(Collection<Genotype> genotypes) {
        // todo -- we should check for uniqueness of genotypes
        return subContextFromGenotypes(new HashSet<Genotype>(genotypes), getAttributes());
    }

    public VariantContext subContextFromGenotypes(Collection<Genotype> genotypes, Map<Object, Object> attributes) {
        return new VariantContext(getLocation(), allelesOfGenotypes(genotypes), genotypes, attributes, getNegLog10PError(), getFilters());
    }

    /** helper routnine for subcontext */
    private Set<Allele> allelesOfGenotypes(Collection<Genotype> genotypes) {
        Set<Allele> alleles = new HashSet<Allele>();

        boolean addedref = false;
        for ( Genotype g : genotypes ) {
            for ( Allele a : g.getAlleles() ) {
                addedref = addedref || a.isReference();
                if ( a.isCalled() )
                    alleles.add(a);
            }
        }
        if ( ! addedref ) alleles.add(getReference());

        return alleles;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Filter
    //
    // ---------------------------------------------------------------------------------------------------------

    public Set<Object> getFilters() {
        return filters;
    }

    public boolean isFiltered() {
        return filters.size() > 0;
    }

    public boolean isNotFiltered() {
        return ! isFiltered();
    }

    public void addFilter(Object filter) {
        if ( filter == null ) throw new IllegalArgumentException("BUG: Attempting to add null filter " + this);
        if ( getFilters().contains(filter) ) throw new IllegalArgumentException("BUG: Attempting to add duplicate filter " + filter + " at " + this);
        filters.add(filter);
    }

    public void addFilters(Collection<? extends Object> filters) {
        if ( filters == null ) throw new IllegalArgumentException("BUG: Attempting to add null filters at" + this);
        for ( Object f : filters )
            addFilter(f);
    }

    public void clearFilters() {
        filters.clear();
    }

    public void setFilters(Collection<? extends Object> filters) {
        clearFilters();
        addFilters(filters);
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
        MIXED,

        // special types
        UNDETERMINED // do not use this value, it is reserved for the VariantContext itself
    }

    /**
     * convenience method for switching over the allele type
     *
     * @return the AlleleType of this allele
     **/
    public Type getType() {
        if ( type == Type.UNDETERMINED )
            determineType();

        return type;
    }

    /**
     * convenience method for SNPs
     *
     * @return true if this is a SNP, false otherwise
     */
    public boolean isSNP() { return getType() == Type.SNP; }

    public BaseUtils.BaseSubstitutionType getSNPSubstitutionType() {
        if ( ! isSNP() || ! isBiallelic() ) throw new IllegalStateException("Requested SNP substitution type for bialleic non-SNP " + this);
        return BaseUtils.SNPSubstitutionType(getReference().getBases()[0], getAlternateAllele(0).getBases()[0]);
    }

    public boolean isTransition()       { return getSNPSubstitutionType() == BaseUtils.BaseSubstitutionType.TRANSITION; } 
    public boolean isTransversion()     { return getSNPSubstitutionType() == BaseUtils.BaseSubstitutionType.TRANSVERSION; }

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
     * @return true if the alleles indicate a simple insertion (i.e., the reference allele is Null)
     */
    public boolean isInsertion() {
        return getType() == Type.INDEL && getReference().isNull();
    }

    /**
     * @return true if the alleles indicate a simple deletion (i.e., a single alt allele that is Null)
     */
    public boolean isDeletion() {
        return getType() == Type.INDEL && ! isInsertion();
    }

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

    public Allele getAllele(String allele) {
        return getAllele(allele.getBytes());
    }

    public Allele getAllele(byte[] allele) {
        for ( Allele a : getAlleles() ) {
            if ( a.basesMatch(allele) ) {
                return a;
            }
        }

        return null;    // couldn't find anything
    }

    public boolean hasAllele(Allele allele) {
        for ( Allele a : getAlleles() ) {
            if ( a.equals(allele) )
                return true;
        }

        return false;
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

    public Allele getAlternateAllele(int i) {
        int n = 0;

        for ( Allele allele : alleles ) {
            if ( allele.isNonReference() && n++ == i )
                return allele;
        }

        throw new IllegalArgumentException("Requested " + i + " alternative allele but there are only " + n + " alternative alleles " + this);
    }


    public void setAlleles(Collection<Allele> alleles) {
        this.alleles.clear();
        for ( Allele a : alleles )
            addAllele(a);
    }

    public void addAllele(Allele allele) {
        addAllele(allele, false);
    }

    public void addAllele(Allele allele, boolean allowDuplicates) {
        type = Type.UNDETERMINED;

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

    public boolean hasSingleSample() { return genotypes.size() == 1; }

    /**
     * @return set of all Genotypes associated with this context
     */
    public Map<String, Genotype> getGenotypes() { return genotypes; }

    public Map<String, Genotype> getGenotypes(String sampleName) {
        return getGenotypes(Arrays.asList(sampleName));
    }

    public Map<String, Genotype> getGenotypes(Collection<String> sampleNames) {
        HashMap<String, Genotype> map = new HashMap<String, Genotype>();

        for ( String name : sampleNames ) {
            if ( map.containsKey(name) ) throw new IllegalArgumentException("Duplicate names detected in requested samples " + sampleNames);
            map.put(name, getGenotype(name));
        }

        return map;
    }

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
        int n = 0;

        for ( Genotype g : getGenotypes().values() ) {
            n += g.isNoCall() ? 0 : g.getPloidy();
        }

        return n;
    }

    /**
     * Returns the number of chromosomes carrying allele A in the genotypes
     *
     * @param a
     * @return
     */
    public int getChromosomeCount(Allele a) {
        int n = 0;

        for ( Genotype g : getGenotypes().values() ) {
            n += g.getAlleles(a).size();
        }

        return n;
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
            addGenotype(g.getSampleName(), g);
        }
    }

    public void setGenotypes(Map<String, Genotype> genotypes) {
        this.genotypes.clear();

        for ( Map.Entry<String, Genotype> elt : genotypes.entrySet() ) {
            addGenotype(elt.getKey(), elt.getValue());
        }
    }

    public void addGenotypes(Map<String, Genotype> genotypes) {
        for ( Map.Entry<String, Genotype> g : genotypes.entrySet() )
            addGenotype(g.getKey(), g.getValue());
    }


    public void addGenotypes(Collection<Genotype> genotypes) {
        for ( Genotype g : genotypes )
            addGenotype(g);
    }

    public void addGenotype(Genotype genotype) {
        addGenotype(genotype.getSampleName(), genotype, false);
    }

    public void addGenotype(String sampleName, Genotype genotype) {
        addGenotype(sampleName, genotype, false);
    }

    public void addGenotype(String sampleName, Genotype genotype, boolean allowOverwrites) {
        if ( hasGenotype(sampleName) && ! allowOverwrites )
            throw new StingException("Attempting to overwrite sample->genotype binding: " + sampleName + " this=" + this);

        if ( ! sampleName.equals(genotype.getSampleName()) )
            throw new StingException("Sample name doesn't equal genotype.getSample(): " + sampleName + " genotype=" + genotype);

        this.genotypes.put(sampleName, genotype);
    }

    public void removeGenotype(String sampleName) {
        if ( ! this.genotypes.containsKey(sampleName) )
            throw new IllegalArgumentException("Sample name isn't contained in genotypes " + sampleName + " genotypes =" + genotypes);

        this.genotypes.remove(sampleName);
    }

    public void removeGenotype(Genotype genotype) {
        removeGenotype(genotype.getSampleName());
    }

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
        try {
            validateAlleles();
            validateGenotypes();
        } catch ( IllegalArgumentException e ) {
            if ( throwException )
                throw e;
            else
                return false;
        }

        return true;
    }

    private void validateAlleles() {
        // check alleles
        boolean alreadySeenRef = false, alreadySeenNull = false;
        for ( Allele allele : alleles ) {
            // make sure there's only one reference allele
            if ( allele.isReference() ) {
                if ( alreadySeenRef ) throw new IllegalArgumentException("BUG: Received two reference tagged alleles in VariantContext " + alleles + " this=" + this);
                alreadySeenRef = true;
            }

            if ( allele.isNoCall() ) {
                throw new IllegalArgumentException("BUG: Cannot add a no call allele to a variant context " + alleles + " this=" + this);
            }

            // make sure there's only one null allele
            if ( allele.isNull() ) {
                if ( alreadySeenNull ) throw new IllegalArgumentException("BUG: Received two null alleles in VariantContext " + alleles + " this=" + this);
                alreadySeenNull = true;
            }
        }

        // make sure there's one reference allele
        if ( ! alreadySeenRef )
            throw new IllegalArgumentException("No reference allele found in VariantContext");

//        if ( getType() == Type.INDEL ) {
//            if ( getReference().length() != (getLocation().size()-1) ) {
        if ( (getReference().isNull() && getLocation().size() != 1 ) ||
                (getReference().isNonNull() && getReference().length() != getLocation().size()) ) {
            throw new IllegalStateException("BUG: GenomeLoc " + getLocation() + " has a size == " + getLocation().size() + " but the variation reference allele has length " + getReference().length() + " this = " + this);
        }
    }

    private void validateGenotypes() {
        if ( this.genotypes == null ) throw new IllegalStateException("Genotypes is null");

        for ( Map.Entry<String, Genotype> elt : this.genotypes.entrySet() ) {
            String name = elt.getKey();
            Genotype g = elt.getValue();

            if ( ! name.equals(g.getSampleName()) ) throw new IllegalStateException("Bound sample name " + name + " does not equal the name of the genotype " + g.getSampleName());

            for ( Allele gAllele : g.getAlleles() ) {
                if ( ! hasAllele(gAllele) && gAllele.isCalled() )
                    throw new IllegalStateException("Allele in genotype " + gAllele + " not in the variant context " + alleles);
            }
        }
    }



    // ---------------------------------------------------------------------------------------------------------
    //
    // utility routines
    //
    // ---------------------------------------------------------------------------------------------------------

    private void determineType() {
        if ( type == Type.UNDETERMINED ) {
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
                getLocation(), this.getType(), this.getAlleles(), this.getAttributes(), this.getGenotypes().values());
    }

    // todo -- move to utils
    /**
     * @param allele  the allele to be queried
     *
     * @return the frequency of the given allele in this context
     */
//    public double getAlleleFrequency(Allele allele) {
//        int alleleCount = 0;
//        int totalCount = 0;
//
//        for ( Genotype g : getGenotypes().values() ) {
//            for ( Allele a : g.getAlleles() ) {
//                totalCount++;
//                if ( allele.equals(a) )
//                    alleleCount++;
//            }
//        }
//
//        return totalCount == 0 ? 0.0 : (double)alleleCount / (double)totalCount;
//    }
}