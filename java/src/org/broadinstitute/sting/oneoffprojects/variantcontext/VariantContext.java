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
 * == High-level overview ==
 *
 * The VariantContext object is a single general class system for representing genetic variation data composed of:
 *
 * * Allele: representing single genetic haplotypes (A, T, ATC, -)
 * * Genotype: an assignment of alleles for each chromosome of a single named sample at a particular locus
 * * VariantContext: an abstract class holding all segregating alleles at a locus as well as genotypes
 *    for multiple individuals containing alleles at that locus
 *
 * The class system works by defining segregating alleles, creating a variant context representing the segregating
 * information at a locus, and potentially creating and associating genotypes with individuals in the context.
 *
 * All of the classes are highly validating -- call validate() if you modify them -- so you can rely on the
 * self-consistency of the data once you have a VariantContext in hand.  The system has a rich set of assessor
 * and manipulator routines, as well as more complex static support routines in VariantContextUtils.
 *
 * The VariantContext (and Genotype) objects are attributed (supporting addition of arbitrary key/value pairs) and
 * filtered (can represent a variation that is viewed as suspect).
 *
 * VariantContexts are dynamically typed, so whether a VariantContext is a SNP, Indel, or NoVariant depends
 * on the properties of the alleles in the context.  See the detailed documentation on the Type parameter below.
 *
 * It's also easy to create subcontexts based on selected genotypes.
 *
 * == Working with Variant Contexts ==
 * === Some example data ===
 *
 * Allele A, Aref, T, Tref;
 * Allele del, delRef, ATC, ATCref;
 *
 * A [ref] / T at 10
 * GenomeLoc snpLoc = GenomeLocParser.createGenomeLoc("chr1", 10, 10);
 *
 * - / ATC [ref] from 20-23
 * GenomeLoc delLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 22);
 *
 *  // - [ref] / ATC immediately after 20
 * GenomeLoc insLoc = GenomeLocParser.createGenomeLoc("chr1", 20, 20);
 *
 * === Alleles ===
 *
 * See the documentation in the Allele class itself
 *
 * What are they?
 *
 * Alleles can be either reference or non-reference
 *
 * Example alleles used here:
 *
 *   del = new Allele("-");
 *   A = new Allele("A");
 *   Aref = new Allele("A", true);
 *   T = new Allele("T");
 *   ATC = new Allele("ATC");
 *
 * === Creating variant contexts ===
 *
 * ==== By hand ====
 *
 * Here's an example of a A/T polymorphism with the A being reference:
 *
 * <pre>
 * VariantContext vc = new VariantContext(snpLoc, Arrays.asList(Aref, T));
 * </pre>
 *
 * If you want to create a non-variant site, just put in a single reference allele
 *
 * <pre>
 * VariantContext vc = new VariantContext(snpLoc, Arrays.asList(Aref));
 * </pre>
 *
 * A deletion is just as easy:
 *
 * <pre>
 * VariantContext vc = new VariantContext(delLoc, Arrays.asList(ATCref, del));
 * </pre>
 *
 * The only 2 things that distinguishes between a insertion and deletion are the reference allele
 * and the location of the variation.  An insertion has a Null reference allele and at least
 * one non-reference Non-Null allele.  Additionally, the location of the insertion is immediately after
 * a 1-bp GenomeLoc (at say 20).
 *
 * <pre>
 * VariantContext vc = new VariantContext(insLoc, Arrays.asList(delRef, ATC));
 * </pre>
 *
 * ==== Converting rods and other data structures to VCs ====
 *
 * You can convert many common types into VariantContexts using the general function:
 *
 * <pre>
 * VariantContextAdaptors.convertToVariantContext(myObject)
 * </pre>
 *
 * dbSNP and VCFs, for example, can be passed in as myObject and a VariantContext corresponding to that
 * object will be returned.  A null return type indicates that the type isn't yet supported.  This is the best
 * and easiest way to create contexts using RODs.
 *
 *
 * === Working with genotypes ===
 *
 * <pre>
 * List<Allele> alleles = Arrays.asList(Aref, T);
 * VariantContext vc = new VariantContext(snpLoc, alleles);
 *
 * Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "g1", 10);
 * Genotype g2 = new Genotype(Arrays.asList(Aref, T), "g2", 10);
 * Genotype g3 = new Genotype(Arrays.asList(T, T), "g3", 10);
 * vc.addGenotypes(Arrays.asList(g1, g2, g3));
 * </pre>
 *
 * At this point we have 3 genotypes in our context, g1-g3.
 *
 * You can assess a good deal of information about the genotypes through the VariantContext:
 *
 * <pre>
 * vc.hasGenotypes()
 * vc.isMonomorphic()
 * vc.isPolymorphic()
 * vc.getSampleNames().size()
 *
 * vc.getGenotypes()
 * vc.getGenotypes().get("g1")
 * vc.hasGenotype("g1")
 *
 * vc.getChromosomeCount()
 * vc.getChromosomeCount(Aref)
 * vc.getChromosomeCount(T)
 * </pre>
 *
 * === NO_CALL alleles ===
 *
 * The system allows one to create Genotypes carrying special NO_CALL alleles that aren't present in the
 * set of context alleles and that represent undetermined alleles in a genotype:
 *
 * Genotype g4 = new Genotype(Arrays.asList(Allele.NO_CALL, Allele.NO_CALL), "NO_DATA_FOR_SAMPLE", 10);
 *
 *
 * === subcontexts ===
 * It's also very easy get subcontext based only the data in a subset of the genotypes:
 *
 * <pre>
 * VariantContext vc12 = vc.subContextFromGenotypes(Arrays.asList(g1,g2));
 * VariantContext vc1 = vc.subContextFromGenotypes(Arrays.asList(g1));
 * </pre>
 */
public class VariantContext extends AttributedObject {
    /** The location of this VariantContext */
    private GenomeLoc loc;

    /** The type (cached for performance reasons) of this context */
    private Type type = Type.UNDETERMINED;

    /** A set of the alleles segregating in this context */
    private Set<Allele> alleles = new HashSet<Allele>();

    /** A mapping from sampleName -> genotype objects for all genotypes associated with this context */
    private Map<String, Genotype> genotypes = new HashMap<String, Genotype>();


    // ---------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // ---------------------------------------------------------------------------------------------------------
    public VariantContext(GenomeLoc loc) {
        super();

        if ( loc == null ) { throw new StingException("GenomeLoc cannot be null"); }
        this.loc = loc;
    }                                              A

    // todo Add Allele... alleles syntax

    // todo Make root of VariantContext and Genotype immutatable, but extend the system to have MutableVariantContext
    // and MutableGenotype, providing mutation operations.  Then provide .freeze() method on mutable objects to
    // make them immutable

    // todo -- - I'm personally not a huge fan of blanket types like Type.UNDETERMINED, especially when the types are
    // todo public and you have to add 'DON'T USE THIS' to the Javadoc.  Couldn't you use null to represent this type instead?

    // todo wrap collections that are returned directly so that users can't modify them.  return Collections.unmodifiableSet()

    // todo -- rename I'm uncomfortable with the entire AttributedObject system.  The name Object is too general for the application.  The name AttributedObject suggests that it can be used to add attributes to any object, much like Lisp property lists.  However, it contains a negLog10PError member, making the class only applicable to variants.

    // todo move all of attribute object attributes into Map<> and make special filter value for printing out values when
    // emitting VC -> VCF or whatever

    // todo -- Map<String,Object> in attributed object instead of Map<Object,Object>
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
        addGenotypes(genotypes);
        validate();
    }

    public VariantContext(GenomeLoc loc, Collection<Allele> alleles,
                          Collection<Genotype> genotypes, Map<Object, Object> attributes,
                          double negLog10PError, Collection<Object> filters) {
        this(loc);
        setAlleles(alleles);
        addGenotypes(genotypes);
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
        addGenotypes(genotypes);
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

    /**
     * Returns a context identical to this (i.e., filter, qual are all the same) but containing only the Genotype
     * genotype and alleles in genotype.  This is the right way to test if a single genotype is actually
     * variant or not.
     *
     * @param genotype
     * @return
     */
    public VariantContext subContextFromGenotypes(Genotype genotype) {
        return subContextFromGenotypes(Arrays.asList(genotype));
    }


    /**
     * Returns a context identical to this (i.e., filter, qual are all the same) but containing only the Genotypes
     * genotypes and alleles in these genotypes.  This is the right way to test if a single genotype is actually
     * variant or not.
     *
     * @param genotypes
     * @return
     */
    public VariantContext subContextFromGenotypes(Collection<Genotype> genotypes) {
        // todo -- we should check for uniqueness of genotypes
        return new VariantContext(getLocation(), allelesOfGenotypes(genotypes), genotypes, getAttributes(), getNegLog10PError(), getFilters());
    }

    /**
     * helper routnine for subcontext
     * @param genotypes
     * @return
     */
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
     * Also supports NO_VARIATION type, used to indicate that the site isn't polymorphic in the population
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
     * Determines (if necessary) and returns the type of this variation by examining the alleles it contains.
     *
     * @return the type of this VariantContext
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

    /** If this is a BiAlleic SNP, is it a transition? */
    public boolean isTransition()       { return getSNPSubstitutionType() == BaseUtils.BaseSubstitutionType.TRANSITION; } 

    /** If this is a BiAlleic SNP, is it a transversion? */
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

    /** Private helper routine that grabs the reference allele but doesn't through an error if there's no such allele */
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

    /**
     * @return The number of segregating alleles in this context
     */
    public int getNAlleles() {
        return alleles.size();
    }

    /**
     * @return The allele sharing the same bases as this String.  A convenience method; better to use byte[]
     */
    public Allele getAllele(String allele) {
        return getAllele(allele.getBytes());
    }

    /**
     * @return The allele sharing the same bases as this byte[], or null if no such allele is present.
     */
    public Allele getAllele(byte[] allele) {
        for ( Allele a : getAlleles() ) {
            if ( a.basesMatch(allele) ) {
                return a;
            }
        }

        return null;    // couldn't find anything
    }

    /**
     * @return True if this context contains Allele allele, or false otherwise
     */
    public boolean hasAllele(Allele allele) {
        return hasAllele(allele, false);
    }

    public boolean hasAllele(Allele allele, boolean ignoreRefState) {
        for ( Allele a : getAlleles() ) {
            if ( a.equals(allele, ignoreRefState) )
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

    /**
     * @param i -- the ith allele (from 0 to n - 2 for a context with n alleles including a reference allele)
     * @return the ith non-reference allele in this context
     * @throws IllegalArgumentException if i is invalid
     */
    public Allele getAlternateAllele(int i) {
        int n = 0;

        for ( Allele allele : alleles ) {
            if ( allele.isNonReference() && n++ == i )
                return allele;
        }

        throw new IllegalArgumentException("Requested " + i + " alternative allele but there are only " + n + " alternative alleles " + this);
    }

    /**
     * Sets the alleles segregating in this context to the collect of alleles.  Each of which must be unique according
     * to equals() in Allele.  Validate() should be called when you are done modifying the context.
     *
     * @param alleles
     */
    public void setAlleles(Collection<Allele> alleles) {
        this.alleles.clear();
        for ( Allele a : alleles )
            addAllele(a);
    }

    /**
     * Adds allele to the segregating allele list in this context to the collection of alleles.  The new
     * allele must be be unique according to equals() in Allele.
     * Validate() should be called when you are done modifying the context.
     *
     * @param allele
     */
    public void addAllele(Allele allele) {
        final boolean allowDuplicates = false;  // used to be a parameter

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
     * @return true if the context has associated genotypes
     */
    public boolean hasGenotypes() { return genotypes.size() > 0; }

    /**
     * @return set of all Genotypes associated with this context
     */
    public Map<String, Genotype> getGenotypes() { return genotypes; }

    /**
     * Returns a map from sampleName -> Genotype for the genotype associated with sampleName.  Returns a map
     * for consistency with the multi-get function.
     *
     * @param sampleName
     * @return
     * @throws IllegalArgumentException if sampleName isn't bound to a genotype
     */
    public Map<String, Genotype> getGenotypes(String sampleName) {
        return getGenotypes(Arrays.asList(sampleName));
    }

    /**
     * Returns a map from sampleName -> Genotype for each sampleName in sampleNames.  Returns a map
     * for consistency with the multi-get function.
     *
     * @param sampleNames a unique list of sample names
     * @return
     * @throws IllegalArgumentException if sampleName isn't bound to a genotype
     */
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
     * Genotype-specific functions -- are the genotypes monomorphic w.r.t. to the alleles segregating at this
     * site?  That is, is the number of alternate alleles among all fo the genotype == 0?
     *
     * @return
     */
    public boolean isMonomorphic() {
        return ! isVariant() || getChromosomeCount(getReference()) == getChromosomeCount();
    }

    /**
     * Genotype-specific functions -- are the genotypes polymorphic w.r.t. to the alleles segregating at this
     * site?  That is, is the number of alternate alleles among all fo the genotype > 0?
     *
     * @return
     */
    public boolean isPolymorphic() {
        return ! isMonomorphic();
    }

    public void clearGenotypes() {
        this.genotypes.clear();
    }

    /**
     * Adds this single genotype to the context, not allowing duplicate genotypes to be added
     * @param genotype
     */
    public void addGenotypes(Genotype genotype) {
        putGenotype(genotype.getSampleName(), genotype, false);
    }

    /**
     * Adds these genotypes to the context, not allowing duplicate genotypes to be added
     * @param genotypes
     */
    public void addGenotypes(Collection<Genotype> genotypes) {
        for ( Genotype g : genotypes ) {
            addGenotype(g);
        }
    }

    /**
     * Adds these genotype to the context, not allowing duplicate genotypes to be added.
     * @param genotypes
     */
    public void addGenotypes(Map<String, Genotype> genotypes) {

        for ( Map.Entry<String, Genotype> elt : genotypes.entrySet() ) {
            addGenotype(elt.getValue());
        }
    }

    /**
     * Adds these genotypes to the context.
     *
     * @param genotypes
     */
    public void putGenotypes(Map<String, Genotype> genotypes) {
        for ( Map.Entry<String, Genotype> g : genotypes.entrySet() )
            putGenotype(g.getKey(), g.getValue());
    }

    /**
     * Adds these genotypes to the context.
     *
     * @param genotypes
     */
    public void putGenotypes(Collection<Genotype> genotypes) {
        for ( Genotype g : genotypes )
            putGenotype(g);
    }

    /**
     * Adds this genotype to the context, throwing an error if it's already bound.
     *
     * @param genotype
     */
    public void addGenotype(Genotype genotype) {
        addGenotype(genotype.getSampleName(), genotype);
    }

    /**
     * Adds this genotype to the context, throwing an error if it's already bound.
     *
     * @param genotype
     */
    public void addGenotype(String sampleName, Genotype genotype) {
        putGenotype(sampleName, genotype, false);
    }

    /**
     * Adds this genotype to the context.
     *
     * @param genotype
     */
    public void putGenotype(Genotype genotype) {
        putGenotype(genotype.getSampleName(), genotype);
    }

    /**
     * Adds this genotype to the context.
     *
     * @param genotype
     */
    public void putGenotype(String sampleName, Genotype genotype) {
        putGenotype(sampleName, genotype, true);
    }

    private void putGenotype(String sampleName, Genotype genotype, boolean allowOverwrites) {
        if ( hasGenotype(sampleName) && ! allowOverwrites )
            throw new StingException("Attempting to overwrite sample->genotype binding: " + sampleName + " this=" + this);

        if ( ! sampleName.equals(genotype.getSampleName()) )
            throw new StingException("Sample name doesn't equal genotype.getSample(): " + sampleName + " genotype=" + genotype);

        this.genotypes.put(sampleName, genotype);
    }

    /**
     * Removes the binding from sampleName to genotype.  If this doesn't exist, throws an IllegalArgumentException
     * @param sampleName
     */
    public void removeGenotype(String sampleName) {
        if ( ! this.genotypes.containsKey(sampleName) )
            throw new IllegalArgumentException("Sample name isn't contained in genotypes " + sampleName + " genotypes =" + genotypes);

        this.genotypes.remove(sampleName);
    }

    /**
     * Removes genotype from the context.  If this doesn't exist, throws an IllegalArgumentException
     * @param genotype
     */
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

    private boolean validate(boolean throwException) {
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
}