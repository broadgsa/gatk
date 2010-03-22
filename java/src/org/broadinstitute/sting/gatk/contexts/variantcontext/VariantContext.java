package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;

import java.util.*;

/**
 * Class VariantContext
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
 * By default, VariantContexts are immutable.  In order to access (in the rare circumstances where you need them)
 * setter routines, you need to create MutableVariantContexts and MutableGenotypes.
 *
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
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref, T));
 * </pre>
 *
 * If you want to create a non-variant site, just put in a single reference allele
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, snpLoc, Arrays.asList(Aref));
 * </pre>
 *
 * A deletion is just as easy:
 *
 * <pre>
 * VariantContext vc = new VariantContext(name, delLoc, Arrays.asList(ATCref, del));
 * </pre>
 *
 * The only 2 things that distinguishes between a insertion and deletion are the reference allele
 * and the location of the variation.  An insertion has a Null reference allele and at least
 * one non-reference Non-Null allele.  Additionally, the location of the insertion is immediately after
 * a 1-bp GenomeLoc (at say 20).
 *
 * <pre>
 * VariantContext vc = new VariantContext("name", insLoc, Arrays.asList(delRef, ATC));
 * </pre>
 *
 * ==== Converting rods and other data structures to VCs ====
 *
 * You can convert many common types into VariantContexts using the general function:
 *
 * <pre>
 * VariantContextAdaptors.convertToVariantContext(name, myObject)
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
 * Genotype g1 = new Genotype(Arrays.asList(Aref, Aref), "g1", 10);
 * Genotype g2 = new Genotype(Arrays.asList(Aref, T), "g2", 10);
 * Genotype g3 = new Genotype(Arrays.asList(T, T), "g3", 10);
 * VariantContext vc = new VariantContext(snpLoc, alleles, Arrays.asList(g1, g2, g3));
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
 *
 * @author depristo
 */
public class VariantContext {
    protected InferredGeneticContext commonInfo = null;
    public final static double NO_NEG_LOG_10PERROR = InferredGeneticContext.NO_NEG_LOG_10PERROR;

    /** The location of this VariantContext */
    private GenomeLoc loc;

    /** The type (cached for performance reasons) of this context */
    protected Type type = null;

    /** A set of the alleles segregating in this context */
    protected Set<Allele> alleles = null;

    /** A mapping from sampleName -> genotype objects for all genotypes associated with this context */
    protected Map<String, Genotype> genotypes = null;

    /** Counts for each of the possible Genotype types in this context */
    protected int[] genotypeCounts = null;

    protected final static Map<String, Genotype> NO_GENOTYPES = Collections.unmodifiableMap(new HashMap<String, Genotype>());

    // ---------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // ---------------------------------------------------------------------------------------------------------

    // todo move all of attribute object attributes into Map<> and make special filter value for printing out values when
    // emitting VC -> VCF or whatever
    
    /**
     * the complete constructor.  Makes a complete VariantContext from its arguments
     */
    public VariantContext(String name, GenomeLoc loc, Collection<Allele> alleles, Map<String, Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        if ( loc == null ) { throw new StingException("GenomeLoc cannot be null"); }
        this.loc = loc;
        this.commonInfo = new InferredGeneticContext(name, negLog10PError, filters, attributes);

        if ( alleles == null ) { throw new StingException("Alleles cannot be null"); }
        // we need to make this a LinkedHashSet in case the user prefers a given ordering of alleles
        this.alleles = Collections.unmodifiableSet(alleleCollectionToSet(new LinkedHashSet<Allele>(), alleles));

        if ( genotypes == null ) { genotypes = NO_GENOTYPES; }
        this.genotypes = Collections.unmodifiableMap(genotypes);
        calculateGenotypeCounts();

        validate();
    }

    /**
     * Create a new VariantContext
     *
     * @param name
     * @param loc
     * @param alleles
     * @param genotypes
     * @param negLog10PError
     * @param filters
     * @param attributes
     */
    public VariantContext(String name, GenomeLoc loc, Collection<Allele> alleles, Collection<Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        this(name, loc, alleles, genotypes != null ? genotypeCollectionToMap(new HashMap<String, Genotype>(), genotypes) : null, negLog10PError, filters, attributes);
    }

    /**
     * Create a new variant context without genotypes and no Perror, no filters, and no attributes
     * @param name
     * @param loc
     * @param alleles
     */
    public VariantContext(String name, GenomeLoc loc, Collection<Allele> alleles) {
        this(name, loc, alleles, NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    /**
     * Create a new variant context without genotypes and no Perror, no filters, and no attributes
     * @param name
     * @param loc
     * @param alleles
     */
    public VariantContext(String name, GenomeLoc loc, Collection<Allele> alleles, Collection<Genotype> genotypes) {
        this(name, loc, alleles, genotypes, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    /**
     * Copy constructor
     *
     * @param other the VariantContext to copy
     */
    public VariantContext(VariantContext other) {
        this(other.getName(), other.getLocation(), other.getAlleles(), other.getGenotypes(), other.getNegLog10PError(), other.getFilters(), other.getAttributes());
    }


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
        return new VariantContext(getName(), getLocation(), allelesOfGenotypes(genotypes), genotypes, getNegLog10PError(), getFilters(), getAttributes());
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
    }

    /**
     * Determines (if necessary) and returns the type of this variation by examining the alleles it contains.
     *
     * @return the type of this VariantContext
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
    // get routines to access context info fields
    //
    // ---------------------------------------------------------------------------------------------------------
    public String getName()                     { return commonInfo.getName(); }
    public Set<String> getFilters()             { return commonInfo.getFilters(); }
    public boolean isFiltered()                 { return commonInfo.isFiltered(); }
    public boolean isNotFiltered()              { return commonInfo.isNotFiltered(); }
    public boolean hasNegLog10PError()          { return commonInfo.hasNegLog10PError(); }
    public double getNegLog10PError()           { return commonInfo.getNegLog10PError(); }
    public double getPhredScaledQual()          { return commonInfo.getPhredScaledQual(); }

    public Map<String, Object>  getAttributes()  { return commonInfo.getAttributes(); }
    public boolean hasAttribute(String key)     { return commonInfo.hasAttribute(key); }
    public Object getAttribute(String key)      { return commonInfo.getAttribute(key); }

    public Object getAttribute(String key, Object defaultValue) {
        return commonInfo.getAttribute(key, defaultValue);
    }

    public String getAttributeAsString(String key)                        { return commonInfo.getAttributeAsString(key); }
    public String getAttributeAsString(String key, String defaultValue)   { return commonInfo.getAttributeAsString(key, defaultValue); }
    public int getAttributeAsInt(String key)                              { return commonInfo.getAttributeAsInt(key); }
    public int getAttributeAsInt(String key, int defaultValue)            { return commonInfo.getAttributeAsInt(key, defaultValue); }
    public double getAttributeAsDouble(String key)                        { return commonInfo.getAttributeAsDouble(key); }
    public double getAttributeAsDouble(String key, double  defaultValue)  { return commonInfo.getAttributeAsDouble(key, defaultValue); }


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

    /** Private helper routine that grabs the reference allele but doesn't throw an error if there's no such allele */
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
        return Allele.getMatchingAllele(getAlleles(), allele);
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

        return Collections.unmodifiableSet(altAlleles);
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

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with genotypes
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the number of samples in the context
     */
   public int getNSamples() { return genotypes.size(); }

    /**
     * @return true if the context has associated genotypes
     */
    public boolean hasGenotypes() { return genotypes.size() > 0; }

    public boolean hasGenotypes(Collection<String> sampleNames) {
        for ( String name : sampleNames ) {
            if ( ! genotypes.containsKey(name) )
                return false;
        }
        return true;
    }

    /**
     * @return set of all Genotypes associated with this context
     */
    public Map<String, Genotype> getGenotypes() { return genotypes; }

    public List<Genotype> getGenotypesSortedByName() { return Utils.sorted(genotypes); }

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

    public Genotype getGenotype(int ith) {
        return getGenotypesSortedByName().get(ith);
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
     * @return true if it's monomorphic
     */
    public boolean isMonomorphic() {
        return ! isVariant() || getChromosomeCount(getReference()) == getChromosomeCount();
    }

    /**
     * Genotype-specific functions -- are the genotypes polymorphic w.r.t. to the alleles segregating at this
     * site?  That is, is the number of alternate alleles among all fo the genotype > 0?
     *
     * @return true if it's polymorphic
     */
    public boolean isPolymorphic() {
        return ! isMonomorphic();
    }

    private void calculateGenotypeCounts() {
        genotypeCounts = new int[Genotype.Type.values().length];

        for ( Genotype g : getGenotypes().values() ) {
            if ( g.isNoCall() )
                genotypeCounts[Genotype.Type.NO_CALL.ordinal()]++;
            else if ( g.isHomRef() )
                genotypeCounts[Genotype.Type.HOM_REF.ordinal()]++;
            else if ( g.isHet() )
                genotypeCounts[Genotype.Type.HET.ordinal()]++;
            else if ( g.isHomVar() )
                genotypeCounts[Genotype.Type.HOM_VAR.ordinal()]++;
            else
                throw new StingException("Genotype of unknown type: " + g);
        }
    }

    /**
     * Genotype-specific functions -- how many no-calls are there in the genotypes?
     *
     * @return number of no calls
     */
    public int getNoCallCount() {
        return genotypeCounts[Genotype.Type.NO_CALL.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many hom ref calls are there in the genotypes?
     *
     * @return number of hom ref calls
     */
    public int getHomRefCount() {
        return genotypeCounts[Genotype.Type.HOM_REF.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many het calls are there in the genotypes?
     *
     * @return number of het calls
     */
    public int getHetCount() {
        return genotypeCounts[Genotype.Type.HET.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many hom var calls are there in the genotypes?
     *
     * @return number of hom var calls
     */
    public int getHomVarCount() {
        return genotypeCounts[Genotype.Type.HOM_VAR.ordinal()];
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
        return String.format("[VC %s @ %s of type=%s alleles=%s attr=%s GT=%s",
                getName(), getLocation(), this.getType(),
                Utils.sorted(this.getAlleles()), Utils.sortedString(this.getAttributes()), this.getGenotypesSortedByName());
    }

    // protected basic manipulation routines
    private static Set<Allele> alleleCollectionToSet(Set<Allele> dest, Collection<Allele> alleles) {
        for ( Allele a : alleles ) {
            for ( Allele b : dest ) {
                if ( a.basesMatch(b) )
                    throw new IllegalArgumentException("Duplicate allele added to VariantContext: " + a);
            }
            
            dest.add(a);
        }

        return dest;
    }

    private static Map<String, Genotype> genotypeCollectionToMap(Map<String, Genotype> dest, Collection<Genotype> genotypes) {
        for ( Genotype g : genotypes ) {
            if ( dest.containsKey(g.getSampleName() ) )
                throw new IllegalArgumentException("Duplicate genotype added to VariantContext: " + g);
            dest.put(g.getSampleName(), g);
        }

        return dest;
    }

}