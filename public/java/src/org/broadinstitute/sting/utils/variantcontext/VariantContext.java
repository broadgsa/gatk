package org.broadinstitute.sting.utils.variantcontext;

import org.broad.tribble.Feature;
import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

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
 * vc.getSamples().size()
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
public class VariantContext implements Feature { // to enable tribble intergration
    protected InferredGeneticContext commonInfo = null;
    public final static double NO_NEG_LOG_10PERROR = InferredGeneticContext.NO_NEG_LOG_10PERROR;
    public final static String UNPARSED_GENOTYPE_MAP_KEY = "_UNPARSED_GENOTYPE_MAP_";
    public final static String UNPARSED_GENOTYPE_PARSER_KEY = "_UNPARSED_GENOTYPE_PARSER_";
    public final static String ID_KEY = "ID";

    private final Byte REFERENCE_BASE_FOR_INDEL;

    public final static Set<String> PASSES_FILTERS = Collections.unmodifiableSet(new LinkedHashSet<String>());

    /** The location of this VariantContext */
    protected String contig;
    protected long start;
    protected long stop;

    /** The type (cached for performance reasons) of this context */
    protected Type type = null;

    /** A set of the alleles segregating in this context */
    final protected List<Allele> alleles;

    /** A mapping from sampleName -> genotype objects for all genotypes associated with this context */
    protected Map<String, Genotype> genotypes = null;

    /** Counts for each of the possible Genotype types in this context */
    protected int[] genotypeCounts = null;

    public final static Map<String, Genotype> NO_GENOTYPES = Collections.unmodifiableMap(new HashMap<String, Genotype>());

    // a fast cached access point to the ref / alt alleles for biallelic case
    private Allele REF = null;

    // set to the alt allele when biallelic, otherwise == null
    private Allele ALT = null;

    // were filters applied?
    private boolean filtersWereAppliedToContext;

    // ---------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // ---------------------------------------------------------------------------------------------------------


    /**
     * the complete constructor.  Makes a complete VariantContext from its arguments
     * This is the only constructor that is able to create indels! DO NOT USE THE OTHER ONES.
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles         alleles
     * @param genotypes       genotypes map
     * @param negLog10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     * @param referenceBaseForIndel   padded reference base
     */
    public VariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, Map<String, Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes, Byte referenceBaseForIndel) {
        this(source, contig, start, stop, alleles, genotypes, negLog10PError, filters, attributes, referenceBaseForIndel, false);
    }

    /**
     * the complete constructor.  Makes a complete VariantContext from its arguments
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles         alleles
     * @param genotypes       genotypes map
     * @param negLog10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     */
    public VariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, Map<String, Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        this(source, contig, start, stop, alleles, genotypes, negLog10PError, filters, attributes, null, false);
    }

    /**
     * Makes a VariantContext from its arguments without parsing the genotypes.
     * Note that this constructor assumes that if there is genotype data, then it's been put into
     * the attributes with the UNPARSED_GENOTYPE_MAP_KEY and that the codec has been added with the
     * UNPARSED_GENOTYPE_PARSER_KEY.  It doesn't validate that this is the case because it's possible
     * that there is no genotype data.
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles         alleles
     * @param negLog10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     * @param referenceBaseForIndel   padded reference base
     */
    public VariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, Byte referenceBaseForIndel) {
        this(source, contig, start, stop, alleles, NO_GENOTYPES, negLog10PError, filters, attributes, referenceBaseForIndel, true);
    }

    /**
     * Create a new VariantContext
     *
     * @param source         source
     * @param contig         the contig
     * @param start          the start base (one based)
     * @param stop           the stop reference base (one based)
     * @param alleles        alleles
     * @param genotypes      genotypes set
     * @param negLog10PError qual
     * @param filters        filters: use null for unfiltered and empty set for passes filters
     * @param attributes     attributes
     */
    public VariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, Collection<Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        this(source, contig, start, stop, alleles, genotypes != null ? genotypeCollectionToMap(new TreeMap<String, Genotype>(), genotypes) : null, negLog10PError, filters, attributes, null, false);
    }

    /**
     * Create a new variant context without genotypes and no Perror, no filters, and no attributes
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles alleles
     */
    public VariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles) {
        this(source, contig, start, stop, alleles, NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null, null, false);
    }

    /**
     * Create a new variant context with genotypes but without Perror, filters, and attributes
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles   alleles
     * @param genotypes genotypes
     */
    public VariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, Collection<Genotype> genotypes) {
        this(source, contig, start, stop, alleles, genotypes, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    /**
     * Copy constructor
     *
     * @param other the VariantContext to copy
     */
    public VariantContext(VariantContext other) {
        this(other.getSource(), other.getChr(), other.getStart(), other.getEnd() , other.getAlleles(), other.getGenotypes(), other.getNegLog10PError(), other.filtersWereApplied() ? other.getFilters() : null, other.getAttributes(), other.REFERENCE_BASE_FOR_INDEL, false);
    }

    /**
     * the actual constructor.  Private access only
     *
     * @param source          source
     * @param contig          the contig
     * @param start           the start base (one based)
     * @param stop            the stop reference base (one based)
     * @param alleles         alleles
     * @param genotypes       genotypes map
     * @param negLog10PError  qual
     * @param filters         filters: use null for unfiltered and empty set for passes filters
     * @param attributes      attributes
     * @param referenceBaseForIndel   padded reference base
     * @param genotypesAreUnparsed    true if the genotypes have not yet been parsed
     */
    private VariantContext(String source, String contig, long start, long stop,
                           Collection<Allele> alleles, Map<String, Genotype> genotypes,
                           double negLog10PError, Set<String> filters, Map<String, ?> attributes,
                           Byte referenceBaseForIndel, boolean genotypesAreUnparsed) {
        if ( contig == null ) { throw new IllegalArgumentException("Contig cannot be null"); }
        this.contig = contig;
        this.start = start;
        this.stop = stop;

        if ( !genotypesAreUnparsed && attributes != null ) {
            if ( attributes.containsKey(UNPARSED_GENOTYPE_MAP_KEY) )
                attributes.remove(UNPARSED_GENOTYPE_MAP_KEY);
            if ( attributes.containsKey(UNPARSED_GENOTYPE_PARSER_KEY) )
                attributes.remove(UNPARSED_GENOTYPE_PARSER_KEY);
        }

        this.commonInfo = new InferredGeneticContext(source, negLog10PError, filters, attributes);
        filtersWereAppliedToContext = filters != null;
        REFERENCE_BASE_FOR_INDEL = referenceBaseForIndel;

        if ( alleles == null ) { throw new IllegalArgumentException("Alleles cannot be null"); }

        // we need to make this a LinkedHashSet in case the user prefers a given ordering of alleles
        this.alleles = makeAlleles(alleles);


        if ( genotypes == null ) { genotypes = NO_GENOTYPES; }
        this.genotypes = Collections.unmodifiableMap(genotypes);

        // cache the REF and ALT alleles
        int nAlleles = alleles.size();
        for ( Allele a : alleles ) {
            if ( a.isReference() ) {
                REF = a;
            } else if ( nAlleles == 2 ) { // only cache ALT when biallelic
                ALT = a;
            }
        }

        validate();
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Partial-cloning routines (because Variant Context is immutable).
    //    Note that we don't call vc.getGenotypes() because that triggers the lazy loading.
    //    Also note that we need to create a new attributes map because it's unmodifiable and the constructor may try to modify it.
    //
    // ---------------------------------------------------------------------------------------------------------

    public static VariantContext modifyGenotypes(VariantContext vc, Map<String, Genotype> genotypes) {
        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, new HashMap<String, Object>(vc.getAttributes()), vc.getReferenceBaseForIndel(), false);
    }

    public static VariantContext modifyLocation(VariantContext vc, String chr, int start, int end) {
        return new VariantContext(vc.getSource(), chr, start, end, vc.getAlleles(), vc.genotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, new HashMap<String, Object>(vc.getAttributes()), vc.getReferenceBaseForIndel(), true);
    }

    public static VariantContext modifyFilters(VariantContext vc, Set<String> filters) {
        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd() , vc.getAlleles(), vc.genotypes, vc.getNegLog10PError(), filters, new HashMap<String, Object>(vc.getAttributes()), vc.getReferenceBaseForIndel(), true);
    }

    public static VariantContext modifyAttributes(VariantContext vc, Map<String, Object> attributes) {
        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), vc.genotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, attributes, vc.getReferenceBaseForIndel(), true);
    }

    public static VariantContext modifyReferencePadding(VariantContext vc, Byte b) {
        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), vc.genotypes, vc.getNegLog10PError(), vc.filtersWereApplied() ? vc.getFilters() : null, vc.getAttributes(), b, true);
    }

    public static VariantContext modifyPErrorFiltersAndAttributes(VariantContext vc, double negLog10PError, Set<String> filters, Map<String, Object> attributes) {
        return new VariantContext(vc.getSource(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), vc.genotypes, negLog10PError, filters, attributes, vc.getReferenceBaseForIndel(), true);
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
     * @param genotype genotype
     * @return vc subcontext
     */
    public VariantContext subContextFromGenotypes(Genotype genotype) {
        return subContextFromGenotypes(Arrays.asList(genotype));
    }


    /**
     * Returns a context identical to this (i.e., filter, qual are all the same) but containing only the Genotypes
     * genotypes and alleles in these genotypes.  This is the right way to test if a single genotype is actually
     * variant or not.
     *
     * @param genotypes genotypes
     * @return vc subcontext
     */
    public VariantContext subContextFromGenotypes(Collection<Genotype> genotypes) {
        return subContextFromGenotypes(genotypes, allelesOfGenotypes(genotypes)) ;
    }

    /**
     * Returns a context identical to this (i.e., filter, qual are all the same) but containing only the Genotypes
     * genotypes.  Also, the resulting variant context will contain the alleles provided, not only those found in genotypes
     *
     * @param genotypes genotypes
     * @param alleles the set of allele segregating alleles at this site.  Must include those in genotypes, but may be more
     * @return vc subcontext
     */
    public VariantContext subContextFromGenotypes(Collection<Genotype> genotypes, Collection<Allele> alleles) {
        return new VariantContext(getSource(), contig, start, stop, alleles, genotypes != null ? genotypeCollectionToMap(new TreeMap<String, Genotype>(), genotypes) : null, getNegLog10PError(), filtersWereApplied() ? getFilters() : null, getAttributes(), getReferenceBaseForIndel());
    }


    /**
     * helper routine for subcontext
     * @param genotypes genotypes
     * @return allele set
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
        MNP,    // a multi-nucleotide polymorphism
        INDEL,
        SYMBOLIC,
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


    /**
     * convenience method for variants
     *
     * @return true if this is a variant allele, false if it's reference
     */
    public boolean isVariant() { return getType() != Type.NO_VARIATION; }

    /**
     * convenience method for point events
     *
     * @return true if this is a SNP or ref site, false if it's an indel or mixed event
     */
    public boolean isPointEvent() { return isSNP() || !isVariant(); }

    /**
     * convenience method for indels
     *
     * @return true if this is an indel, false otherwise
     */
    public boolean isIndel() { return getType() == Type.INDEL; }

    /**
     * @return true if the alleles indicate a simple insertion (i.e., the reference allele is Null)
     */
    public boolean isSimpleInsertion() {
        // can't just call !isSimpleDeletion() because of complex indels
        return getType() == Type.INDEL && getReference().isNull() && isBiallelic();
    }

    /**
     * @return true if the alleles indicate a simple deletion (i.e., a single alt allele that is Null)
     */
    public boolean isSimpleDeletion() {
        // can't just call !isSimpleInsertion() because of complex indels
        return getType() == Type.INDEL && getAlternateAllele(0).isNull() && isBiallelic();
    }

    /**
     * @return true if the alleles indicate neither a simple deletion nor a simple insertion
     */
    public boolean isComplexIndel() {
        return isIndel() && !isSimpleDeletion() && !isSimpleInsertion();
    }

    public boolean isSymbolic() {
        return getType() == Type.SYMBOLIC;
    }

    public boolean isMNP() {
        return getType() == Type.MNP;
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

    public boolean hasID() {
        return commonInfo.hasAttribute(ID_KEY);
    }

    public String getID() {
        return (String)commonInfo.getAttribute(ID_KEY);
    }

    public boolean hasReferenceBaseForIndel() {
        return REFERENCE_BASE_FOR_INDEL != null;
    }

    // the indel base that gets stripped off for indels
    public Byte getReferenceBaseForIndel() {
        return REFERENCE_BASE_FOR_INDEL;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // get routines to access context info fields
    //
    // ---------------------------------------------------------------------------------------------------------
    public String getSource()                   { return commonInfo.getName(); }
    public Set<String> getFilters()             { return commonInfo.getFilters(); }
    public boolean isFiltered()                 { return commonInfo.isFiltered(); }
    public boolean isNotFiltered()              { return commonInfo.isNotFiltered(); }
    public boolean filtersWereApplied()         { return filtersWereAppliedToContext; }
    public boolean hasNegLog10PError()          { return commonInfo.hasNegLog10PError(); }
    public double getNegLog10PError()           { return commonInfo.getNegLog10PError(); }
    public double getPhredScaledQual()          { return commonInfo.getPhredScaledQual(); }

    public Map<String, Object>  getAttributes() { return commonInfo.getAttributes(); }
    public boolean hasAttribute(String key)     { return commonInfo.hasAttribute(key); }
    public Object getAttribute(String key)      { return commonInfo.getAttribute(key); }

    public Object getAttribute(String key, Object defaultValue) {
        return commonInfo.getAttribute(key, defaultValue);
    }

    public String getAttributeAsString(String key, String defaultValue)   { return commonInfo.getAttributeAsString(key, defaultValue); }
    public int getAttributeAsInt(String key, int defaultValue)            { return commonInfo.getAttributeAsInt(key, defaultValue); }
    public double getAttributeAsDouble(String key, double  defaultValue)  { return commonInfo.getAttributeAsDouble(key, defaultValue); }
    public boolean getAttributeAsBoolean(String key, boolean  defaultValue)  { return commonInfo.getAttributeAsBoolean(key, defaultValue); }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with alleles
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the reference allele for this context
     */
    public Allele getReference() {
        Allele ref = REF;
        if ( ref == null )
            throw new IllegalStateException("BUG: no reference allele found at " + this);
        return ref;
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
        if ( allele == REF || allele == ALT ) // optimization for cached cases
            return true;

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
    public List<Allele> getAlleles() { return alleles; }

    /**
     * Gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference allele.  There are no constraints imposed on the ordering of alleles
     * in the set.
     *
     * @return the set of alternate alleles
     */
    public List<Allele> getAlternateAlleles() {
        return alleles.subList(1, alleles.size());
    }

    /**
     * Gets the sizes of the alternate alleles if they are insertion/deletion events, and returns a list of their sizes
     *
     * @return a list of indel lengths ( null if not of type indel or mixed )
     */
    public List<Integer> getIndelLengths() {
        if ( getType() != Type.INDEL && getType() != Type.MIXED ) {
            return null;
        }

        List<Integer> lengths = new ArrayList<Integer>();
        for ( Allele a : getAlternateAlleles() ) {
            lengths.add(a.length() - getReference().length());
        }

        return lengths;
    }

    /**
     * @param i -- the ith allele (from 0 to n - 2 for a context with n alleles including a reference allele)
     * @return the ith non-reference allele in this context
     * @throws IllegalArgumentException if i is invalid
     */
    public Allele getAlternateAllele(int i) {
        return alleles.get(i+1);
    }

    /**
     * @param  other  VariantContext whose alternate alleles to compare against
     * @return true if this VariantContext has the same alternate alleles as other,
     *         regardless of ordering. Otherwise returns false.
     */
    public boolean hasSameAlternateAllelesAs ( VariantContext other ) {
        List<Allele> thisAlternateAlleles = getAlternateAlleles();
        List<Allele> otherAlternateAlleles = other.getAlternateAlleles();

        if ( thisAlternateAlleles.size() != otherAlternateAlleles.size() ) {
            return false;
        }

        for ( Allele allele : thisAlternateAlleles ) {
            if ( ! otherAlternateAlleles.contains(allele) ) {
                return false;
            }
        }

        return true;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Working with genotypes
    //
    // ---------------------------------------------------------------------------------------------------------

    private void loadGenotypes() {
        if ( !hasAttribute(UNPARSED_GENOTYPE_MAP_KEY) ) {
            if ( genotypes == null )
                genotypes = NO_GENOTYPES;
            return;
        }

        Object parserObj = getAttribute(UNPARSED_GENOTYPE_PARSER_KEY);
        if ( parserObj == null || !(parserObj instanceof VCFParser) )
            throw new IllegalStateException("There is no VCF parser stored to unparse the genotype data");
        VCFParser parser = (VCFParser)parserObj;

        Object mapObj = getAttribute(UNPARSED_GENOTYPE_MAP_KEY);
        if ( mapObj == null )
            throw new IllegalStateException("There is no mapping string stored to unparse the genotype data");

        genotypes = parser.createGenotypeMap(mapObj.toString(), new ArrayList<Allele>(alleles), getChr(), getStart());

        commonInfo.removeAttribute(UNPARSED_GENOTYPE_MAP_KEY);
        commonInfo.removeAttribute(UNPARSED_GENOTYPE_PARSER_KEY);

        validateGenotypes();
    }

    /**
     * @return the number of samples in the context
     */
    public int getNSamples() {
        loadGenotypes();
        return genotypes.size();
    }

    /**
     * @return true if the context has associated genotypes
     */
    public boolean hasGenotypes() {
        loadGenotypes();
        return genotypes.size() > 0;
    }

    public boolean hasGenotypes(Collection<String> sampleNames) {
        loadGenotypes();
        for ( String name : sampleNames ) {
            if ( ! genotypes.containsKey(name) )
                return false;
        }
        return true;
    }

    /**
     * @return set of all Genotypes associated with this context
     */
    public Map<String, Genotype> getGenotypes() {
        loadGenotypes();
        return genotypes;
    }

    public List<Genotype> getGenotypesSortedByName() {
        loadGenotypes();
        Collection<Genotype> types = new TreeMap<String,Genotype>(genotypes).values();
        return new ArrayList<Genotype>(types);
    }

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
            final Genotype g = getGenotype(name);
            if ( g != null ) {
                map.put(name, g);
            }
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
     * @return chromosome count
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
     * @param a allele
     * @return chromosome count
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
        return ! isVariant() || (hasGenotypes() && getHomRefCount() + getNoCallCount() == getNSamples());
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
        if ( genotypeCounts == null ) {
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
                    genotypeCounts[Genotype.Type.MIXED.ordinal()]++;
            }
        }
    }

    /**
     * Genotype-specific functions -- how many no-calls are there in the genotypes?
     *
     * @return number of no calls
     */
    public int getNoCallCount() {
        calculateGenotypeCounts();
        return genotypeCounts[Genotype.Type.NO_CALL.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many hom ref calls are there in the genotypes?
     *
     * @return number of hom ref calls
     */
    public int getHomRefCount() {
        calculateGenotypeCounts();
        return genotypeCounts[Genotype.Type.HOM_REF.ordinal()];
    }

    /**
     * Genotype-specific functions -- how many het calls are there in the genotypes?
     *
     * @return number of het calls
     */
    public int getHetCount() {
        calculateGenotypeCounts();
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

    /**
     * Genotype-specific functions -- how many mixed calls are there in the genotypes?
     *
     * @return number of mixed calls
     */
    public int getMixedCount() {
        return genotypeCounts[Genotype.Type.MIXED.ordinal()];
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // validation: extra-strict validation routines for paranoid users
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Run all extra-strict validation tests on a Variant Context object
     *
     * @param reference        the true reference allele
     * @param paddedRefBase    the reference base used for padding indels
     * @param rsIDs            the true dbSNP IDs
     */
    public void extraStrictValidation(Allele reference, Byte paddedRefBase, Set<String> rsIDs) {
        // validate the reference
        validateReferenceBases(reference, paddedRefBase);

        // validate the RS IDs
        validateRSIDs(rsIDs);

        // validate the altenate alleles
        validateAlternateAlleles();

        // validate the AN and AC fields
        validateChromosomeCounts();

        // TODO: implement me
        //checkReferenceTrack();
    }

    public void validateReferenceBases(Allele reference, Byte paddedRefBase) {
        // don't validate if we're a complex event
        if ( !isComplexIndel() && !reference.isNull() && !reference.basesMatch(getReference()) ) {
            throw new TribbleException.InternalCodecException(String.format("the REF allele is incorrect for the record at position %s:%d, fasta says %s vs. VCF says %s", getChr(), getStart(), reference.getBaseString(), getReference().getBaseString()));
        }

        // we also need to validate the padding base for simple indels
        if ( hasReferenceBaseForIndel() && !getReferenceBaseForIndel().equals(paddedRefBase) ) {
            throw new TribbleException.InternalCodecException(String.format("the padded REF base is incorrect for the record at position %s:%d, fasta says %s vs. VCF says %s", getChr(), getStart(), (char)paddedRefBase.byteValue(), (char)getReferenceBaseForIndel().byteValue()));
        }
    }

    public void validateRSIDs(Set<String> rsIDs) {
        if ( rsIDs != null && hasAttribute(VariantContext.ID_KEY) ) {
            for ( String id : getID().split(VCFConstants.ID_FIELD_SEPARATOR) ) {
                if ( id.startsWith("rs") && !rsIDs.contains(id) )
                    throw new TribbleException.InternalCodecException(String.format("the rsID %s for the record at position %s:%d is not in dbSNP", id, getChr(), getStart()));
            }
        }
    }

    public void validateAlternateAlleles() {
        if ( !hasGenotypes() )
            return;

        List<Allele> reportedAlleles = getAlleles();
        Set<Allele> observedAlleles = new HashSet<Allele>();
        observedAlleles.add(getReference());
        for ( Genotype g : getGenotypes().values() ) {
            if ( g.isCalled() )
                observedAlleles.addAll(g.getAlleles());
        }

        if ( reportedAlleles.size() != observedAlleles.size() )
            throw new TribbleException.InternalCodecException(String.format("the ALT allele(s) for the record at position %s:%d do not match what is observed in the per-sample genotypes", getChr(), getStart()));

        int originalSize = reportedAlleles.size();
        // take the intersection and see if things change
        observedAlleles.retainAll(reportedAlleles);
        if ( observedAlleles.size() != originalSize )
            throw new TribbleException.InternalCodecException(String.format("the ALT allele(s) for the record at position %s:%d do not match what is observed in the per-sample genotypes", getChr(), getStart()));
    }

    public void validateChromosomeCounts() {
        Map<String, Object> observedAttrs = calculateChromosomeCounts();

        // AN
        if ( hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) ) {
            int reportedAN = Integer.valueOf(getAttribute(VCFConstants.ALLELE_NUMBER_KEY).toString());
            int observedAN = (Integer)observedAttrs.get(VCFConstants.ALLELE_NUMBER_KEY);
            if ( reportedAN != observedAN )
                throw new TribbleException.InternalCodecException(String.format("the Allele Number (AN) tag is incorrect for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedAN, observedAN));
        }

        // AC
        if ( hasAttribute(VCFConstants.ALLELE_COUNT_KEY) ) {
            List<Integer> observedACs = (List<Integer>)observedAttrs.get(VCFConstants.ALLELE_COUNT_KEY);
            if ( getAttribute(VCFConstants.ALLELE_COUNT_KEY) instanceof List ) {
                Collections.sort(observedACs);
                List reportedACs = (List)getAttribute(VCFConstants.ALLELE_COUNT_KEY);
                Collections.sort(reportedACs);
                if ( observedACs.size() != reportedACs.size() )
                    throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag doesn't have the correct number of values for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedACs.size(), observedACs.size()));
                for (int i = 0; i < observedACs.size(); i++) {
                    if ( Integer.valueOf(reportedACs.get(i).toString()) != observedACs.get(i) )
                        throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag is incorrect for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedACs.get(i), observedACs.get(i)));
                }
            } else {
                if ( observedACs.size() != 1 )
                    throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag doesn't have enough values for the record at position %s:%d", getChr(), getStart()));
                int reportedAC = Integer.valueOf(getAttribute(VCFConstants.ALLELE_COUNT_KEY).toString());
                if ( reportedAC != observedACs.get(0) )
                    throw new TribbleException.InternalCodecException(String.format("the Allele Count (AC) tag is incorrect for the record at position %s:%d, %d vs. %d", getChr(), getStart(), reportedAC, observedACs.get(0)));
            }
        }
    }

    private Map<String, Object> calculateChromosomeCounts() {
        Map<String, Object> attributes = new HashMap<String, Object>();

        attributes.put(VCFConstants.ALLELE_NUMBER_KEY, getChromosomeCount());
        ArrayList<Double> alleleFreqs = new ArrayList<Double>();
        ArrayList<Integer> alleleCounts = new ArrayList<Integer>();

        // if there are alternate alleles, record the relevant tags
        if ( getAlternateAlleles().size() > 0 ) {
            for ( Allele allele : getAlternateAlleles() ) {
                alleleCounts.add(getChromosomeCount(allele));
                alleleFreqs.add((double)getChromosomeCount(allele) / (double)getChromosomeCount());
            }
        }
        // otherwise, set them to 0
        else {
            alleleCounts.add(0);
            alleleFreqs.add(0.0);
        }

        attributes.put(VCFConstants.ALLELE_COUNT_KEY, alleleCounts);
        attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqs);
        return attributes;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // validation: the normal validation routines are called automatically upon creation of the VC
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * To be called by any modifying routines
     */
    private boolean validate() {
        return validate(true);
    }

    private boolean validate(boolean throwException) {
        try {
            validateReferencePadding();
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

    private void validateReferencePadding() {
        if (hasSymbolicAlleles()) // symbolic alleles don't need padding...
            return;

        boolean needsPadding = (getReference().length() == getEnd() - getStart()); // off by one because padded base was removed

        if ( needsPadding && !hasReferenceBaseForIndel() )
            throw new ReviewedStingException("Badly formed variant context at location " + getChr() + ":" + getStart() + "; no padded reference base was provided.");
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
        long length = (stop - start) + 1;
        if ( (getReference().isNull() && length != 1 ) ||
                (getReference().isNonNull() && (length - getReference().length()  > 1))) {
            throw new IllegalStateException("BUG: GenomeLoc " + contig + ":" + start + "-" + stop + " has a size == " + length + " but the variation reference allele has length " + getReference().length() + " this = " + this);
        }
    }

    private void validateGenotypes() {
        if ( this.genotypes == null ) throw new IllegalStateException("Genotypes is null");

        for ( Map.Entry<String, Genotype> elt : this.genotypes.entrySet() ) {
            String name = elt.getKey();
            Genotype g = elt.getValue();

            if ( ! name.equals(g.getSampleName()) ) throw new IllegalStateException("Bound sample name " + name + " does not equal the name of the genotype " + g.getSampleName());

            if ( g.isAvailable() ) {
                for ( Allele gAllele : g.getAlleles() ) {
                    if ( ! hasAllele(gAllele) && gAllele.isCalled() )
                        throw new IllegalStateException("Allele in genotype " + gAllele + " not in the variant context " + alleles);
                }
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
            switch ( getNAlleles() ) {
                case 0:
                    throw new IllegalStateException("Unexpected error: requested type of VariantContext with no alleles!" + this);
                case 1:
                    // note that this doesn't require a reference allele.  You can be monomorphic independent of having a
                    // reference allele
                    type = Type.NO_VARIATION;
                    break;
                default:
                    determinePolymorphicType();
            }
        }
    }

    private void determinePolymorphicType() {
        type = null;

        // do a pairwise comparison of all alleles against the reference allele
        for ( Allele allele : alleles ) {
            if ( allele == REF )
                continue;

            // find the type of this allele relative to the reference
            Type biallelicType = typeOfBiallelicVariant(REF, allele);

            // for the first alternate allele, set the type to be that one
            if ( type == null ) {
                type = biallelicType;
            }
            // if the type of this allele is different from that of a previous one, assign it the MIXED type and quit
            else if ( biallelicType != type ) {
                type = Type.MIXED;
                return;
            }
        }
    }

    private static Type typeOfBiallelicVariant(Allele ref, Allele allele) {
        if ( ref.isSymbolic() )
            throw new IllegalStateException("Unexpected error: encountered a record with a symbolic reference allele");

        if ( allele.isSymbolic() )
            return Type.SYMBOLIC;

        if ( ref.length() == allele.length() ) {
            if ( allele.length() == 1 )
                return Type.SNP;
            else
                return Type.MNP;
        }

        // Important note: previously we were checking that one allele is the prefix of the other.  However, that's not an
        // appropriate check as can be seen from the following example:
        // REF = CTTA and ALT = C,CT,CA
        // This should be assigned the INDEL type but was being marked as a MIXED type because of the prefix check.
        // In truth, it should be absolutely impossible to return a MIXED type from this method because it simply
        // performs a pairwise comparison of a single alternate allele against the reference allele (whereas the MIXED type
        // is reserved for cases of multiple alternate alleles of different types).  Therefore, if we've reached this point
        // in the code (so we're not a SNP, MNP, or symbolic allele), we absolutely must be an INDEL.
        return Type.INDEL;

        // old incorrect logic:
        // if (oneIsPrefixOfOther(ref, allele))
        //     return Type.INDEL;
        // else
        //     return Type.MIXED;
    }

    public String toString() {
        return String.format("[VC %s @ %s of type=%s alleles=%s attr=%s GT=%s",
                getSource(), contig + ":" + (start - stop == 0 ? start : start + "-" + stop), this.getType(),
                ParsingUtils.sortList(this.getAlleles()), ParsingUtils.sortedString(this.getAttributes()), this.getGenotypesSortedByName());
    }

    // protected basic manipulation routines
    private static List<Allele> makeAlleles(Collection<Allele> alleles) {
        final List<Allele> alleleList = new ArrayList<Allele>(alleles.size());

        boolean sawRef = false;
        for ( final Allele a : alleles ) {
            for ( final Allele b : alleleList ) {
                if ( a.equals(b, true) )
                    throw new IllegalArgumentException("Duplicate allele added to VariantContext: " + a);
            }

            // deal with the case where the first allele isn't the reference
            if ( a.isReference() ) {
                if ( sawRef )
                    throw new IllegalArgumentException("Alleles for a VariantContext must contain at most one reference allele: " + alleles);
                alleleList.add(0, a);
                sawRef = true;
            }
            else
                alleleList.add(a);
        }

        if ( alleleList.isEmpty() )
            throw new IllegalArgumentException("Cannot create a VariantContext with an empty allele list");

        if ( alleleList.get(0).isNonReference() )
            throw new IllegalArgumentException("Alleles for a VariantContext must contain at least one reference allele: " + alleles);

        return alleleList;
    }

    public static Map<String, Genotype> genotypeCollectionToMap(Map<String, Genotype> dest, Collection<Genotype> genotypes) {
        for ( Genotype g : genotypes ) {
            if ( dest.containsKey(g.getSampleName() ) )
                throw new IllegalArgumentException("Duplicate genotype added to VariantContext: " + g);
            dest.put(g.getSampleName(), g);
        }

        return dest;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // tribble integration routines -- not for public consumption
    //
    // ---------------------------------------------------------------------------------------------------------
    public String getChr() {
        return contig;
    }

    public int getStart() {
        return (int)start;
    }

    public int getEnd() {
        return (int)stop;
    }

    private boolean hasSymbolicAlleles() {
        for (Allele a: getAlleles()) {
            if (a.isSymbolic()) {
                return true;
            }
        }
        return false;
    }

    public static VariantContext createVariantContextWithPaddedAlleles(VariantContext inputVC, boolean refBaseShouldBeAppliedToEndOfAlleles) {

        // see if we need to pad common reference base from all alleles
        boolean padVC;

        // We need to pad a VC with a common base if the length of the reference allele is less than the length of the VariantContext.
        // This happens because the position of e.g. an indel is always one before the actual event (as per VCF convention).
        long locLength = (inputVC.getEnd() - inputVC.getStart()) + 1;
        if (inputVC.hasSymbolicAlleles())
            padVC = true;
        else if (inputVC.getReference().length() == locLength)
            padVC = false;
        else if (inputVC.getReference().length() == locLength-1)
            padVC = true;
        else throw new IllegalArgumentException("Badly formed variant context at location " + String.valueOf(inputVC.getStart()) +
                    " in contig " + inputVC.getChr() + ". Reference length must be at most one base shorter than location size");

        // nothing to do if we don't need to pad bases
        if (padVC) {

            if ( !inputVC.hasReferenceBaseForIndel() )
                throw new ReviewedStingException("Badly formed variant context at location " + inputVC.getChr() + ":" + inputVC.getStart() + "; no padded reference base is available.");

            Byte refByte = inputVC.getReferenceBaseForIndel();

            List<Allele> alleles = new ArrayList<Allele>();
            Map<String, Genotype> genotypes = new TreeMap<String, Genotype>();

            Map<String, Genotype> inputGenotypes = inputVC.getGenotypes();

            for (Allele a : inputVC.getAlleles()) {
                // get bases for current allele and create a new one with trimmed bases
                if (a.isSymbolic()) {
                    alleles.add(a);
                } else {
                    String newBases;
                    if ( refBaseShouldBeAppliedToEndOfAlleles )
                        newBases = a.getBaseString() + new String(new byte[]{refByte});
                    else
                        newBases = new String(new byte[]{refByte}) + a.getBaseString();
                    alleles.add(Allele.create(newBases,a.isReference()));
                }
            }

            // now we can recreate new genotypes with trimmed alleles
            for (String sample : inputVC.getSampleNames()) {
                Genotype g = inputGenotypes.get(sample);

                List<Allele> inAlleles = g.getAlleles();
                List<Allele> newGenotypeAlleles = new ArrayList<Allele>();
                for (Allele a : inAlleles) {
                    if (a.isCalled()) {
                        if (a.isSymbolic()) {
                            newGenotypeAlleles.add(a);
                        } else {
                            String newBases;
                            if ( refBaseShouldBeAppliedToEndOfAlleles )
                                newBases = a.getBaseString() + new String(new byte[]{refByte});
                            else
                                newBases = new String(new byte[]{refByte}) + a.getBaseString();
                            newGenotypeAlleles.add(Allele.create(newBases,a.isReference()));
                        }
                    }
                    else {
                        // add no-call allele
                        newGenotypeAlleles.add(Allele.NO_CALL);
                    }
                }
                genotypes.put(sample, new Genotype(sample, newGenotypeAlleles, g.getNegLog10PError(),
                        g.getFilters(),g.getAttributes(),g.isPhased()));

            }

            // Do not change the filter state if filters were not applied to this context
            Set<String> inputVCFilters = inputVC.filtersWereAppliedToContext ? inputVC.getFilters() : null;
            return new VariantContext(inputVC.getSource(), inputVC.getChr(), inputVC.getStart(), inputVC.getEnd(), alleles, genotypes, inputVC.getNegLog10PError(), inputVCFilters, inputVC.getAttributes(),refByte);
        }
        else
            return inputVC;

    }

    public ArrayList<Allele> getTwoAllelesWithHighestAlleleCounts() {
        // first idea: get two alleles with highest AC
        int maxAC1 = 0, maxAC2=0,maxAC1ind =0, maxAC2ind = 0;
        int i=0;
        int[] alleleCounts = new int[this.getAlleles().size()];
        ArrayList<Allele> alleleArray = new ArrayList<Allele>();
        for (Allele a:this.getAlleles()) {
            int ac = this.getChromosomeCount(a);
            if (ac >=maxAC1) {
                maxAC1 = ac;
                maxAC1ind = i;
            }
            alleleArray.add(a);
            alleleCounts[i++] = ac;
        }
        // now get second best allele
        for (i=0; i < alleleCounts.length; i++) {
            if (i == maxAC1ind)
                continue;
            if (alleleCounts[i] >= maxAC2) {
                maxAC2 = alleleCounts[i];
                maxAC2ind = i;
            }
        }

        Allele alleleA, alleleB;
        if (alleleArray.get(maxAC1ind).isReference()) {
            alleleA = alleleArray.get(maxAC1ind);
            alleleB = alleleArray.get(maxAC2ind);
        }
        else if  (alleleArray.get(maxAC2ind).isReference()) {
            alleleA = alleleArray.get(maxAC2ind);
            alleleB = alleleArray.get(maxAC1ind);
        } else {
            alleleA = alleleArray.get(maxAC1ind);
            alleleB = alleleArray.get(maxAC2ind);
        }
        ArrayList<Allele> a = new ArrayList<Allele>();
        a.add(alleleA);
        a.add(alleleB);
        return a;
    }
    public Allele getAltAlleleWithHighestAlleleCount() {
        // first idea: get two alleles with highest AC
        Allele best = null;
        int maxAC1 = 0;
        for (Allele a:this.getAlternateAlleles()) {
            int ac = this.getChromosomeCount(a);
            if (ac >=maxAC1) {
                maxAC1 = ac;
                best = a;
            }

        }
        return best;
    }

    public int[] getGLIndecesOfAllele(Allele inputAllele) {
        int[] idxVector = new int[3];
        int numAlleles =  this.getAlleles().size();

        int idxDiag = numAlleles;
        int incr = numAlleles - 1;
        int k=1;
        for (Allele a: getAlternateAlleles()) {
            // multi-allelic approximation, part 1: Ideally
            // for each alt allele compute marginal (suboptimal) posteriors -
            // compute indices for AA,AB,BB for current allele - genotype likelihoods are a linear vector that can be thought of
            // as a row-wise upper triangular matrix of likelihoods.
            // So, for example, with 2 alt alleles, likelihoods have AA,AB,AC,BB,BC,CC.
            // 3 alt alleles: AA,AB,AC,AD BB BC BD CC CD DD

            int idxAA = 0;
            int idxAB = k++;
            // yy is always element on the diagonal.
            // 2 alleles: BBelement 2
            // 3 alleles: BB element  3. CC element 5
            // 4 alleles:
            int idxBB = idxDiag;

            if (a.equals(inputAllele)) {
                idxVector[0] = idxAA;
                idxVector[1] = idxAB;
                idxVector[2]  = idxBB;
                break;
            }
            idxDiag += incr--;
        }
        return idxVector;
    }
}
