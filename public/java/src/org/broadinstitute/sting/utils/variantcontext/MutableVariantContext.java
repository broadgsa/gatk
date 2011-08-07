package org.broadinstitute.sting.utils.variantcontext;


import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Mutable version of VariantContext
 *
 * @author depristo
 */
public class MutableVariantContext extends VariantContext {
    // ---------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // ---------------------------------------------------------------------------------------------------------

    public MutableVariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, Collection<Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        super(source, contig, start, stop, alleles, genotypes, negLog10PError, filters, attributes);
    }

    public MutableVariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, Map<String, Genotype> genotypes, double negLog10PError, Set<String> filters, Map<String, ?> attributes) {
        super(source, contig, start, stop, alleles, genotypes, negLog10PError, filters, attributes);
    }

    public MutableVariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles) {
        super(source, contig, start, stop, alleles, NO_GENOTYPES, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    public MutableVariantContext(String source, String contig, long start, long stop, Collection<Allele> alleles, Collection<Genotype> genotypes) {
        super(source, contig, start, stop, alleles, genotypes, InferredGeneticContext.NO_NEG_LOG_10PERROR, null, null);
    }

    public MutableVariantContext(VariantContext parent) {
        super(parent.getSource(), parent.contig, parent.start, parent.stop, parent.getAlleles(), parent.getGenotypes(), parent.getNegLog10PError(), parent.getFilters(), parent.getAttributes(), parent.getReferenceBaseForIndel());
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

        type = null;

        for ( Allele a : alleles ) {
            if ( a.basesMatch(allele) && ! allowDuplicates )
                throw new IllegalArgumentException("Duplicate allele added to VariantContext" + this);
        }

        // we are a novel allele
        alleles.add(allele);
    }

    public void clearGenotypes() {
        genotypes = new TreeMap<String, Genotype>();
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
            throw new IllegalStateException("Attempting to overwrite sample->genotype binding: " + sampleName + " this=" + this);

        if ( ! sampleName.equals(genotype.getSampleName()) )
            throw new IllegalStateException("Sample name doesn't equal genotype.getSample(): " + sampleName + " genotype=" + genotype);

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

    // todo -- add replace genotype routine

    // ---------------------------------------------------------------------------------------------------------
    //
    // InferredGeneticContext mutation operators
    //
    // ---------------------------------------------------------------------------------------------------------

    public void setSource(String source)                { commonInfo.setName(source); }
    public void addFilter(String filter)                { commonInfo.addFilter(filter); }
    public void addFilters(Collection<String> filters)  { commonInfo.addFilters(filters); }
    public void clearFilters()                          { commonInfo.clearFilters(); }
    public void setFilters(Collection<String> filters)  { commonInfo.setFilters(filters); }
    public void setAttributes(Map<String, ?> map)       { commonInfo.setAttributes(map); }
    public void clearAttributes()                       { commonInfo.clearAttributes(); }
    public void putAttribute(String key, Object value)  { commonInfo.putAttribute(key, value); }
    public void removeAttribute(String key)             { commonInfo.removeAttribute(key); }
    public void putAttributes(Map<String, ?> map)       { commonInfo.putAttributes(map); }
    public void setNegLog10PError(double negLog10PError) { commonInfo.setNegLog10PError(negLog10PError); }
    public void putAttribute(String key, Object value, boolean allowOverwrites) { commonInfo.putAttribute(key, value, allowOverwrites); }
    public void setID(String id) { putAttribute(ID_KEY, id, true); }
}