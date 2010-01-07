package org.broadinstitute.sting.gatk.contexts;

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

    private static final String UNIQUIFIED_SUFFIX = ".unique";

    private Set<Allele> alleles;

    private Set<Genotype> genotypes;

    private Allele reference;

    private GenomeLoc loc;

    private HashMap<Object, Object> attributes;


    public VariantContext(VariationRod rod) {

        // TODO -- VariationRod should eventually require classes to implement toVariationContext()
        // TODO --  (instead of using a temporary adapter class)

        loc = rod.getLocation();
        reference = new Allele(Allele.AlleleType.REFERENCE, rod.getReference());

        // TODO -- populate the alleles and genotypes through an adapter
        alleles = new HashSet<Allele>();
        genotypes = new HashSet<Genotype>();

        attributes = new HashMap<Object, Object>();
    }

    protected VariantContext(GenomeLoc loc, Allele reference, Set<Allele> alleles, Set<Genotype> genotypes, HashMap<Object, Object> attributes) {
        this.loc = loc;
        this.reference = reference;
        this.alleles = new HashSet<Allele>(alleles);
        this.genotypes = new HashSet<Genotype>(genotypes);
        this.attributes = new HashMap<Object, Object>(attributes);
    }

    /**
     * @param other  another variant context
     *
     * throws an exception if there is a collision such that the same sample exists in both contexts
     * @return a context representing the merge of this context and the other provided context
     */
    public VariantContext merge(VariantContext other) {
        return merge(other, false);
    }

    /**
     * @param other            another variant context
     * @param uniquifySamples  if true and there is a collision such that the same sample exists in both contexts,
     *                           the samples will be uniquified(based on their sources);
     *                           otherwise, an exception will be thrown
     *
     * @return a context representing the merge of this context and the other provided context
     */
    public VariantContext merge(VariantContext other, boolean uniquifySamples) {
        if ( !loc.equals(other.getLocation()) )
            throw new IllegalArgumentException("The locations must be identical for two contexts to be merged");

        Set<String> samples = getSampleNames();
        Set<Genotype> Gs = new HashSet<Genotype>(genotypes);

        for ( Genotype g : other.getGenotypes() ) {
            if ( samples.contains(g.getSample()) ) {
                if ( uniquifySamples )
                    g.setSample(g.getSample() + UNIQUIFIED_SUFFIX);
                else
                    throw new IllegalStateException("The same sample name exists in both contexts when attempting to merge");
            }
            Gs.add(g);
        }

        HashMap<Object, Object> attrs = new HashMap<Object, Object>(attributes);
        attrs.putAll(other.getAttributes());

        return createNewContext(Gs, attrs);
    }

    /**
     * @return the location of this context
     */
    public GenomeLoc getLocation() { return loc; }

    /**
     * @return the reference allele for this context
     */
    public Allele getReference() { return reference; }

    /**
     * @return true if the context is variant (i.e. contains a non-reference allele)
     */
    public boolean isVariant() {
        for ( Allele allele : alleles ) {
            if ( allele.isVariant() )
                return true;
        }
        return false;
    }

    /**
     * @return true if the context is strictly bi-allelic
     */
    public boolean isBiallelic() {
        return getAlternateAlleles().size() == 1;
    }

    /**
     * @return true if the context represents point alleles only (i.e. no indels or structural variants)
     */
    public boolean isPointAllele() {
        for ( Allele allele : alleles ) {
            if ( allele.isVariant() && !allele.isSNP() )
                return false;
        }
        return true;
    }

    /**
     * @return the set of all sample names in this context
     */
    public Set<String> getSampleNames() {
        Set<String> samples = new TreeSet<String>();
        for ( Genotype g : genotypes )
            samples.add(g.getSample());
        return samples;
    }

    /**
     * @return true if the context represents variants with associated genotypes
     */
    public boolean hasGenotypes() { return genotypes.size() > 0; }

    /**
     * @return set of all Genotypes associated with this context
     */
    public Set<Genotype> getGenotypes() { return genotypes; }

    /**
     * @param sample  the sample name
     *
     * @return the Genotype associated with the given sample in this context or null if the sample is not in this context
     */
    public Genotype getGenotype(String sample) {
        for ( Genotype g : genotypes ) {
            if ( g.getSample().equals(sample) )
                return g;
        }
        return null;
    }

    /**
     * @return set of all subclasses within this context
     */
    public Set<Object> getSubclasses() {
        Set<Object> subclasses = new HashSet<Object>();
        for ( Genotype g : genotypes )
            subclasses.addAll(g.getAttributes().keySet());
        return subclasses;
    }

    /**
     * @param subclass  the name of a subclass of variants to select
     *
     * @return a subset of this context which selects based on the given subclass
     */
    public VariantContext select(String subclass) {
        HashSet<Genotype> Gs = new HashSet<Genotype>();
        for ( Genotype g : genotypes ) {
            if ( g.getAttribute(subclass) != null )
                Gs.add(g);
        }
        return createNewContext(Gs, attributes);
    }

    /**
     * @param expr  a jexl expression describing how to filter this context
     *
     * @return a subset of this context which is filtered based on the given expression
     */
    public VariantContext filter(String expr) {
        HashSet<Genotype> Gs = new HashSet<Genotype>();
        try {
            Expression filterExpression = ExpressionFactory.createExpression(expr);

            for ( Genotype g : genotypes ) {
                JexlContext jContext = JexlHelper.createContext();
                jContext.setVars(g.getAttributes());
                if ( (Boolean)filterExpression.evaluate(jContext) )
                    Gs.add(g);
            }

        } catch (Exception e) {
            throw new StingException("JEXL error in VariantContext: " + e.getMessage());
        }

        return createNewContext(Gs, attributes);
    }

    /**
     * @return a set of new variant contexts, one for each sample from this context
     */
    public Set<VariantContext> splitBySample() {
        Set<VariantContext> contexts = new HashSet<VariantContext>();
        for ( Genotype g : genotypes ) {
            HashSet<Genotype> gAsSet = new HashSet<Genotype>();
            gAsSet.add(g);
            contexts.add(createNewContext(gAsSet, attributes));
        }
        return contexts;
    }

    /**
     * @param Gs    the set of genotypes from which to create a new context
     * @param attrs the attributes for the new context
     *
     * @return a new context based on the given genotypes
     */
    private VariantContext createNewContext(Set<Genotype> Gs, HashMap<Object, Object> attrs) {
        HashSet<Allele> As = new HashSet<Allele>();
        for ( Genotype g : Gs )
            As.addAll(g.getAlleles());
        return new VariantContext(loc, reference, As, Gs, attrs);
    }

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
            if ( !allele.equals(reference) )
                altAlleles.add(allele);
        }
        return altAlleles;
    }

    /**
     * Sets the given attribute
     *
     * @param key    the attribute key
     * @param value  the attribute value
     */
    public void setAttribute(Object key, Object value) {
        attributes.put(key, value);
    }

    /**
     * @param key    the attribute key
     *
     * @return the attribute value for the given key (or null if not set)
     */
    public Object getAttribute(Object key) {
        return attributes.get(key);
    }

    /**
     * @return the attribute map
     */
    public Map<Object, Object> getAttributes() {
        return attributes;
    }

}