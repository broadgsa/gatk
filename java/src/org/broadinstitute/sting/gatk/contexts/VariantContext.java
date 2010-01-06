package org.broadinstitute.sting.gatk.contexts;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.*;
import org.apache.commons.jexl.*;


/**
 * @author ebanks
 *         <p/>
 *         Interface VariantContext
 *         <p/>
 *         This class represents a context that unifies one or more variants
 */
public interface VariantContext {

    /**
     * @return the location of this context
     */
    public GenomeLoc getLocation();

    /**
     * @return the reference base or bases for this context, as a string
     */
    public String getReference();

    /**
     * @return true if the context is variant (i.e. contains a non-reference allele)
     */
    public boolean isVariant();

    /**
     * @return true if the context is strictly bi-allelic
     */
    public boolean isBiallelic();

    /**
     * @return true if the context represents point mutations only (i.e. no indels or structural variants)
     */
    public boolean isPointMutation();

    /**
     * @param sample  the sample name
     *
     * @return true if the context contains data for the given sample
     */
    public boolean hasSample(String sample);

    /**
     * @param sample  the sample name
     *
     * @return true if the context contains data for only a single instance of the given sample
     */
    public boolean hasSingleInstanceOfSample(String sample);

    /**
     * @return the number of samples in this context.
     *      If a given sample has more than one instance in this context, each one is counted separately
     */
    public int size();

    /**
     * @return the number of unique samples in this context
     */
    public int numSamples();

    /**
     * @return the set of all unique sample names in this context
     */
    public Set<String> getSampleNames();

    /**
     * @return true if the context represents variants with associated genotypes
     */
    public boolean hasGenotypes();

    /**
     * @return set of all Genotypes associated with this context
     */
    public Set<Genotype> getGenotypes();

    /**
     * @param sample  the sample name
     *
     * @return set of all Genotypes associated with the given sample in this context
     */
    public Set<Genotype> getGenotypes(String sample);

    /**
     * @param sample  the sample name
     *
     * @return the Genotype associated with the given sample in this context or null if the sample is not in this context
     *     throws an exception if multiple instances of the sample exist in the context
     */
    public Genotype getGenotype(String sample);    

    /**
     * @param sample  the sample name
     *
     * @return the Variation associated with the given sample in this context or null if the sample is not in this context
     */
    public Variation toVariation(String sample);

    /**
     * @return set of all subclasses within this context
     */
    public Set<String> getSubclasses();

    /**
     * @param subclass  the name of a subclass of variants to select
     *
     * @return a subset of this context which selects based on the given subclass
     */
    public VariantContext select(String subclass);

    /**
     * @param expr  a jexl expression describing how to filter this context
     *
     * @return a subset of this context which is filtered based on the given expression
     */
    public VariantContext filter(String expr);

    /**
     * @return a set of variant contexts where each one represents (all instances of) a single sample from this context
     */
    public Set<VariantContext> splitBySample();

    /**
     * get the alternate allele frequency of this variant context, if the context is variable.
     * If this context is not variable or a frequency cannot be provided, this method should return 0.
     *
     * WARNING: This method is valid only for bi-allelic contexts; the contract is to check isBiallelic()
     * before calling this method
     *
     * @return double the minor allele frequency
     */
    public double getNonRefAlleleFrequency();

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the postive number space log based error estimate
     */
    public double getNegLog10PError();

    /**
     * gets the alleles.  This method should return all of the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles. If the reference base is not an allele in this varation
     * it will not be in the list (i.e. there is no guarantee that the reference base is in the list).
     *
     * @return an alternate allele list
     */
    public List<String> getAlleles();

    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    public List<String> getAlternateAlleleList();

}