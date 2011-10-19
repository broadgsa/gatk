package org.broadinstitute.sting.utils.variantcontext;

import java.util.*;

/**
 * This class emcompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
public class MutableGenotype extends Genotype {
    public MutableGenotype(Genotype parent) {
        super(parent.getSampleName(), parent.getAlleles(), parent.getNegLog10PError(), parent.getFilters(), parent.getAttributes(), parent.isPhased());
    }

    public MutableGenotype(String sampleName, Genotype parent) {
        super(sampleName, parent.getAlleles(), parent.getNegLog10PError(), parent.getFilters(), parent.getAttributes(), parent.isPhased());
    }


    public MutableGenotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean genotypesArePhased) {
        super(sampleName, alleles, negLog10PError, filters, attributes, genotypesArePhased);
    }

    public MutableGenotype(String sampleName, List<Allele> alleles, double negLog10PError) {
        super(sampleName, alleles, negLog10PError);
    }

    public MutableGenotype(String sampleName, List<Allele> alleles) {
        super(sampleName, alleles);
    }

    public Genotype unmodifiableGenotype() {
        return new Genotype(getSampleName(), getAlleles(), getNegLog10PError(), getFilters(), getAttributes(), isPhased());
    }


    /**
     *
     * @param alleles list of alleles
     */
    public void setAlleles(List<Allele> alleles) {
        this.alleles = new ArrayList<Allele>(alleles);
        validate();
    }

    public void setPhase(boolean isPhased) {
        super.isPhased = isPhased;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // InferredGeneticContext mutation operators
    //
    // ---------------------------------------------------------------------------------------------------------
    public void setName(String name)                    { commonInfo.setName(name); }
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

}