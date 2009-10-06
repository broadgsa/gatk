package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.genotype.Genotype;

import java.util.ArrayList;
import java.util.List;
/**
 * loc ref alt EM_alt_freq discovery_likelihood discovery_null discovery_prior discovery_lod EM_N n_ref n_het n_hom
 * chr1:1104840 A N 0.000000 -85.341265 -85.341265 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104841 C N 0.000000 -69.937928 -69.937928 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104842 A N 0.000000 -84.816002 -84.816002 0.000000 0.000000 324.000000 162 0 0
 *
 */
public class PooledEMSNPROD extends TabularROD implements SNPCallFromGenotypes, VariationRod {
    public PooledEMSNPROD(final String name) {
        super(name);
    }
    
    //GenomeLoc getLocation();
    public String getRefBasesFWD() { return this.get("ref"); }
    public char getRefSnpFWD() throws IllegalStateException { return getRefBasesFWD().charAt(0); }
    public String getAltBasesFWD() { return this.get("alt"); }
    public char getAltSnpFWD() throws IllegalStateException { return getAltBasesFWD().charAt(0); }
    public boolean isReference()   { return getVariationConfidence() < 0.01; }

    /**
     * gets the alternate base.  Use this method if we're biallelic
     *
     * @return
     */
    @Override
    public String getAlternateBases() {
        return getAltBasesFWD();
    }

    /**
     * gets the alternate bases.  Use this method if the allele count is greater then 2
     *
     * @return
     */
    @Override
    public List<String> getAlternateBaseList() {
        List<String> str = new ArrayList<String>();
        str.add(this.getAltBasesFWD());
        return str;
    }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    @Override
    public double getNonRefAlleleFrequency() {
        return this.getMAF();
    }

    /** @return the VARIANT_TYPE of the current variant */
    @Override
    public VARIANT_TYPE getType() {
        if (isSNP()) {
            return VARIANT_TYPE.SNP;
        }
        return VARIANT_TYPE.REFERENCE;
    }

    public boolean isSNP()         { return ! isReference(); }
    public boolean isInsertion()   { return false; }
    public boolean isDeletion()    { return false; }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    @Override
    public String getReference() {
        return this.get("ref");
    }

    public boolean isIndel()       { return false; }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     * of
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getAlternativeBaseForSNP() {
        return this.getAltSnpFWD();
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        return this.getRefSnpFWD();
    }

    public double getMAF()         { return Double.parseDouble(this.get("EM_alt_freq")); }
    public double getHeterozygosity() { return 2 * getMAF() * (1 - getMAF()); }
    public boolean isGenotype()    { return false; }
    public double getVariationConfidence() { return Double.parseDouble(this.get("lod")); }
    public double getConsensusConfidence() { return -1; }
    public List<String> getGenotype() throws IllegalStateException { throw new IllegalStateException(); }
    public int getPloidy() throws IllegalStateException { return 2; }
    public boolean isBiallelic()   { return true; }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        return this.getVariationConfidence();
    }

    public int length() { return 1; }

    // SNPCallFromGenotypes interface
    public int nIndividuals()      { return Integer.parseInt(this.get("EM_N")); }
    public int nHomRefGenotypes()  { return Integer.parseInt(this.get("n_ref")); }
    public int nHetGenotypes()     { return Integer.parseInt(this.get("n_het")); }
    public int nHomVarGenotypes()  { return Integer.parseInt(this.get("n_hom")); }
    public List<Genotype> getGenotypes() { return null; }
}
