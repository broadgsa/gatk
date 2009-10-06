package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser; 
import org.broadinstitute.sting.utils.genotype.Genotype;
import java.util.Arrays;
import java.util.List;

public class SangerSNPROD extends TabularROD implements SNPCallFromGenotypes {
    public SangerSNPROD(final String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        loc = GenomeLocParser.createGenomeLoc(this.get("0"), Long.parseLong(this.get("1")));
        return loc;
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    @Override
    public String getReference() {
         return String.valueOf(getRefBasesFWD().charAt(0));
    }

    public String getRefBasesFWD() { return this.get("2"); }
    public char getRefSnpFWD() throws IllegalStateException { return getRefBasesFWD().charAt(0); }
    public String getAltBasesFWD() { return this.get("3"); }
    public char getAltSnpFWD() throws IllegalStateException { return getAltBasesFWD().charAt(0); }
    public boolean isReference()   { return getVariationConfidence() < 0.01; }

    /**
     * get the frequency of this variant, if we're a variant.  If we're reference this method
     * should return 0.
     *
     * @return double with the stored frequency
     */
    @Override
    public double getNonRefAlleleFrequency() {
       return this.getMAF();
    }

    /**
     * A convenience method, for switching over the variation type
     *
     * @return the VARIANT_TYPE of the current variant
     */
    @Override
    public VARIANT_TYPE getType() {
        if (isReference()) return VARIANT_TYPE.REFERENCE;
        else return VARIANT_TYPE.SNP;
    }

    public boolean isSNP()         { return ! isReference(); }
    public boolean isInsertion()   { return false; }
    public boolean isDeletion()    { return false; }
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

    public double getMAF()         { return -1; }
    public double getHeterozygosity() { return -1; }
    public boolean isGenotype()    { return false; }
    public double getVariationConfidence() { return -1; }
    public double getConsensusConfidence() { return -1; }
    public List<String> getGenotype() throws IllegalStateException { throw new IllegalStateException(); }
    public int getPloidy() throws IllegalStateException { return 2; }
    public boolean isBiallelic()   { return true; }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the postive number space log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        return this.getVariationConfidence();
    }

    /**
     * gets the alternate base.  Use this method if we're biallelic
     *
     * @return
     */
    @Override
    public String getAlternateBases() {
        return this.get("3");
    }

    /**
     * gets the alternate bases.  Use this method if the allele count is greater then 2 (not biallelic)
     *
     * @return
     */
    @Override
    public List<String> getAlternateBaseList() {
        return Arrays.asList(this.get("3"));
    }

    public int length() { return 1; }

    // SNPCallFromGenotypes interface
    public int nIndividuals()      { return -1; }
    public int nHomRefGenotypes()  { return -1; }
    public int nHetGenotypes()     { return -1; }
    public int nHomVarGenotypes()  { return -1; }
    public List<Genotype> getGenotypes() { return null; }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getLocation().getContig() + "\t" + getLocation().getStart() + "\t");
        sb.append(getRefBasesFWD() + "\t" + getAltBasesFWD());
        return sb.toString();
    }
}
