package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.Arrays;
import java.util.List;

public class SequenomROD extends TabularROD implements Variation {

    public SequenomROD(String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        long pos = Long.parseLong(this.get("1"));
        return GenomeLocParser.createGenomeLoc(this.get("0"), pos);
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    @Override
    public String getReference() {
        return "";
    }

    public List<String> getFWDAlleles() {
        return Arrays.asList(this.get("2"));
    }

    public String getAltBasesFWD() { return getFWDAlleles().get(0); }

    /**
     * get the frequency of this variant, if we're a variant.  If we're reference this method
     * should return 0.
     *
     * @return double with the stored frequency
     */
    @Override
    public double getNonRefAlleleFrequency() {
        return 0.0;
    }

    /**
     * A convenience method, for switching over the variation type
     *
     * @return the VARIANT_TYPE of the current variant
     */
    @Override
    public VARIANT_TYPE getType() {
        return VARIANT_TYPE.SNP;
    }

    public boolean isSNP() { return true; }
    public boolean isReference() { return false; }
    public boolean isInsertion() { return false; }
    public boolean isDeletion() { return false; }
    public boolean isIndel() { return false; }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     * of
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getAlternativeBaseForSNP() {
        return getAltBasesFWD().charAt(0);
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        return 'N';
    }

    public boolean isBiallelic() { return true; }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the postive number space log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        return 0.0;
    }

    /**
     * gets the alternate bases.  Use this method if we're biallelic
     *
     * @return
     */
    @Override
    public String getAlternateBases() {
        return getAltBasesFWD();
    }

    /**
     * gets the alternate bases.  Use this method if the allele count is greater then 2 (not biallelic)
     *
     * @return
     */
    @Override
    public List<String> getAlternateBaseList() {
        throw new StingException("SequenomRod is not biallelic");
    }

    public boolean isHom() { return false; }
    public boolean isHet() { return false; }
    public double getHeterozygosity() { return 0.0; }
    public double getMAF() { return 0.0; }
    public int getPloidy() { return 2; }
    public int length() { return 1; }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getLocation().getContig() + "\t" + getLocation().getStart() + "\t" + getFWDAlleles().get(0));
        return sb.toString();
    }
 }
