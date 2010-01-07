package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Arrays;
import java.util.List;

public class SimpleIndelROD extends TabularROD implements VariationRod {

    private boolean KGENOMES_FORMAT = false, checkedFormat = false;

    public SimpleIndelROD(String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        long pos = Long.parseLong(this.get("1"));
        return GenomeLocParser.createGenomeLoc(this.get("0"), pos, (isDeletion() ? pos+length() : pos+1));
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    public String getReference() {
        return String.valueOf(getRef());
    }

    public List<String> getFWDAlleles() {
        if ( is1KGFormat() )
            return Arrays.asList(this.get("4"));

        String str = this.get("3");
        return Arrays.asList(str.substring(1, str.indexOf(":")));    
    }

    public String getFWDRefBases() { return ""; }
    public String getAltBasesFWD() { return getFWDAlleles().get(0); }
    public String getRefBasesFWD() { return ""; }
    public char getRefSnpFWD() { throw new IllegalStateException("I'm an indel, not a SNP"); }
    public char getAltSnpFWD() { throw new IllegalStateException("I'm an indel, not a SNP"); }
    public char getRef() { return 'N'; }
    public List<String> getGenotype() { return getFWDAlleles(); }
    public boolean isGenotype() { return false; }
    public boolean isPointGenotype() { return false; }
    public boolean isIndelGenotype() { return true; }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    public double getNonRefAlleleFrequency() {
        return 0.0;
    }

    /** @return the VARIANT_TYPE of the current variant */
    public VARIANT_TYPE getType() {
        return isInsertion() ? VARIANT_TYPE.INSERTION : VARIANT_TYPE.DELETION;
    }

    public boolean isSNP() { return false; }
    public boolean isReference() { return false; }

    public boolean isInsertion() {
        if ( is1KGFormat() )
            return this.get("3").equals("I");
        return this.get("3").charAt(0) == '+';
    }
    public boolean isDeletion() {
        if ( is1KGFormat() )
            return this.get("3").equals("D");
        return this.get("3").charAt(0) == '-';
    }
    public boolean isIndel() { return true; }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     * of
     *
     * @return a char, representing the alternate base
     */
    public char getAlternativeBaseForSNP() {
        return getAltSnpFWD();
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    public char getReferenceForSNP() {
        return getRefSnpFWD();
    }

    public double getVariantConfidence() { return 0.0; }
    public double getVariationConfidence() { return 0.0; }
    public double getConsensusConfidence() { return 0.0; }
    public boolean isBiallelic() { return true; }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    public double getNegLog10PError() {
        return getVariationConfidence();
    }

    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    public List<String> getAlternateAlleleList() {
        List<String> ret = getAlleleList();
        for (String val : ret) {
            if (val.equals(this.getReference())) ret.remove(val);
        }
        return ret;
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles.
     *
     * @return an alternate allele list
     */
    public List<String> getAlleleList() {
        return this.getFWDAlleles();
    }

    public boolean isHom() { return false; }
    public boolean isHet() { return false; }
    public double getHeterozygosity() { return 0.0; }
    public double getMAF() { return 0.0; }
    public int getPloidy() { return 2; }
    public int length() {
        if ( is1KGFormat() )
             return Integer.parseInt(this.get("2"));
        return getFWDAlleles().get(0).length();
    }

    public boolean allowIncompleteRecords() {
        return true;
    }

    public String getSamplesString() {
        return (is1KGFormat() && this.get("5") != null ? this.get("5") : "");
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getLocation().getContig() + "\t" + getLocation().getStart() + "\t");
        sb.append(length() + "\t" + (isInsertion() ? "I" : "D") + "\t" + getFWDAlleles().get(0));
        String samples = getSamplesString();
        if ( samples.length() > 0 )
            sb.append("\t" + samples);
        return sb.toString();
    }

    private boolean is1KGFormat() {
        if ( !checkedFormat ) {
            checkedFormat = true;
            KGENOMES_FORMAT = this.get("3").equals("D") || this.get("3").equals("I");
        }
        return KGENOMES_FORMAT;
    }
}