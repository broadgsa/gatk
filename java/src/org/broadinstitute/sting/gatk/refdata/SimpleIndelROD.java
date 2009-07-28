package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.*;

public class SimpleIndelROD extends TabularROD implements Genotype, AllelicVariant {

    private boolean KGENOMES_FORMAT = false, checkedFormat = false;

    public SimpleIndelROD(String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        return GenomeLocParser.createGenomeLoc(this.get("0"), Long.parseLong(this.get("1")));
    }

    public List<String> getFWDAlleles() {
        if ( is1KGFormat() )
            return Arrays.asList(this.get("4"));

        String str = this.get("3");
        return Arrays.asList(str.substring(0, str.indexOf(":")));    
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
    public double getVariantConfidence() { return 0.0; }
    public double getVariationConfidence() { return 0.0; }
    public double getConsensusConfidence() { return 0.0; }
    public boolean isBiallelic() { return true; }
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

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getLocation().getContig() + "\t" + getLocation().getStart() + "\t");
        String indel = getFWDAlleles().get(0);
        sb.append((indel.length()-1) + "\t" + (isInsertion() ? "I" : "D") + "\t" + indel.substring(1));
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