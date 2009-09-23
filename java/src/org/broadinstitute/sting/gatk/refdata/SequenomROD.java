package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.*;

public class SequenomROD extends TabularROD implements AllelicVariant {

    public SequenomROD(String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        long pos = Long.parseLong(this.get("1"));
        return GenomeLocParser.createGenomeLoc(this.get("0"), pos);
    }

    public List<String> getFWDAlleles() {
        return Arrays.asList(this.get("2"));
    }

    public String getFWDRefBases() { return ""; }
    public String getAltBasesFWD() { return getFWDAlleles().get(0); }
    public String getRefBasesFWD() { return ""; }
    public char getRefSnpFWD() { return 'N'; }
    public char getAltSnpFWD() { return getAltBasesFWD().charAt(0); }
    public char getRef() { return 'N'; }
    public List<String> getGenotype() { return getFWDAlleles(); }
    public boolean isGenotype() { return false; }
    public boolean isPointGenotype() { return true; }
    public boolean isIndelGenotype() { return false; }
    public boolean isSNP() { return true; }
    public boolean isReference() { return false; }
    public boolean isInsertion() { return false; }
    public boolean isDeletion() { return false; }
    public boolean isIndel() { return false; }
    public double getVariantConfidence() { return 0.0; }
    public double getVariationConfidence() { return 0.0; }
    public double getConsensusConfidence() { return 0.0; }
    public boolean isBiallelic() { return true; }
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