package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.*;

public class SimpleIndelROD extends TabularROD implements Genotype {

    public SimpleIndelROD(String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        return GenomeLocParser.createGenomeLoc(this.get("0"), Long.parseLong(this.get("1")));
    }

    public List<String> getFWDAlleles() {
        String str = this.get("3");
        return Arrays.asList(str.substring(0, str.indexOf(":")));    
    }

    public String getFWDRefBases() { return ""; }
    public char getRef() { return 'N'; }
    public boolean isPointGenotype() { return false; }
    public boolean isIndelGenotype() { return true; }
    public boolean isSNP() { return false; }
    public boolean isReference() { return false; }
    public boolean isInsertion() { return this.get("3").charAt(0) == '+'; }
    public boolean isDeletion() { return this.get("3").charAt(0) == '-'; }
    public boolean isIndel() { return false; }
    public double getVariantConfidence() { return 0.0; }
    public double getConsensusConfidence() { return 0.0; }
    public boolean isBiallelic() { return true; }
    public boolean isHom() { return false; }
    public boolean isHet() { return false; }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getLocation().getContig() + "\t" + getLocation().getStart() + "\t");
        String indel = getFWDAlleles().get(0);
        sb.append((indel.length()-1) + "\t" + (isInsertion() ? "I" : "D") + "\t" + indel.substring(1));
        return sb.toString();
    }
 }