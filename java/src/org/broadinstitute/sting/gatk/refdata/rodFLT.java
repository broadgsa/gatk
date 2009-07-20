package org.broadinstitute.sting.gatk.refdata;

import java.util.*;

import org.broadinstitute.sting.utils.*;

/**
 * 1	469	C	S	51	52	0.00	44	1	C	0	M
 * 1	492	C	Y	130	51	0.00	46	1	C	0	T
 * 1	20786	G	K	5	6	0.00	46	1	T	38	G
 *
 */
public class rodFLT extends TabularROD implements SNPCallFromGenotypes {
    public rodFLT(final String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        loc = GenomeLocParser.createGenomeLoc(this.get("0"), Long.parseLong(this.get("1")));
        return loc;
    }
    public String getRefBasesFWD() { return this.get("2"); }
    public char getRefSnpFWD() throws IllegalStateException { return getRefBasesFWD().charAt(0); }
    public String getAltBasesFWD() { return new String(BaseUtils.iupacToBases(this.get("3").charAt(0))); }
    public char getAltSnpFWD() throws IllegalStateException {
        char[] bases = BaseUtils.iupacToBases(this.get("3").charAt(0));
        if ( bases[0] != getRefSnpFWD() )
            return bases[0];
        else
            return bases[1];
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(loc.getContig() + "\t" + loc.getStart() + "\t");
        sb.append(getRefSnpFWD() + "\t-1\t-1\t" + getAltBasesFWD());
        for (int i=0; i < 12; i++)
            sb.append("\t0");
        return sb.toString();
    }

    public boolean isReference()   { return false; }
    public boolean isSNP()         { return true; }
    public boolean isInsertion()   { return false; }
    public boolean isDeletion()    { return false; }
    public boolean isIndel()       { return false; }
    public double getMAF()         { return 0.0; }
    public double getHeterozygosity() { return 0.0; }
    public boolean isGenotype()    { return false; }
    public double getVariationConfidence() {
        double value = Double.parseDouble(this.get("6"));
        return (value > 0.0 ? value : 1.0);
    }
    public double getConsensusConfidence() { return -1; }
    public List<String> getGenotype() throws IllegalStateException { throw new IllegalStateException(); }
    public int getPloidy() throws IllegalStateException { return 2; }
    public boolean isBiallelic()   { return true; }

    // SNPCallFromGenotypes interface
    public int nIndividuals()      { return -1; }
    public int nHomRefGenotypes()  { return -1; }
    public int nHetGenotypes()     { return -1; }
    public int nHomVarGenotypes()  { return -1; }
    public List<Genotype> getGenotypes() { return null; }
}