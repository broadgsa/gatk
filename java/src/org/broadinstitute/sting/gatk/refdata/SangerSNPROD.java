package org.broadinstitute.sting.gatk.refdata;

import java.util.*;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

public class SangerSNPROD extends TabularROD implements SNPCallFromGenotypes {
    public SangerSNPROD(final String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        loc = GenomeLocParser.createGenomeLoc(this.get("0"), Long.parseLong(this.get("1")));
        return loc;
    }
    public String getRefBasesFWD() { return this.get("2"); }
    public char getRefSnpFWD() throws IllegalStateException { return getRefBasesFWD().charAt(0); }
    public String getAltBasesFWD() { return this.get("3"); }
    public char getAltSnpFWD() throws IllegalStateException { return getAltBasesFWD().charAt(0); }
    public boolean isReference()   { return getVariationConfidence() < 0.01; }
    public boolean isSNP()         { return ! isReference(); }
    public boolean isInsertion()   { return false; }
    public boolean isDeletion()    { return false; }
    public boolean isIndel()       { return false; }
    public double getMAF()         { return -1; }
    public double getHeterozygosity() { return -1; }
    public boolean isGenotype()    { return false; }
    public double getVariationConfidence() { return -1; }
    public double getConsensusConfidence() { return -1; }
    public List<String> getGenotype() throws IllegalStateException { throw new IllegalStateException(); }
    public int getPloidy() throws IllegalStateException { return 2; }
    public boolean isBiallelic()   { return true; }
    public int length() { return 1; }

    // SNPCallFromGenotypes interface
    public int nIndividuals()      { return -1; }
    public int nHomRefGenotypes()  { return -1; }
    public int nHetGenotypes()     { return -1; }
    public int nHomVarGenotypes()  { return -1; }
    public List<Genotype> getGenotypes() { return null; }
}
