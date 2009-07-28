package org.broadinstitute.sting.gatk.refdata;

import java.util.*;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

/**
 * loc ref alt EM_alt_freq discovery_likelihood discovery_null discovery_prior discovery_lod EM_N n_ref n_het n_hom
 * chr1:1104840 A N 0.000000 -85.341265 -85.341265 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104841 C N 0.000000 -69.937928 -69.937928 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104842 A N 0.000000 -84.816002 -84.816002 0.000000 0.000000 324.000000 162 0 0
 *
 */
public class KGenomesSNPROD extends TabularROD implements SNPCallFromGenotypes {
    public KGenomesSNPROD(final String name) {
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
    public double getMAF()         { return Double.parseDouble(this.get("4")); }
    public double getHeterozygosity() { return 2 * getMAF() * (1 - getMAF()); }
    public boolean isGenotype()    { return false; }
    public double getVariationConfidence() { return Double.parseDouble(this.get("8")); }
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
    public int length() { return 1; }
}