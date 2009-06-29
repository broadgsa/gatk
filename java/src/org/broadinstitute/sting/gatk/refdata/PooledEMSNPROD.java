package org.broadinstitute.sting.gatk.refdata;

import java.util.*;
import java.util.regex.MatchResult;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.xReadLines;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.apache.log4j.Logger;

/**
 * loc ref alt EM_alt_freq discovery_likelihood discovery_null discovery_prior discovery_lod EM_N n_ref n_het n_hom
 * chr1:1104840 A N 0.000000 -85.341265 -85.341265 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104841 C N 0.000000 -69.937928 -69.937928 0.000000 0.000000 324.000000 162 0 0
 * chr1:1104842 A N 0.000000 -84.816002 -84.816002 0.000000 0.000000 324.000000 162 0 0
 *
 */
public class PooledEMSNPROD extends TabularROD implements SNPCallFromGenotypes {
    public PooledEMSNPROD(final String name) {
        super(name);
    }
    
    //GenomeLoc getLocation();
    public String getRefBasesFWD() { return this.get("ref"); }
    public char getRefSnpFWD() throws IllegalStateException { return getRefBasesFWD().charAt(0); }
    public String getAltBasesFWD() { return this.get("alt"); }
    public char getAltSnpFWD() throws IllegalStateException { return getAltBasesFWD().charAt(0); }
    public boolean isReference()   { return getVariationConfidence() < 0.01; }
    public boolean isSNP()         { return ! isReference(); }
    public boolean isInsertion()   { return false; }
    public boolean isDeletion()    { return false; }
    public boolean isIndel()       { return false; }
    public double getMAF()         { return Double.parseDouble(this.get("EM_alt_freq")); }
    public double getHeterozygosity() { return 2 * getMAF() * (1 - getMAF()); }
    public boolean isGenotype()    { return false; }
    public double getVariationConfidence() { return Double.parseDouble(this.get("discovery_lod")); }
    public double getConsensusConfidence() { return -1; }
    public List<String> getGenotype() throws IllegalStateException { throw new IllegalStateException(); }
    public int getPloidy() throws IllegalStateException { return 2; }
    public boolean isBiallelic()   { return true; }

    // SNPCallFromGenotypes interface
    public int nIndividuals()      { return Integer.parseInt(this.get("EM_N")); }
    public int nHomRefGenotypes()  { return Integer.parseInt(this.get("n_ref")); }
    public int nHetGenotypes()     { return Integer.parseInt(this.get("n_het")); }
    public int nHomVarGenotypes()  { return Integer.parseInt(this.get("n_hom")); }
    public List<Genotype> getGenotypes() { return null; }
}
