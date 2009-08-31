package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;


public class CleanedOutSNPROD extends TabularROD {

    private static final String REAL_SNP_STRING = "SAME_SNP";
    private static final String FALSE_SNP_STRING = "NOT_SNP";

    public CleanedOutSNPROD(String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        return GenomeLocParser.parseGenomeLoc(this.get("0"));
    }

    public boolean isRealSNP() { return this.get("1").equals(REAL_SNP_STRING); }
 }