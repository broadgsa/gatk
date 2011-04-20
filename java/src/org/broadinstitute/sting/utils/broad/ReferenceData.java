/*
 * Copyright (c) 2011, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.broad;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Tracks data related to reference files at the Broad.
 */
public enum ReferenceData {
    /**
     * HG18 reference data
     */
    HG18("hg18"),

    /**
     * HG19 reference data
     */
    HG19("hg19");

    private static final String REFSEQ_DIR = "/humgen/gsa-hpprojects/GATK/data/Annotations/refseq/";
    private static final String DBSNP_DIR = "/humgen/gsa-hpprojects/GATK/data/";

    private final String name;
    private final String reference;
    private final String refseq;
    private final Map<Integer,String> dbsnps;

    ReferenceData(String name) {
        this.name = name;
        Map<Integer,String> dbsnps = new TreeMap<Integer,String>();
        if ("hg18".equals(name)) {
            this.reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
            this.refseq = REFSEQ_DIR + "refGene-big-table-hg18.txt";
            dbsnps.put(129, DBSNP_DIR + "dbsnp_129_hg18.rod");
        } else if ("hg19".equals(name)) {
            this.reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
            this.refseq = REFSEQ_DIR + "refGene-big-table-hg19.txt";
            dbsnps.put(129, DBSNP_DIR + "dbsnp_129_b37.vcf");
            dbsnps.put(132, DBSNP_DIR + "dbsnp_132_b37.vcf");
        } else
            throw new UnsupportedOperationException("Unknown reference: " + name);
        this.dbsnps = Collections.unmodifiableMap(dbsnps);
    }

    /**
     * Returns the name of the reference.
     * @return the name of the reference.
     */
    public String getName() {
        return name;
    }

    /**
     * Returns the path to the fasta.
     * @return the path to the fasta.
     */
    public String getReference() {
        return reference;
    }

    /**
     * Returns the path to the refseq table.
     * @return the path to the refseq table.
     */
    public String getRefseq() {
        return refseq;
    }

    /**
     * Returns the dbsnp versions available.
     * @return the dbsnp versions available.
     */
    public Set<Integer> getDbsnpVersions() {
        return dbsnps.keySet();
    }

    /**
     * Returns the dbsnp path for the version.
     * @param version version from getDbsnpVersions()
     * @return the dbsnp path for the version.
     */
    public String getDbsnp(int version) {
        return dbsnps.get(version);
    }

    /**
     * Returns the dbsnp type for the version, "VCF" or "DBSNP".
     * @param version version from getDbsnpVersions()
     * @return the dbsnp type for the version, "VCF" or "DBSNP".
     */
    public String getDbsnpType(int version) {
        String dbsnp = getDbsnp(version);
        if (dbsnp == null)
            return null;
        return dbsnp.toLowerCase().endsWith(".vcf") ? "VCF" : "DBSNP";
    }

    /**
     * Returns the reference data based on the path or null.
     * @param reference path to the reference
     * @return the reference data based on the path or null.
     */
    public static ReferenceData getByReference(String reference) {
        for (ReferenceData data: ReferenceData.values())
            if (data.reference.equals(reference))
                return data;
        return null;
    }
}
