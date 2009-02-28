/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.importer.genotype;

/**
 * data class for storing snp info
 *
 * @author Doug Voet
 */
public class SNP {
    private final byte referenceIndex;
    private final int position;
    private final byte allele1;
    private final byte allele2;

    public SNP(byte chromosome, int position, byte allele1, byte allele2) {
        this.referenceIndex = chromosome;
        this.position = position;
        this.allele1 = allele1;
        this.allele2 = allele2;
    }

    public byte getReferenceIndex() { return referenceIndex; }
    public int getPosition() { return position; }
    public byte getAllele1() { return allele1; }
    public byte getAllele2() { return allele2; }
}
