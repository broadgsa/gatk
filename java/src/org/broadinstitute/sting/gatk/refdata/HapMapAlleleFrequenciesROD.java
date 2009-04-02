package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;
import java.util.Arrays;

import edu.mit.broad.picard.util.SequenceUtil;

/**
 * ReferenceOrderedDatum class to hold HapMap AlleleFrequency Data
 */
public class HapMapAlleleFrequenciesROD extends ReferenceOrderedDatum {
    public GenomeLoc loc;       // genome location of SNP
                                // Reference sequence chromosome or scaffold
                                // Start and stop positions in chrom

    public String rsNumber;     // dbsnp rsNumber for this site

    public String hgBuild;

    public char Strand;         // strand of the supplied alleles

    public char refAllele;
    public char varAllele;

    public double refFreq;
    public double varFreq;


    public String strand;   // maybe we don't need these?
    public String alleles;  // maybe we don't need these?
    public Integer refCounts; // maybe we don't need these?
    public Integer varCounts; // maybe we don't need these?
    public Integer totalCounts; // maybe we don't need these?


    public GenomeLoc getLocation() { return loc; }

    public String toString() {
        //rs11511647      HG18    chr10   62765   +       T/C     T       C       21      97      0.178   0.822   118

        return String.format(
                "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%1.3f\t%1.3f\t%d",
                rsNumber, hgBuild, getContig(), getStart(), strand, alleles, refAllele, varAllele,
                refCounts, varCounts, refFreq, varFreq, totalCounts);
    }

    public String toSimpleString() {
        return String.format("%s:%s:%s:%1.3f", rsNumber, alleles, strand, varFreq);
    }

    public String repl() {
        return toString();
    }

    public void parseLine(final String[] parts) {
        try {
            // rs11511647 <=> HG18 <=> chr10 <=> 62765 <=> + <=> T/C <=> T <=> C <=> 21 <=> 97 <=> 0.178 <=> 0.822 <=> 118

            rsNumber = parts[0];    //rs#
            hgBuild = parts[1];     // build

            String contig = parts[2];   // chrom
            long start = Long.parseLong(parts[3]); // The final is 1 based
            long stop = start;

            strand = parts[4];                   // strand
            alleles = parts[5];                       //alleles
            refAllele = parts[6].charAt(0); // ref_allele
            varAllele = parts[7].charAt(0); // var_allele
            refCounts = Integer.parseInt(parts[8]); // CEU_ref
            varCounts = Integer.parseInt(parts[9]); // CEU_var
            refFreq = Double.parseDouble(parts[10]); // CEU_ref_freq
            varFreq = Double.parseDouble(parts[11]); // CEU_var_freq
            totalCounts = Integer.parseInt(parts[12]); // CEU_var

            loc = new GenomeLoc(contig, start, stop);

        } catch ( RuntimeException e ) {
            System.out.printf("  Exception caught during parsing HapMap Allele Freq %s%n", Utils.join(" <=> ", parts));
            throw e;
        }
    }

    public double getVarAlleleFreq() { return this.varFreq; }

    public List<String> getAllelesFWD() {
        List<String> alleleList;
        if ( onFwdStrand() )
            alleleList = Arrays.asList(alleles.split("/"));
        else
            alleleList = Arrays.asList(SequenceUtil.reverseComplement(alleles).split("/"));

        return alleleList;
    }

    public boolean onFwdStrand() {
        return strand.equals("+");
    }



}
