package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.AlleleFrequencyWalker;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Mar 18, 2009
 * Time: 5:28:58 PM
 * To change this template use File | Settings | File Templates.
 */


public class AlleleFrequencyMetrics {

    long dbsnp_tp=0;
    long dbsnp_fp=0;
    long num_snps=0;
    long num_loci=0;

    public void calculateMetrics(List<ReferenceOrderedDatum> rodData, AlleleFrequencyWalker.AlleleFrequencyEstimate alleleFreq) {

        boolean is_dbSNP_SNP = false;

        for ( ReferenceOrderedDatum datum : rodData ) {
            if ( datum != null && datum instanceof rodDbSNP) {
                rodDbSNP dbsnp = (rodDbSNP)datum;
                if (dbsnp.isSNP()) is_dbSNP_SNP = true;
            }
        }

        if (alleleFreq.getQstar() > 0.0 && alleleFreq.getLogOddsVarRef() >= 5.0) { // we confidently called it a SNP!
            if (is_dbSNP_SNP) {
                dbsnp_tp += 1;
            }else{
                dbsnp_fp += 1;
            }
        }

        if (alleleFreq.getQstar() > 0.0) num_snps++;
        num_loci++;

    }

    public void printMetrics() {
        System.out.println("\nAllele Frequency Metrics:\n");
        System.out.printf("Precision of LOD >= 5 SNPs w.r.t dbSNP: %.2f\n", (float)dbsnp_tp / (dbsnp_fp + dbsnp_tp) * 100);
        System.out.printf("\\--TP: %d\n", dbsnp_tp);
        System.out.printf("\\--FP: %d\n", dbsnp_fp);
        System.out.println();
        System.out.printf("SNPs (LOD > 0): %d\n", num_snps);
        System.out.printf("Total loci: %d\n", num_loci);
        System.out.printf("SNPs / loci: 1/%.0f\n", (float)num_loci/num_snps);
        System.out.println();
        



    }
    
}
