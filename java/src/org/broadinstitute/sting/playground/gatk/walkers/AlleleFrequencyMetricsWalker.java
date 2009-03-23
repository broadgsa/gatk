package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.BasicLociWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Mar 18, 2009
 * Time: 5:28:58 PM
 * To change this template use File | Settings | File Templates.
 */

public class AlleleFrequencyMetricsWalker extends BasicLociWalker<AlleleFrequencyEstimate, String> 
{

    long dbsnp_hits=0;
    long num_variants=0;
    long num_loci_total=0;
    long num_loci_confident=0;
    double LOD_cutoff = 5;

    AlleleFrequencyWalker caller;

    public AlleleFrequencyEstimate map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) 
    {
        AlleleFrequencyEstimate alleleFreq = caller.map(rodData, ref, context);

        boolean is_dbSNP_SNP = false;

        for ( ReferenceOrderedDatum datum : rodData ) 
        {
            if ( datum != null && datum instanceof rodDbSNP) 
            {
                rodDbSNP dbsnp = (rodDbSNP)datum;
                if (dbsnp.isSNP()) is_dbSNP_SNP = true;
            }
        }

        num_loci_total += 1;

        if (Math.abs(alleleFreq.LOD) >= LOD_cutoff) { num_loci_confident += 1; }

        if (alleleFreq.getQstar() > 0.0 && alleleFreq.getLOD() >= LOD_cutoff)
        { 
            // Confident variant.
           
            num_variants += 1;

            if (is_dbSNP_SNP) 
            {
                dbsnp_hits += 1;
            }
        }

        return alleleFreq;
    }

    public void printMetrics() 
    {
        if (num_loci_total == 0) { return; }

        System.out.printf("\n");
        System.out.printf("METRICS Allele Frequency Metrics (LOD >= %.0f)\n", LOD_cutoff);
        System.out.printf("METRICS -------------------------------------------------\n");
        System.out.printf("METRICS Total loci                         : %d\n", num_loci_total);
        System.out.printf("METRICS Total called with confidence       : %d (%.2f%%)\n", num_loci_confident, 100.0 * (float)num_loci_confident / (float)num_loci_total);
        if (num_variants != 0)
        {
	        System.out.printf("METRICS Number of Variants                 : %d (%.2f%%) (1/%d)\n", num_variants, 100.0 * (float)num_variants / (float)num_loci_confident, num_loci_confident / num_variants);
	        System.out.printf("METRICS Fraction of variant sites in dbSNP : %.2f%%\n", 100.0 * (float)dbsnp_hits / (float)num_variants);
        }
        System.out.println();
    }

    public void onTraversalDone() 
    {
        printMetrics();
    }

    public String reduceInit() 
    { 
        caller = new AlleleFrequencyWalker();
        return ""; 
    }

    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
    {
        if ((alleleFreq.LOD >= 5) || (alleleFreq.LOD <= -5))
        {
	        System.out.print(String.format("RESULT %s %c %c %f %f %f %d\n", 
	                                        alleleFreq.location,
	                                        alleleFreq.ref, 
	                                        alleleFreq.alt, 
	                                        alleleFreq.qhat, 
	                                        alleleFreq.qstar, 
	                                        alleleFreq.LOD, 
	                                        alleleFreq.depth));
        }

        if (this.num_loci_total % 10000 == 0) { printMetrics(); }

        return "null";
    }

    
}
