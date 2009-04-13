package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;

import java.util.List;
import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Apr 1, 2009
 * Time: 5:53:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class AlleleMetrics {

    long dbsnp_hits=0;
    long num_variants=0;
    long num_loci_total=0;
    long num_loci_confident=0;
    double LOD_cutoff = 5;
    long hapmap_genotype_correct = 0;
    long hapmap_genotype_incorrect = 0;
    long hapmap_refvar_correct = 0;
    long hapmap_refvar_incorrect = 0;

    protected PrintStream out;

    public AlleleMetrics(String MetricsOutputFile) {
        try
        {
            /*if ( MetricsOutputFile.equals("-") )
                this.out = out;
            else*/
            this.out = new PrintStream(MetricsOutputFile);
        }
        catch (Exception e)
        {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    public void nextPosition(AlleleFrequencyEstimate alleleFreq, RefMetaDataTracker tracker) {
        num_loci_total += 1;

        boolean is_dbSNP_SNP = false;
        boolean has_hapmap_chip_genotype = false;
        rodGFF hapmap_chip_genotype = null;

        for ( ReferenceOrderedDatum datum : tracker.getAllRods() )
        {
            if ( datum != null )
            {
                if ( datum instanceof rodDbSNP)
                {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    if (dbsnp.isSNP()) is_dbSNP_SNP = true;
                }

                if ( datum instanceof rodGFF )
                {
                    has_hapmap_chip_genotype = true;
                    hapmap_chip_genotype = (rodGFF)datum;
                }
            }
        }

        if (Math.abs(alleleFreq.lodVsRef) >= LOD_cutoff) { num_loci_confident += 1; }

        if (alleleFreq.qstar > 0.0 && alleleFreq.lodVsRef >= LOD_cutoff)
        {
            // Confident variant.

            num_variants += 1;

            if (is_dbSNP_SNP)
            {
                dbsnp_hits += 1;
            }
        }

        if (has_hapmap_chip_genotype) {

            // convert hapmap call to mixture of ref/nonref
            String hapmap_genotype = hapmap_chip_genotype.getFeature();
            long refs=0, alts=0;
            double hapmap_q;

            for (char c : hapmap_genotype.toCharArray()) {
                if (c == alleleFreq.ref) { refs++; }
                if (c == alleleFreq.alt) { alts++; }
            }

            if (refs+alts > 0) {
                hapmap_q = (float)alts / (refs+alts);
            }else{
                hapmap_q = -1;
            }

            // Hapmap debug info
            //out.format("HAPMAP DEBUG %.2f %.2f %.2f ", hapmap_q, alleleFreq.qstar, alleleFreq.lodVsRef);
            //String called_genotype = alleleFreq.asString();
            //out.format("%s %s %c %c", hapmap_genotype, called_genotype, alleleFreq.ref, alleleFreq.alt);

            if (alleleFreq.lodVsNextBest >= LOD_cutoff) {

                // Calculate genotyping performance - did we get the correct genotype of the N+1 choices?
                if (hapmap_q != -1 && hapmap_q == alleleFreq.qstar) {
                    hapmap_genotype_correct++;
                }else{
                    hapmap_genotype_incorrect++;
                    //System.out.printf(" INCORRECT GENOTYPE    Bases: %s", AlleleFrequencyWalker.getBases(context));
                    //out.printf(" INCORRECT GENOTYPE");
                    //AlleleFrequencyWalker.print_base_qual_matrix(AlleleFrequencyWalker.getOneBaseQuals(context));
                }
            }

            if (alleleFreq.lodVsRef >= LOD_cutoff || -1 * alleleFreq.lodVsRef >= LOD_cutoff) {

                // Now calculate ref / var performance - did we correctly classify the site as
                // reference or variant without regard to genotype; i.e. het/hom "miscalls" don't matter here
                boolean hapmap_var = hapmap_q != 0.0;
                boolean called_var = alleleFreq.qstar != 0.0;
                if (hapmap_q != -1 && hapmap_var != called_var) {
                    hapmap_refvar_incorrect++;
                    //out.printf(" INCORRECT REFVAR CALL");
                }else{
                    hapmap_refvar_correct++;
                }
            }

            //out.print("\n");
        }
    }

    public void printMetrics()
    {
        if (num_loci_total == 0) { return; }

        out.printf("\n");
        out.printf("Allele Frequency Metrics (LOD >= %.0f)\n", LOD_cutoff);
        out.printf("-------------------------------------------------\n");
        out.printf("Total loci                            : %d\n", num_loci_total);
        out.printf("Total called with confidence          : %d (%.2f%%)\n", num_loci_confident, 100.0 * (float)num_loci_confident / (float)num_loci_total);
        if (num_variants != 0)
        {
	        out.printf("Number of variants                    : %d (%.2f%%) (1/%d)\n", num_variants, 100.0 * (float)num_variants / (float)num_loci_confident, num_loci_confident / num_variants);
            out.printf("Fraction of variant sites in dbSNP    : %.2f%%\n", 100.0 * (float)dbsnp_hits / (float)num_variants);
            out.printf("-------------------------------------------------\n");
            out.printf("-- Hapmap Genotyping performance --\n");
            out.printf("Num. conf. calls at Hapmap chip sites : %d\n", hapmap_genotype_correct + hapmap_genotype_incorrect);
            out.printf("Conf. calls at chip sites correct     : %d\n", hapmap_genotype_correct);
            out.printf("Conf. calls at chip sites incorrect   : %d\n", hapmap_genotype_incorrect);
            out.printf("%% of confident calls that are correct : %.2f%%\n", 100.0 * (float) hapmap_genotype_correct / (float)(hapmap_genotype_correct + hapmap_genotype_incorrect));
            out.printf("-------------------------------------------------\n");
            out.printf("-- Hapmap Reference/Variant performance --\n");
            out.printf("Num. conf. calls at Hapmap chip sites : %d\n", hapmap_refvar_correct + hapmap_refvar_incorrect);
            out.printf("Conf. calls at chip sites correct     : %d\n", hapmap_refvar_correct);
            out.printf("Conf. calls at chip sites incorrect   : %d\n", hapmap_refvar_incorrect);
            out.printf("%% of confident calls that are correct : %.2f%%\n", 100.0 * (float) hapmap_refvar_correct / (float)(hapmap_refvar_correct + hapmap_refvar_incorrect));
        }
        out.println();
    }

    public void printMetricsAtLocusIntervals(int loci_interval) {
        if (num_loci_total % loci_interval == 0) printMetrics();
    }


}
