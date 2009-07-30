package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.utils.genotype.calls.GenotypeCall;
import org.broadinstitute.sting.utils.genotype.calls.SSGGenotypeCall;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

public class AlleleMetrics {
    private double LOD_cutoff = 5;
    
    private long dbsnp_hits = 0;
    private long num_variants = 0;
    private long num_loci_total = 0;
    private long num_loci_confident = 0;
    private long hapmap_genotype_correct = 0;
    private long hapmap_genotype_incorrect = 0;
    private long hapmap_refvar_correct = 0;
    private long hapmap_refvar_incorrect = 0;

    private final double dbl_cmp_precision = 0.0001;

    protected PrintStream out;

    public AlleleMetrics(String metricsOutputPath) {
        initialize(new File(metricsOutputPath), LOD_cutoff);
    }

    public AlleleMetrics(String metricsOutputPath, double lodThreshold) {
        initialize(new File(metricsOutputPath), lodThreshold);
    }

    public AlleleMetrics(File metricsOutputFile) {
        initialize(metricsOutputFile, LOD_cutoff);
    }

    public AlleleMetrics(File metricsOutputFile, double lodThreshold) {
        initialize(metricsOutputFile, lodThreshold);
    }

    private void initialize(File metricsOutputFile, double lodThresold) {
        try {
            this.out = new PrintStream(metricsOutputFile);
        } catch (FileNotFoundException e) {
            System.err.format("Unable to open file '%s'. Perhaps the parent directory does not exist or is read-only.", metricsOutputFile.getAbsolutePath());
            System.exit(-1);
        }

        this.LOD_cutoff = lodThresold;
    }

    public void nextPosition(GenotypeCall cl, RefMetaDataTracker tracker) {
        SSGGenotypeCall call = (SSGGenotypeCall)cl;
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
        double result = call.getBestRef();
        if (Math.abs(call.getBestNext()) >= LOD_cutoff) { num_loci_confident += 1; }

        if (call.isVariation() && result >= LOD_cutoff)
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
            String str = call.getBases();
            char alt = str.charAt(0);
            if (str.charAt(0) == call.getReferencebase()) alt = str.charAt(1);
            for (char c : hapmap_genotype.toCharArray()) {
                if (c == call.getReferencebase()) { refs++; }
                if (c == alt) { alts++; }
            }

            if (refs+alts > 0) {
                hapmap_q = ((double) alts) / ((double) (refs+alts));
            }else{
                hapmap_q = -1.0;
            }

            // Hapmap debug info
            //out.format("HAPMAP DEBUG %.2f %.2f %.2f ", hapmap_q, alleleFreq.qstar, alleleFreq.lodVsRef);
            //String called_genotype = alleleFreq.asString();
            //out.format("%s %s %c %c", hapmap_genotype, called_genotype, alleleFreq.ref, alleleFreq.alt);

            //System.out.printf("DBG %f %s\n", LOD_cutoff, alleleFreq.asTabularString());
            if (call.getBestNext() >= LOD_cutoff) {
                // TODO : should this all be commented out?
                /*
                System.out.printf("DBG %f %f %f %f\n",
                                        hapmap_q,
                                        alleleFreq.qhat,
                                        alleleFreq.qstar,
                                        alleleFreq.lodVsNextBest);
                */

                // Calculate genotyping performance - did we get the correct genotype of the N+1 choices?
                //if (hapmap_q != -1 && hapmap_q == alleleFreq.qstar) {
                /*if (Math.abs(hapmap_q - -1.0) > dbl_cmp_precision && Math.abs(hapmap_q - alleleFreq.qstar) <= dbl_cmp_precision) {
                    hapmap_genotype_correct++;
                }else{
                    hapmap_genotype_incorrect++;
                    //System.out.printf(" INCORRECT GENOTYPE    Bases: %s", AlleleFrequencyWalker.getBases(context));
                    //out.printf(" INCORRECT GENOTYPE");
                    //AlleleFrequencyWalker.print_base_qual_matrix(AlleleFrequencyWalker.getOneBaseQuals(context));
                }*/
            }

            if (result >= LOD_cutoff || -1 * result >= LOD_cutoff) {

                // Now calculate ref / var performance - did we correctly classify the site as
                // reference or variant without regard to genotype; i.e. het/hom "miscalls" don't matter here
                boolean hapmap_var = hapmap_q != 0.0;
                boolean called_var = call.isVariation();
                //if (hapmap_q != -1 && hapmap_var != called_var) {
                if (Math.abs(hapmap_q - -1.0) > dbl_cmp_precision && hapmap_var != called_var) {
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
