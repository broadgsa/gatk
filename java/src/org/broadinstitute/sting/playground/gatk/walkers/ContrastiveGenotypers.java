package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.PrintStream;
import java.io.File;

public class ContrastiveGenotypers extends LocusWalker<Integer, Integer> {
    private SingleSampleGenotyper caller_1b;
    private SingleSampleGenotyper caller_4b;

    String[] genotypes = { "AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT" };
    int[][] goodAltBreakdown = new int[10][2];
    int[][] badAltBreakdown = new int[10][2];

    public void initialize() {
        caller_1b = new SingleSampleGenotyper();
        caller_1b.METRICS_FILE = new File("/dev/stdout");
        caller_1b.METRICS_INTERVAL = 6000;
        caller_1b.NO_SECONDARY_BASES = false;
        //caller_1b.retest = false;
        //caller_1b.qHom = 0.04;
        //caller_1b.qHet = 0.49;
        //caller_1b.qHomNonRef = 0.97;
        caller_1b.LOD_THRESHOLD = 5.0;
        caller_1b.initialize();
        caller_1b.reduceInit();

        caller_4b = new SingleSampleGenotyper();
        caller_4b.METRICS_FILE = new File("/dev/stdout");
        caller_4b.METRICS_INTERVAL = caller_1b.METRICS_INTERVAL;
        caller_4b.NO_SECONDARY_BASES = true;
        //caller_4b.retest = true;
        //caller_4b.qHom = 0.04;
        //caller_4b.qHet = 0.49;
        //caller_4b.qHomNonRef = 0.97;
        caller_4b.LOD_THRESHOLD = 5.0;
        caller_4b.initialize();
        caller_4b.reduceInit();
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ref = Character.toUpperCase(ref);

        rodGFF hmi = getHapmapInfo(tracker);
        boolean isInDbSNP = isInDbSNP(tracker);

        if (hmi != null) {
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
            String bases = pileup.getBases();
            String bases2 = pileup.getSecondaryBasePileup();

            int genotypeIndex = getGenotypeIndex(hmi.getFeature());

            for (int pileupIndex = 0; pileupIndex < bases.length(); pileupIndex++) {
                if (hmi.getFeature().charAt(0) != bases.charAt(pileupIndex) && hmi.getFeature().charAt(1) != bases.charAt(pileupIndex)) {
                    if (hmi.getFeature().charAt(0) == bases2.charAt(pileupIndex) || hmi.getFeature().charAt(1) == bases2.charAt(pileupIndex)) {
                        badAltBreakdown[genotypeIndex][0]++;
                    } else {
                        badAltBreakdown[genotypeIndex][1]++;
                    }
                } else {
                    if (hmi.getFeature().charAt(0) == bases2.charAt(pileupIndex) || hmi.getFeature().charAt(1) == bases2.charAt(pileupIndex)) {
                        goodAltBreakdown[genotypeIndex][0]++;
                    } else {
                        goodAltBreakdown[genotypeIndex][1]++;
                    }
                }
            }
        }

        AlleleFrequencyEstimate call_1b = caller_1b.map(tracker, ref, context);
        AlleleFrequencyEstimate call_4b = caller_4b.map(tracker, ref, context);

        //if (call_1b.qhat != call_4b.qhat && call_4b.lodVsNextBest >= 5.0) {
        //    printDebuggingInfo(ref, context, call_1b, call_4b, hmi, System.out);
        //}

        String calltype_1b = "homref";
        if (isHet(call_1b)) { calltype_1b = "het"; }
        if (isHomNonRef(call_1b)) { calltype_1b = "homnonref"; }

        String calltype_4b = "homref";
        if (isHet(call_4b)) { calltype_4b = "het"; }
        if (isHomNonRef(call_4b)) { calltype_4b = "homnonref"; }

        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        String bases2 = pileup.getSecondaryBasePileup();

        out.printf("%s %c %s %s %4.4f %4.4f %s %s %4.4f %4.4f %b %b %s %s %s\n",
                          context.getLocation().toString(),
                          ref,
                          calltype_1b,
                          call_1b.genotype(),
                          call_1b.lodVsRef,
                          call_1b.lodVsNextBest,
                          calltype_4b,
                          call_4b.genotype(),
                          call_4b.lodVsRef,
                          call_4b.lodVsNextBest,
                          isInDbSNP,
                          hmi != null,
                          (hmi != null) ? hmi.getFeature() : "NN",
                          bases,
                          bases2
                          );

        caller_1b.metricsOut.nextPosition(call_1b, tracker);
        caller_4b.metricsOut.nextPosition(call_4b, tracker);

        /*
        if (caller_1b.metricsOut.num_loci_total % caller_1b.METRICS_INTERVAL == 0) {
            System.out.print("1-Base Metrics");
            caller_1b.metrics.printMetricsAtLocusIntervals(caller_1b.METRICS_INTERVAL);
            System.out.print("4-Base Metrics");
            caller_4b.metrics.printMetricsAtLocusIntervals(caller_1b.METRICS_INTERVAL);
            System.out.println("--------------");
        }
        */

        return null;
    }

    int getGenotypeIndex(String genotype) {
        for (int genotypeIndex = 0; genotypeIndex < genotypes.length; genotypeIndex++) {
            if (genotypes[genotypeIndex].matches(genotype)) {
                return genotypeIndex;
            }
        }

        return -1;
    }

    private boolean isHomRef(AlleleFrequencyEstimate freq) { return MathUtils.compareDoubles(freq.qhat, 0.0) == 0; }
    private boolean isHet(AlleleFrequencyEstimate freq) { return MathUtils.compareDoubles(freq.qhat, 0.5) == 0; }
    private boolean isHomNonRef(AlleleFrequencyEstimate freq) { return MathUtils.compareDoubles(freq.qhat, 1.0) == 0; }

    private boolean isHomRefHapmapSite(char ref, rodGFF hmi) { return (hmi != null && hmi.getFeature().charAt(0) == hmi.getFeature().charAt(1) && hmi.getFeature().charAt(0) == ref); }
    private boolean isHetHapmapSet(char ref, rodGFF hmi) { return hmi != null && !isHomRefHapmapSite(ref, hmi) && !isHomNonRefHapmapSite(ref, hmi); }
    private boolean isHomNonRefHapmapSite(char ref, rodGFF hmi) { return (hmi != null && hmi.getFeature().charAt(0) == hmi.getFeature().charAt(1) && hmi.getFeature().charAt(0) != ref); }

    private rodGFF getHapmapInfo(RefMetaDataTracker tracker) {
        rodGFF hapmap_chip_genotype = null;
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) {
            if ( datum != null && datum instanceof rodGFF ) {
                hapmap_chip_genotype = (rodGFF) datum;
            }
        }

        return hapmap_chip_genotype;
    }

    private boolean isInDbSNP(RefMetaDataTracker tracker) {
        boolean is_dbSNP_SNP = false;

        for ( ReferenceOrderedDatum datum : tracker.getAllRods() )
        {
            if ( datum != null )
            {
                if ( datum instanceof rodDbSNP)
                {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    if (dbsnp.isSNP()) is_dbSNP_SNP = true;
                }
            }
        }

        return is_dbSNP_SNP;
    }

    private void printDebuggingInfo(char ref, LocusContext context, AlleleFrequencyEstimate call_1b, AlleleFrequencyEstimate call_4b, rodGFF hmi, PrintStream lout) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);

        int pRefCount = 0, pAltCount = 0, pNonRefCount = 0;
        int sRefCount = 0, sAltCount = 0, sNonRefCount = 0;
        int sAllelicCondCount = 0, sAllelicCondTotal = 0;
        int[] secondaryBaseCounts = new int[4];

        for (int i = 0; i < pileup.getBases().length(); i++) {
            if   (pileup.getBases().charAt(i) == call_4b.ref) { pRefCount++; }
            else {
                if (pileup.getBases().charAt(i) == call_4b.alt) { pAltCount++; }
                pNonRefCount++;
            }

            if   (pileup.getSecondaryBasePileup().charAt(i) == call_4b.ref) { sRefCount++; }
            else {
                if (pileup.getSecondaryBasePileup().charAt(i) == call_4b.alt) { sAltCount++; }
                sNonRefCount++;
            }

            if (pileup.getBases().charAt(i) != call_4b.ref && pileup.getBases().charAt(i) != call_4b.alt) {
                sAllelicCondTotal++;
                if   (pileup.getSecondaryBasePileup().charAt(i) == call_4b.ref || pileup.getSecondaryBasePileup().charAt(i) == call_4b.alt) {
                    sAllelicCondCount++;
                }
            }

            secondaryBaseCounts[BaseUtils.simpleBaseToBaseIndex(pileup.getSecondaryBasePileup().charAt(i))]++;
        }

        lout.println("Locus: " + context.getLocation());
        lout.println("Depth: " + call_4b.depth);
        if (hmi != null) { lout.println("Hapmap Info: " + hmi.toString()); }
        lout.println("Primary Base pileup:   " + pileup.getBases());
        lout.println("Secondary Base pileup: " + pileup.getSecondaryBasePileup());
        lout.printf("Secondary Base Counts: A:%d C:%d G:%d T:%d\n", secondaryBaseCounts[0], secondaryBaseCounts[1], secondaryBaseCounts[2], secondaryBaseCounts[3]);
        lout.println();

        lout.printf("Pct primary bases that are ref:      %2.2f%%\n", 100.0*((double) pRefCount)/((double) pileup.getBases().length()));
        lout.printf("Pct primary bases that are alt:      %2.2f%%\n", 100.0*((double) pAltCount)/((double) pileup.getBases().length()));
        lout.printf("Pct primary bases that are nonRef:   %2.2f%%\n", 100.0*((double) pNonRefCount)/((double) pileup.getBases().length()));
        lout.println();

        lout.printf("Pct secondary bases that are ref:    %2.2f%%\n", 100.0*((double) sRefCount)/((double) pileup.getBases().length()));
        lout.printf("Pct secondary bases that are alt:    %2.2f%%\n", 100.0*((double) sAltCount)/((double) pileup.getBases().length()));
        lout.printf("Pct secondary bases that are nonRef: %2.2f%%\n", 100.0*((double) sNonRefCount)/((double) pileup.getBases().length()));
        lout.println();

        lout.printf("Pct secondary bases that match one of the alleles when the primary bases don't: %2.2f%%\n", 100.0*((double) sAllelicCondCount)/((double) sAllelicCondTotal));

        lout.printf("1-base result: ref=%c alt=%c q=%2.2f lodBestVsRef=%4.4f lodBestVsNextBest=%4.4f\n", call_1b.ref, call_1b.alt, call_1b.qhat, call_1b.lodVsRef, call_1b.lodVsNextBest);
        lout.printf("4-base result: ref=%c alt=%c q=%2.2f lodBestVsRef=%4.4f lodBestVsNextBest=%4.4f\n", call_4b.ref, call_4b.alt, call_4b.qhat, call_4b.lodVsRef, call_4b.lodVsNextBest);
        lout.println();

        /*
        String[] names = { "Homozygous reference", "Heterozygous", "Homozygous non-reference" };
        for (int i = 0; i < call_4b.posteriors.length; i++) {
            lout.println(names[i] + " posterior: " + call_4b.posteriors[i]);
        }
        */

        //lout.println("\nProb pileup:\n" + pileup.getProbDistPileup());
        lout.println("\n------------------------------------------------------------------------------\n");
    }


    public Integer reduceInit() { return null; }
    public Integer reduce(Integer value, Integer sum) { return null; }

    public void onTraversalDone(Integer result) {
        for (int genotypeIndex = 0; genotypeIndex < genotypes.length; genotypeIndex++) {
            System.out.printf("%s: good %d %d ( %3.3f %3.3f ); bad %d %d ( %3.3f %3.3f )\n",
                              genotypes[genotypeIndex],
                              
                              goodAltBreakdown[genotypeIndex][0],
                              goodAltBreakdown[genotypeIndex][1],
                              ((double) goodAltBreakdown[genotypeIndex][0])/((double) (goodAltBreakdown[genotypeIndex][0] + goodAltBreakdown[genotypeIndex][1])),
                              ((double) goodAltBreakdown[genotypeIndex][1])/((double) (goodAltBreakdown[genotypeIndex][0] + goodAltBreakdown[genotypeIndex][1])),

                              badAltBreakdown[genotypeIndex][0],
                              badAltBreakdown[genotypeIndex][1],
                              ((double) badAltBreakdown[genotypeIndex][0])/((double) (badAltBreakdown[genotypeIndex][0] + badAltBreakdown[genotypeIndex][1])),
                              ((double) badAltBreakdown[genotypeIndex][1])/((double) (badAltBreakdown[genotypeIndex][0] + badAltBreakdown[genotypeIndex][1]))

                              );
        }

        System.out.print("1-Base Metrics");
        caller_1b.metricsOut.printMetricsAtLocusIntervals(1);
        System.out.print("4-Base Metrics");
        caller_4b.metricsOut.printMetricsAtLocusIntervals(1);
        System.out.println("--------------");
    }
}
