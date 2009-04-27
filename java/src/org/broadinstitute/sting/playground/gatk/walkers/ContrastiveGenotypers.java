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

import java.io.PrintStream;

public class ContrastiveGenotypers extends LocusWalker<Integer, Integer> {
    private SingleSampleGenotyper caller_1b;
    private SingleSampleGenotyper caller_4b;

    public void initialize() {
        caller_1b = new SingleSampleGenotyper();
        caller_1b.metricsFileName = "/dev/stdout";
        caller_1b.metricsInterval = 6000;
        caller_1b.printMetrics = false;
        caller_1b.fourBaseMode = false;
        caller_1b.retest = false;
        caller_1b.qHom = 0.04;
        caller_1b.qHet = 0.49;
        caller_1b.qHomNonRef = 0.97;
        caller_1b.lodThreshold = 5.0;
        caller_1b.initialize();
        caller_1b.reduceInit();

        caller_4b = new SingleSampleGenotyper();
        caller_4b.metricsFileName = "/dev/stdout";
        caller_4b.metricsInterval = caller_1b.metricsInterval;
        caller_4b.printMetrics = false;
        caller_4b.fourBaseMode = true;
        caller_4b.retest = false;
        caller_4b.qHom = 0.04;
        caller_4b.qHet = 0.49;
        caller_4b.qHomNonRef = 0.97;
        caller_4b.lodThreshold = 5.0;
        caller_4b.initialize();
        caller_4b.reduceInit();
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ref = Character.toUpperCase(ref);

        rodGFF hmi = getHapmapInfo(tracker);

        AlleleFrequencyEstimate call_1b = caller_1b.map(tracker, ref, context);
        AlleleFrequencyEstimate call_4b = caller_4b.map(tracker, ref, context);

        //if (isHomRefHapmapSite(ref, hmi) && isHomRef(call_1b) && isHet(call_4b)) {
        //    printDebuggingInfo(ref, context, call_1b, call_4b, hmi, System.out);
        //}

        //
        caller_1b.metrics.nextPosition(call_1b, tracker);
        caller_4b.metrics.nextPosition(call_4b, tracker);

        if (caller_1b.metrics.num_loci_total % caller_1b.metricsInterval == 0) {
            System.out.print("1-Base Metrics");
            caller_1b.metrics.printMetricsAtLocusIntervals(caller_1b.metricsInterval);
            System.out.print("4-Base Metrics");
            caller_4b.metrics.printMetricsAtLocusIntervals(caller_1b.metricsInterval);
            System.out.println("--------------");
        }
        //

        return null;
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

    private void printDebuggingInfo(char ref, LocusContext context, AlleleFrequencyEstimate call_1b, AlleleFrequencyEstimate call_4b, rodGFF hmi, PrintStream lout) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);

        int pRefCount = 0, pAltCount = 0, pNonRefCount = 0;
        int sRefCount = 0, sAltCount = 0, sNonRefCount = 0;
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
        }

        lout.println("Locus: " + context.getLocation());
        lout.println("Depth: " + call_4b.depth);
        if (hmi != null) { lout.println("Hapmap Info: " + hmi.toString()); }
        lout.println("Primary Base pileup:   " + pileup.getBases());
        lout.println("Secondary Base pileup: " + pileup.getSecondaryBasePileup());
        lout.println();

        lout.printf("Pct primary bases that are ref:      %2.2f%%\n", 100.0*((double) pRefCount)/((double) pileup.getBases().length()));
        lout.printf("Pct primary bases that are alt:      %2.2f%%\n", 100.0*((double) pAltCount)/((double) pileup.getBases().length()));
        lout.printf("Pct primary bases that are nonRef:   %2.2f%%\n", 100.0*((double) pNonRefCount)/((double) pileup.getBases().length()));
        lout.println();

        lout.printf("Pct secondary bases that are ref:    %2.2f%%\n", 100.0*((double) sRefCount)/((double) pileup.getBases().length()));
        lout.printf("Pct secondary bases that are alt:    %2.2f%%\n", 100.0*((double) sAltCount)/((double) pileup.getBases().length()));
        lout.printf("Pct secondary bases that are nonRef: %2.2f%%\n", 100.0*((double) sNonRefCount)/((double) pileup.getBases().length()));
        lout.println();

        lout.printf("1-base result: ref=%c alt=%c q=%2.2f lodBestVsRef=%4.4f lodBestVsNextBest=%4.4f\n", call_1b.ref, call_1b.alt, call_1b.qhat, call_1b.lodVsRef, call_1b.lodVsNextBest);
        lout.printf("4-base result: ref=%c alt=%c q=%2.2f lodBestVsRef=%4.4f lodBestVsNextBest=%4.4f\n", call_4b.ref, call_4b.alt, call_4b.qhat, call_4b.lodVsRef, call_4b.lodVsNextBest);
        lout.println();

        String[] names = { "Homozygous reference", "Heterozygous", "Homozygous non-reference" };
        for (int i = 0; i < call_4b.posteriors.length; i++) {
            lout.println(names[i] + " posterior: " + call_4b.posteriors[i]);
        }

        //lout.println("\nProb pileup:\n" + pileup.getProbDistPileup());
        lout.println("\n------------------------------------------------------------------------------\n");
    }

    public Integer reduceInit() { return null; }
    public Integer reduce(Integer value, Integer sum) { return null; }

    public void onTraversalDone(Integer result) {
    }
}
