package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.ArrayList;
import java.util.List;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class CallableBasesAnalysis extends BasicVariantAnalysis implements GenotypeAnalysis {
    long all_bases = 0;
    long all_calls = 0;
    final static double[] thresholds = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50, 100};
    long[] discoverable_bases = new long[thresholds.length];
    long[] genotypable_bases = new long[thresholds.length];

    public CallableBasesAnalysis() {
        super("callable_bases");
    }

    public long nSites() {
        return all_bases;
    }

    public long nCalls() {
        return all_calls;
    }

    public long nDiscoverable(int index) {
        return discoverable_bases[index];
    }

    public double percentDiscoverable(int index) {
        return (100.0 * nDiscoverable(index)) / nSites();
    }

    public long nGenotypable(int index) {
        return genotypable_bases[index];
    }

    public double percentGenotypable(int index) {
        return (100.0 * nGenotypable(index)) / nSites();
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        all_bases++;

        if (eval == null)                     // no data here!
            return null;

        // we actually have a record here
        if (!(eval instanceof VariantBackedByGenotype)) {            // evaluation record isn't a genotype, die!
            throw new RuntimeException("Evaluation track isn't backed by a Genotype!");
        }

        all_calls++;
        // For every threshold, updated discoverable and callable
        for (int i = 0; i < thresholds.length; i++) {
            double threshold = thresholds[i];
            DiploidGenotype g = DiploidGenotype.createHomGenotype(ref);
            Genotype genotype = ((VariantBackedByGenotype) eval).getGenotype(g);
            // update discoverable
            if (eval.isSNP() && eval.getNegLog10PError() >= threshold)
                discoverable_bases[i]++;
            if (!eval.isSNP() && genotype.getNegLog10PError() >= threshold)
                discoverable_bases[i]++;
            if (genotype.getNegLog10PError() >= threshold)
                genotypable_bases[i]++;
        }

        return null;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("total_no_sites       %d", nSites()));
        s.add(String.format("total_no_calls       %d", nCalls()));
        s.add(String.format(""));
        s.add(String.format("confidence_threshold\tdiscoverable_bases\tdiscoverable_bases_percent\tgenotypable_bases\tgenotypable_bases_percent"));

        for (int i = 0; i < thresholds.length; i++) {
            double threshold = thresholds[i];
            s.add(String.format("%6.2f\t%d\t%.6f\t%d\t%.6f", threshold, nDiscoverable(i), percentDiscoverable(i), nGenotypable(i), percentGenotypable(i)));
        }

        return s;
    }
}