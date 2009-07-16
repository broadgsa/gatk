package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public class CallableBasesAnalysis extends BasicVariantAnalysis implements GenotypeAnalysis {
    long all_bases = 0;
    long all_calls = 0;
    //final static double[] Qthresholds = { 10, 20, 30, 40, 50, 100, 200, 500, 1000 };
    final static double[] thresholds = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50, 100 };
    long[] discoverable_bases = new long[thresholds.length];
    long[] genotypable_bases = new long[thresholds.length];

    public CallableBasesAnalysis() {
        super("callable_bases");
    }

    public long nSites()                    { return all_bases; }
    public long nCalls()                    { return all_calls; }
    public long nDiscoverable(int index)    { return discoverable_bases[index]; }
    public double percentDiscoverable(int index)     { return (100.0*nDiscoverable(index)) / nSites(); }
    public long nGenotypable(int index)     { return genotypable_bases[index]; }
    public double percentGenotypable(int index)     { return (100.0*nGenotypable(index)) / nSites(); }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        all_bases++;

        if ( eval == null )                     // no data here!
            return null;

        // we actually have a record here
        if ( ! eval.isGenotype() ) {            // evaluation record isn't a genotype, die!
            throw new RuntimeException("Evaluation track isn't an Genotype!");
        }

        all_calls++;
        // For every threshold, updated discoverable and callable
        for ( int i = 0; i < thresholds.length; i++ ) {
            double threshold = thresholds[i];

            // update discoverable
            if ( eval.isSNP() && eval.getVariationConfidence() >= threshold )
                discoverable_bases[i]++;
            if ( eval.isReference() && eval.getConsensusConfidence() >= threshold )
                discoverable_bases[i]++;

            if ( eval.getConsensusConfidence() >= threshold )
                genotypable_bases[i]++;

            //System.out.printf("Updating %s SNP=%b, REF=%b VarConf=%f ConConf=%f where threshold=%f: discoverable = %d, genotypable = %d%n",
            //        eval.getLocation(), eval.isSNP(), eval.isReference(), eval.getVariationConfidence(),
            //        eval.getConsensusConfidence(), threshold, discoverable_bases[i], genotypable_bases[i]);
        }

        return null;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("total_no_sites       %d", nSites()));
        s.add(String.format("total_no_calls       %d", nCalls()));
        s.add(String.format(""));
        s.add(String.format("confidence_threshold\tdiscoverable_bases\tdiscoverable_bases_percent\tgenotypable_bases\tgenotypable_bases_percent"));

        for ( int i = 0; i < thresholds.length; i++ ) {
            double threshold = thresholds[i];
            s.add(String.format("%6.2f\t%d\t%.6f\t%d\t%.6f", threshold, nDiscoverable(i), percentDiscoverable(i), nGenotypable(i), percentGenotypable(i)));
        }

        return s;
    }
}