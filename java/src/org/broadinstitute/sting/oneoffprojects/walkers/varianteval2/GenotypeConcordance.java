package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class GenotypeConcordance extends VariantEvaluator {

    private static final int nGenotypeTypes = Genotype.Type.values().length;

    class SampleStats {

        long[][] concordance = new long[nGenotypeTypes][nGenotypeTypes];

        public String toString() {
            StringBuffer sb = new StringBuffer();
            for (int truth = 0; truth < nGenotypeTypes; truth++ ) {
                // don't print out when truth = no-call
                if ( truth == Genotype.Type.NO_CALL.ordinal() )
                    continue;
                long total = 0;
                for (int called = 0; called < nGenotypeTypes; called++ )
                    total += concordance[truth][called];
                sb.append(String.format("%d %d %.2f ", total, concordance[truth][truth], total == 0 ? 0.0 : (100.0 * (double)concordance[truth][truth] / (double)total)));
                for (int called = 0; called < nGenotypeTypes; called++ ) {
                    if ( called != truth )
                        sb.append(String.format("%d ", concordance[truth][called]));
                }
            }

            return sb.toString();
        }
    }

    class FrequencyStats {
        long nFound = 0;
        long nMissed = 0;

        public String toString() {
            long total = nFound + nMissed;
            return String.format("%d %d %.2f ", nFound, nMissed, total == 0 ? 0.0 : (100.0 * rate(nFound, total)));
        }
    }

    // a mapping from allele count to stats
    private HashMap<Integer, FrequencyStats> alleleCountStats = new HashMap<Integer, FrequencyStats>();

    // a mapping from sample to stats
    private HashMap<String, SampleStats> concordanceStats = null;

    // keep a list of the validation data we saw before the first eval data
    private HashSet<VariantContext> missedValidationData = new HashSet<VariantContext>();


    public GenotypeConcordance(VariantEval2Walker parent) {
        // don't do anything
    }

    public String getName() {
        return "genotypeConcordance";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public boolean enabled() { return true; }

    public String toString() {
        return getName() + ": " + getTableRows();
    }

    private static List<String> SAMPLE_HEADER =
            Arrays.asList("sample",
                    "total_true_ref", "n_ref/ref", "%_ref/ref",
                    "n_ref/no-call", "n_ref/het", "n_ref/hom",
                    "total_true_het", "n_het/het", "%_het/het",
                    "n_het/no-call", "n_het/ref", "n_het/hom",
                    "total_true_hom", "n_hom/hom", "%_hom/hom",
                    "n_hom/no-call", "n_hom/ref", "n_hom/het");

    private static List<String> FREQUENCY_HEADER =
            Arrays.asList("alleleCount", "n_found", "n_missed", "%_found");

    // making it a table
    public List<String> getTableHeader() {
        ArrayList<String> header = new ArrayList<String>();
        header.addAll(SAMPLE_HEADER);
        header.addAll(FREQUENCY_HEADER);
        return header;
    }

    public List<List<String>> getTableRows() {
        ArrayList<List<String>> rows = new ArrayList<List<String>>();

        if ( concordanceStats != null ) {
            for ( Map.Entry<String, SampleStats> sample : concordanceStats.entrySet() )
                rows.add(Arrays.asList(String.format("%s %s", sample.getKey(), sample.getValue().toString()).split(" ")));
        }

        if ( alleleCountStats != null ) {
            for ( Map.Entry<Integer, FrequencyStats> alleleCount : alleleCountStats.entrySet() )
                rows.add(Arrays.asList(String.format("%d %s", alleleCount.getKey(), alleleCount.getValue().toString()).split(" ")));
        }

        return rows;
    }

    public String update2(VariantContext eval, VariantContext validation, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String interesting = null;

        // sanity check that we at least have either eval or validation data
        if ( eval == null && !isValidVC(validation) )
            return interesting;

        if ( concordanceStats == null ) {
            if ( eval != null ) {
                // initialize the concordance table
                createConcordanceTable(eval);
                for ( VariantContext vc : missedValidationData )
                    determineStats(null, vc);
                missedValidationData = null;
            } else {
                missedValidationData.add(validation);
                return interesting;
            }
        }

        determineStats(eval, validation);
               
        return interesting; // we don't capture any interesting sites
    }

    private void determineStats(VariantContext eval, VariantContext validation) {

        // determine concordance for eval data
        if ( eval != null ) {

            for ( String sample : eval.getSampleNames() ) {
                Genotype.Type called = eval.getGenotype(sample).getType();
                Genotype.Type truth;

                if ( !isValidVC(validation) || !validation.hasGenotype(sample) )
                    truth = Genotype.Type.NO_CALL;
                else
                    truth = validation.getGenotype(sample).getType();

                SampleStats stats = concordanceStats.get(sample);
                if ( stats == null )
                    throw new StingException("Sample " + sample + " has not been seen in a previous eval; this analysis module assumes that all samples are present in each variant context");
                stats.concordance[truth.ordinal()][called.ordinal()]++;
            }
        }
        // otherwise, mark no-calls for all samples
        else {

            Genotype.Type called = Genotype.Type.NO_CALL;

            for ( String sample : validation.getSampleNames() ) {
                SampleStats stats = concordanceStats.get(sample);
                if ( stats == null )
                    continue;

                Genotype.Type truth = validation.getGenotype(sample).getType();
                stats.concordance[truth.ordinal()][called.ordinal()]++;
            }
        }

        // determine allele count concordance ()
        if ( isValidVC(validation) && validation.isPolymorphic() ) {
            int trueAlleleCount = 0;
            for ( Allele a : validation.getAlternateAlleles() )
                trueAlleleCount += validation.getChromosomeCount(a);

            if ( !alleleCountStats.containsKey(trueAlleleCount) )
                alleleCountStats.put(trueAlleleCount, new FrequencyStats());
            FrequencyStats stats = alleleCountStats.get(trueAlleleCount);
            if ( eval != null )
                stats.nFound++;
            else
                stats.nMissed++;
        }
    }

    private static boolean isValidVC(VariantContext vc) {
        return (vc != null && !vc.isFiltered());
    }

    private void createConcordanceTable(VariantContext vc) {
        concordanceStats = new HashMap<String, SampleStats>();
        for ( String sample : vc.getSampleNames() )
            concordanceStats.put(sample, new SampleStats());
    }

}