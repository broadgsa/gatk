package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.*;

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
public class GenotypeConcordance extends BasicVariantAnalysis implements GenotypeAnalysis {

    private static final int REF = 0;
    private static final int VAR_HET = 1;
    private static final int VAR_HOM = 2;
    private static final int UNKNOWN = 3;
    private static final int NO_CALL = 3;   // synonym
    private static final String[] TRUTH_NAMES = {"IS_REF", "IS_VAR_HET", "IS_VAR_HOM", "UNKNOWN"};
    private static final String[] CALL_NAMES = {"CALLED_REF", "CALLED_VAR_HET", "CALLED_VAR_HOM", "NO_CALL"};

    private ArrayList<String> rodNames = new ArrayList<String>();
    private HashMap<String, TruthTable> tables = new HashMap<String, TruthTable>();

    public GenotypeConcordance(final String name) {
        super("genotype_concordance");
        rodNames.add(name);
        initialize();
    }

    public GenotypeConcordance(final List<String> names) {
        super("genotype_concordance");
        rodNames.addAll(names);
        initialize();
    }

    private void initialize() {
        for ( String name : rodNames )
            tables.put(name, new TruthTable());
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        // get all of the chip rods at this locus
        HashMap<String, Genotype> chips = new HashMap<String, Genotype>();
        for ( String name : rodNames ) {
            RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData(name, null);
            Variation chip = (rods == null ? null : (Variation)rods.getRecords().get(0));
            if ( chip != null ) {
                if ( !(chip instanceof VariantBackedByGenotype) )
                    throw new StingException("Failure: trying to analyze genotypes using non-genotype truth data");
                chips.put(name, ((VariantBackedByGenotype)chip).getCalledGenotype());
            }
        }

        if ( eval != null && !(eval instanceof VariantBackedByGenotype) )
            throw new StingException("Failure: trying to analyze genotypes of non-genotype data");

        // don't procede if we have no truth data and no call
        if ( eval != null || chips.size() > 0 )
            inc(chips, (eval != null ? ((VariantBackedByGenotype)eval).getGenotypes() : new ArrayList<Genotype>()), ref);
        return null;
    }

    public void inc(Map<String, Genotype> chips, List<Genotype> evals, char ref) {

        // This shouldn't happen, but let's check anyways
        if (BaseUtils.simpleBaseToBaseIndex(ref) == -1)
            return;

        HashMap<String, Genotype> evalHash = makeGenotypeHash(evals);

        for ( String name : rodNames ) {
            Genotype chip = chips.get(name);
            Genotype eval = evalHash.get(name);

            if (chip == null && eval == null)
                continue;

            int truthType = getGenotypeType(chip, ref);
            int callType = getGenotypeType(eval, ref);
            TruthTable table = tables.get(name);

            //System.out.printf("TEST: %d/%d %s vs. %s%n", truthIndex, callIndex, chip, eval);
            table.addEntry(truthType, callType);
        }
    }

    private HashMap<String, Genotype> makeGenotypeHash(List<Genotype> evals) {
        HashMap<String, Genotype> hash = new HashMap<String, Genotype>();
        for ( Genotype eval : evals ) {
            if ( eval instanceof SampleBacked )
                hash.put(((SampleBacked)eval).getSampleName(), eval);
            else if ( rodNames.size() == 1 )
                hash.put(rodNames.get(0), eval);
            else
                throw new StingException("Genotype data has no associated samples but are multi-sample...?");
        }
        return hash;
    }

    private static int getGenotypeType(Genotype g, char ref) {
        int type;

        if ( g == null )
            type = NO_CALL;
        else if ( !g.isVariant(ref) )
            type = REF;
        else if ( g.isHet() )
            type = VAR_HET;
        else
            type = VAR_HOM;

        return type;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();

        for ( String name : rodNames ) {
            TruthTable table = tables.get(name);
            table.addAllStats(s, name);
            s.add("");
        }

        return s;
    }

    private static String cellPercent(int count, int total) {
        StringBuffer sb = new StringBuffer();
        total = Math.max(total, 0);
        sb.append(String.format("%.2f", (100.0 * count) / total));
        sb.append("%");
        return sb.toString();
    }

    private class TruthTable {
        int[][] table = new int[4][4];
        int[] truth_totals = new int[4];
        int[] calls_totals = new int[4];

        public TruthTable() {
            for (int i = 0; i < 4; i++) {
                truth_totals[i] = 0;
                calls_totals[i] = 0;
                for (int j = 0; j < 4; j++)
                    table[i][j] = 0;
            }
        }

        public void addEntry(int truthIndex, int callIndex) {
            table[truthIndex][callIndex]++;
            truth_totals[truthIndex]++;
            calls_totals[callIndex]++;
        }

        public void addAllStats(List<String> s, String name) {
            s.add(String.format("name                 %s", name));
            s.add(String.format("TRUTH_STATE\tCALLED_REF\tCALLED_VAR_HET\tCALLED_VAR_HOM\tNO_CALL\t\tTOTALS\tTRUE_GENOTYPE_CONCORDANCE\tGENOTYPE_SENSITIVITY"));
            for (int i = 0; i < 4; i++) {
                StringBuffer sb = new StringBuffer();
                sb.append(String.format("%15s ", TRUTH_NAMES[i]));
                for (int j = 0; j < 4; j++)
                    sb.append(String.format("%9d ", table[i][j]));
                sb.append(String.format("%9d ", truth_totals[i]));
                if (i == VAR_HET || i == VAR_HOM) {
                    sb.append(String.format("\t%s\t\t", cellPercent(table[i][i], table[i][REF] + table[i][VAR_HET] + table[i][VAR_HOM])));
                    sb.append(String.format("%s", cellPercent(truth_totals[i] - table[i][NO_CALL], truth_totals[i])));
                } else {
                    sb.append("\tN/A\t\t\tN/A");
                }
                s.add(sb.toString());
            }

            addCalledGenotypeConcordance(s);
            addOverallStats(s);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    s.add(String.format("%s_%s_%s %d", TRUTH_NAMES[i], CALL_NAMES[j], "NO_SITES", table[i][j]));
                    s.add(String.format("%s_%s_%s %s", TRUTH_NAMES[i], CALL_NAMES[j], "PERCENT_OF_TRUTH", cellPercent(table[i][j], truth_totals[i])));
                    s.add(String.format("%s_%s_%s %s", TRUTH_NAMES[i], CALL_NAMES[j], "PERCENT_OF_CALLS", cellPercent(table[i][j], calls_totals[j])));
                }
                if (i == VAR_HET || i == VAR_HOM) {
                    s.add(String.format("%s_%s %s", TRUTH_NAMES[i], "TRUE_GENOTYPE_CONCORDANCE", cellPercent(table[i][i], table[i][REF] + table[i][VAR_HET] + table[i][VAR_HOM])));
                    s.add(String.format("%s_%s %s", TRUTH_NAMES[i], "GENOTYPE_SENSITIVITY", cellPercent(truth_totals[i] - table[i][NO_CALL], truth_totals[i])));
                }
            }
        }

        private void addCalledGenotypeConcordance(List<String> s) {
            StringBuilder sb = new StringBuilder();
            sb.append("CALLED_GENOTYPE_CONCORDANCE\t");
            for (int i = 0; i < 4; i++) {
                int nConcordantCallsI = table[i][i];
                String value = "N/A";
                if (i != UNKNOWN)
                    value = String.format("%s\t", cellPercent(nConcordantCallsI, calls_totals[i] - table[UNKNOWN][i]));
                sb.append(value);
            }
            s.add(sb.toString());
        }

        // How many overall calls where made that aren't NO_CALLS or UNKNOWNS?
        private int getNCalled() {
            int n = 0;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    if (i != NO_CALL && j != NO_CALL) n += table[i][j];
            return n;
        }

        private void addOverallStats(List<String> s) {
            int nConcordantRefCalls = table[REF][REF];
            int nConcordantHetCalls = table[VAR_HET][VAR_HET];
            int nConcordantVarHomCalls = table[VAR_HOM][VAR_HOM];
            int nVarCalls = table[VAR_HOM][VAR_HET] + table[VAR_HOM][VAR_HOM] + table[VAR_HET][VAR_HET] + table[VAR_HET][VAR_HOM];
            int nConcordantVarCalls = nConcordantHetCalls + nConcordantVarHomCalls;
            int nConcordantCalls = nConcordantRefCalls + nConcordantVarCalls;
            int nTrueVar = truth_totals[VAR_HET] + truth_totals[VAR_HOM];
            int nCalled = getNCalled();
            s.add(String.format("VARIANT_SENSITIVITY %s", cellPercent(nVarCalls, nTrueVar)));
            s.add(String.format("VARIANT_CONCORDANCE %s", cellPercent(nConcordantVarCalls, nVarCalls)));
            s.add(String.format("OVERALL_CONCORDANCE %s", cellPercent(nConcordantCalls, nCalled)));
        }
    }
}
