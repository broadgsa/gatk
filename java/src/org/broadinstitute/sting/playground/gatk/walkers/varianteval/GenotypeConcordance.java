package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
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
public class GenotypeConcordance extends BasicVariantAnalysis implements GenotypeAnalysis {
    private String dbName;

    private static final int REF = 0;
    private static final int VAR_HET = 1;
    private static final int VAR_HOM = 2;
    private static final int UNKNOWN = 3;
    private static final int NO_CALL = 3;   // synonym
    private static final String[] TRUTH_NAMES = {"IS_REF", "IS_VAR_HET", "IS_VAR_HOM", "UNKNOWN"};
    private static final String[] CALL_NAMES = {"CALLED_REF", "CALLED_VAR_HET", "CALLED_VAR_HOM", "NO_CALL"};

    private int[][] table = new int[4][4];
    private int[] truth_totals = new int[4];
    private int[] calls_totals = new int[4];

    public GenotypeConcordance(final String name) {
        super("genotype_concordance");
        dbName = name;
        for (int i = 0; i < 4; i++) {
            truth_totals[i] = 0;
            calls_totals[i] = 0;
            for (int j = 0; j < 4; j++)
                table[i][j] = 0;
        }
    }

    public void inc(Variation chip, Variation eval, char ref) {
        if ((chip != null && !(chip instanceof VariantBackedByGenotype) || (eval != null && !(eval instanceof VariantBackedByGenotype))))
            throw new StingException("Failure: trying to analyze genotypes of non-genotype data");

	// This shouldn't happen, but let's check anyways
	if ( BaseUtils.simpleBaseToBaseIndex(ref) == -1 )
	    return;

        DiploidGenotype g = DiploidGenotype.createHomGenotype(ref);
        int truthIndex, callIndex;
        if (chip == null)
            truthIndex = UNKNOWN;
        else if (chip.getAlternateBases().equals(g.toString()))
            truthIndex = REF;
        else if (chip.getAlternateBases().charAt(0) != chip.getAlternateBases().charAt(1))
            truthIndex = VAR_HET;
        else
            truthIndex = VAR_HOM;

        if (eval == null)
            callIndex = NO_CALL;
        else if (eval.getAlternateBases().equals(g.toString()))
            callIndex = REF;
        else if (eval.getAlternateBases().charAt(0) != eval.getAlternateBases().charAt(1))
            callIndex = VAR_HET;
        else
            callIndex = VAR_HOM;

        if (chip != null || eval != null) {
            //System.out.printf("TEST: %d/%d %s vs. %s%n", truthIndex, callIndex, chip, eval);
            table[truthIndex][callIndex]++;
            truth_totals[truthIndex]++;
            if (callIndex != NO_CALL) calls_totals[callIndex]++;
        }
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        Variation chip = (Variation) tracker.lookup(dbName, null);
	if ( eval != null || chip != null )
	    inc(chip, eval, ref);
        return null;
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

    /**
     * How many overall calls where made that aren't NO_CALLS or UNKNOWNS?
     */
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

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("name                 %s", dbName));
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
        return s;
    }

    private static String cellPercent(int count, int total) {
        StringBuffer sb = new StringBuffer();
        total = Math.max(total, 0);
        sb.append(String.format("%.2f", (100.0 * count) / total));
        sb.append("%");
        return sb.toString();
    }
}
