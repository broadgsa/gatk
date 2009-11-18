package org.broadinstitute.sting.playground.gatk.walkers.varianteval;


import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.Pair;

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
public class ConcordanceTruthTable {

    private static final int TRUE_POSITIVE = 0;
    private static final int TRUE_NEGATIVE = 1;
    private static final int FALSE_POSITIVE = 2;
    private static final int FALSE_NEGATIVE = 3;
    private static final int VARIANT = 1;
    private static final String[] POOL_HEADERS = {"TRUE_POSITIVE","TRUE_NEGATIVE","FALSE_POSITIVE","FALSE_NEGATIVE"};

    private static final int REF = 0;
    private static final int VAR_HET = 1;
    private static final int VAR_HOM = 2;
    private static final int UNKNOWN = 3;
    private static final int NO_CALL = 3;   // synonym
    private static final String[] TRUTH_NAMES = {"IS_REF", "IS_VAR_HET", "IS_VAR_HOM", "UNKNOWN"};
    private static final String[] CALL_NAMES = {"CALLED_REF", "CALLED_VAR_HET", "CALLED_VAR_HOM", "NO_CALL"};

    private String name = null;
    private boolean singleSampleMode;

    private int[][] table;
    private int[] truth_totals;
    private int[] calls_totals;


    public ConcordanceTruthTable(String name) {
        // there's a specific sample associated with this truth table
        this.name = name;
        singleSampleMode = true;

        table = new int[4][4];
        truth_totals = new int[4];
        calls_totals = new int[4];
        for (int i = 0; i < 4; i++) {
            truth_totals[i] = 0;
            calls_totals[i] = 0;
            for (int j = 0; j < 4; j++)
                table[i][j] = 0;
        }
    }

    public ConcordanceTruthTable(int nSamples) {
        // there's no specific sample associated with this truth table
        singleSampleMode = false;
        name = "pooled_concordance";
        truth_totals = new int[4];
        calls_totals = new int[4];
        for (int i = 0; i < 4; i++) {
            truth_totals[i] = 0;
            calls_totals[i] = 0;
        }

        initializeFrequencyTable(nSamples);
    }

    private void initializeFrequencyTable( int numChips ) {
        // System.out.println("Frequency Table for Pooled Concordance initialized with number of chips = "+numChips);
        table = new int[numChips*2][4];
        for (int i = 0; i < 4; i++) {
            for ( int freq = 0; freq < 2*numChips; freq ++ ) {
                table[freq][i] = 0;
            }
        }

        // System.out.println("Table Size: "+table.length+" by "+table[1].length);
    }

    public void addEntry(List<Pair<Genotype, Genotype>> chipEvals, Variation eval, char ref) {

        // if the table represents a single sample, then we can calculate genotype stats
        if ( singleSampleMode ) {
            for ( Pair<Genotype, Genotype> chipEval : chipEvals ) {

                Genotype chipG = chipEval.first;
                Genotype evalG = chipEval.second;

                if (chipG == null && evalG == null)
                    continue;

                int truthType = getGenotype(chipG, ref);
                int callType = getGenotype(evalG, ref);

                //System.out.printf("TEST: %d/%d %s vs. %s%n", truthIndex, callIndex, chip, eval);
                addGenotypeEntry(truthType, callType);
            }
        } else { // if we cannot associate tables with individuals, then we are working in a pooled context
            // first we need to expand our tables to include frequency information

            Pair<Genotype,Integer> poolVariant = getPooledAlleleFrequency(chipEvals, ref);

            int truthType = getGenotype(poolVariant.getFirst(),ref); // convenience method; now to interpret

            if (truthType == VAR_HET || truthType == VAR_HOM) {
                truthType = VARIANT;
            }

            int callType = getCallIndex(eval,ref);

            addFrequencyEntry( truthType,callType,poolVariant.getSecond() );

        }

        // TODO -- implement me for pooled mode with frequency stats
        // TODO -- You'll want to use eval and the chips from chipEvals (these are the first members of the pair)
        // TODO -- You'll also need to declare (and initialize) the relevant data arrays for the data
        // TODO -- Indexes like TRUE_POSITIVE are defined above for you
    }

    public Pair<Genotype,Integer> getPooledAlleleFrequency( List<Pair<Genotype,Genotype>> chips, char ref) {
        // TODO -- this is actually just a note that I wanted to appear in blue. This method explicitly uses
        // TODO --- the assumption that tri-allelic sites do not really exist, and that if they do the
        // TODO --- site will be marked as such by an 'N' in the reference, so we will not get to this point.
        Genotype variant = null;
        int frequency = 0;
        if ( chips != null ) {
            for ( Pair<Genotype,Genotype> chip : chips ) {
                Genotype c = chip.getFirst();
                if ( c != null ) {
                    if ( c.isVariant(ref) ) {
                        if ( c.isHet() ) {
                            frequency ++;
                        } else { // c is hom
                            frequency += 2;
                        }
                        variant = c;
                    }
                }
            }
        }

        return new Pair<Genotype, Integer>(variant, frequency);

    }

    private void addFrequencyEntry ( int truthIndex, int callIndex, int numTrueSupportingAlleles ) {
        calls_totals[callIndex] ++;
        truth_totals[truthIndex] ++;

        if ( truthIndex == REF && callIndex == REF ) {
            // true negative
            table[numTrueSupportingAlleles][TRUE_NEGATIVE] ++;
            // sanity check - there should never be an entry in
            // [*][TRUE_NEGATIVE] for * > 0
        } else if ( truthIndex == REF && callIndex == VARIANT ) {
            // false positive
            table[numTrueSupportingAlleles][FALSE_POSITIVE] ++;
        } else if ( truthIndex == VARIANT && (callIndex == NO_CALL || callIndex == REF) ) {
            // false negative
            table[numTrueSupportingAlleles][FALSE_NEGATIVE]++;
        } else if ( truthIndex == VARIANT && callIndex == VARIANT ) {
            // true positive
            table[numTrueSupportingAlleles][TRUE_POSITIVE]++;
        } else {
            // something else is going on; wonky site or something. Don't do anything to the table.
        }
    }

    private static int getCallIndex(Variation eval, char ref) {
        int index;

        if ( eval == null ) {
            index = NO_CALL;
        } else if ( ! eval.isSNP() ) {
            index = REF;
        } else {
            index = VARIANT;
        }

        return index;
    }

    private static int getGenotype(Genotype g, char ref) {
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

    private void addGenotypeEntry(int truthIndex, int callIndex) {
        table[truthIndex][callIndex]++;
        truth_totals[truthIndex]++;
        calls_totals[callIndex]++;
    }

    public void addAllStats(List<String> s) {
        if ( singleSampleMode )
            addGenotypeStats(s);
        else
            addFrequencyStats(s);
    }

    private void addFrequencyStats(List<String> s) {

        // TODO -- implement me for pooled mode with frequency stats
        s.add(String.format("name                 %s",name));
        s.add(String.format("TRUTH_ALLELE_FREQUENCY\tERROR_OR_TRUTH_TYPE\tTOTAL\tAS_PRCT_OF_TOTAL_CALLS\tAS_PRCT_OF_CALLS_AT_FREQUENCY"));

        for ( int af = 0; af < table.length; af ++ ) {
            for ( int errorIndex = 0; errorIndex < 4; errorIndex ++ ) {
                StringBuffer sb = new StringBuffer();
                sb.append(String.format("%f ", ((double) af)/ table.length));
                sb.append(String.format("%s ",POOL_HEADERS[errorIndex]));
                sb.append(String.format("%d ", table[af][errorIndex]));
                sb.append(String.format("%s ", percentOfTotal(table,af,errorIndex)));
                sb.append(String.format("%s ", marginalPercent(table[af],errorIndex)));
                s.add(sb.toString());
            }
        }

    }

    private void addGenotypeStats(List<String> s) {
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

    private static String cellPercent(int count, int total) {
        StringBuffer sb = new StringBuffer();
        total = Math.max(total, 0);
        sb.append(String.format("%.2f", (100.0 * count) / total));
        sb.append("%");
        return sb.toString();
    }

    private static String percentOfTotal(int [][] errorMatrix, int alleleFrequency, int errorIndex ) {
        int count = 0;
        for ( int i = 0; i < errorMatrix.length; i ++ ) {
            for ( int j = 0; j < 4; j ++ ) {
                count += errorMatrix[i][j];
            }
        }

        return cellPercent(errorMatrix[alleleFrequency][errorIndex], count );
    }

    private static String marginalPercent ( int[] errors, int errorIndex ) {
        int count = 0;
        for ( int i = 0; i < 4; i ++ ) {
            count += errors[i];
        }

        return cellPercent(errors[errorIndex], count);
    }
}
