package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;
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
public class GenotypeConcordance extends BasicVariantAnalysis implements GenotypeAnalysis {
    private String dbName;

    private static final int TRUTH_REF = 0;
    private static final int TRUTH_VAR_HET = 1;
    private static final int TRUTH_VAR_HOM = 2;
    private static final int TRUTH_UNKNOWN = 3;
    private static final String[] TRUTH_NAMES = {"IS_REF", "IS_VAR_HET", "IS_VAR_HOM", "UNKNOWN"};

    private static final int CALL_REF = 0;
    private static final int CALL_VAR_HET = 1;
    private static final int CALL_VAR_HOM = 2;
    private static final int NO_CALL = 3;
    private static final String[] CALL_NAMES = {"CALLED_REF", "CALLED_VAR_HET", "CALLED_VAR_HOM", "NO_CALL"};

    private int[][] table = new int[4][4];
    private int[] truth_totals = new int[4];
    private int[] calls_totals = new int[4];

    public GenotypeConcordance(final String name) {
        super("genotype_concordance");
        dbName = name;
        for ( int i = 0; i < 4; i++ ) {
            truth_totals[i] = 0;
            calls_totals[i] = 0;
            for ( int j = 0; j < 4; j++ )
                table[i][j] = 0;
        }
    }

    public void inc(AllelicVariant chip, AllelicVariant eval, char ref) {
        if ( (chip != null && !chip.isGenotype()) || (eval != null && !eval.isGenotype()) )
            throw new StingException("Failure: trying to analyze genotypes of non-genotype data");

        int truthIndex, callIndex;
        if ( chip == null )
            truthIndex = TRUTH_UNKNOWN;
        else if ( chip.isReference() && Utils.countOccurrences(ref, chip.getGenotype().get(0)) == chip.getGenotype().get(0).length() )
            truthIndex = TRUTH_REF;
        else if ( isHet(chip) )
            truthIndex = TRUTH_VAR_HET;
        else
            truthIndex = TRUTH_VAR_HOM;

        // todo -- FIXME on countOccurences
        if ( eval == null )
            callIndex = NO_CALL;
        else if ( eval.isReference() && Utils.countOccurrences(ref, eval.getGenotype().get(0)) == eval.getGenotype().get(0).length() )
            callIndex = CALL_REF;
        else if ( isHet(eval) )
            callIndex = CALL_VAR_HET;
        else
            callIndex = CALL_VAR_HOM;

        if ( chip != null || eval != null ) {
            //System.out.printf("TEST: %d/%d %s vs. %s%n", truthIndex, callIndex, chip, eval);
            table[truthIndex][callIndex]++;
            truth_totals[truthIndex]++;
            calls_totals[callIndex]++;
        }
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        AllelicVariant chip = (AllelicVariant)tracker.lookup(dbName, null);
        inc(chip, eval, ref);
        return null;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("name                 %s", dbName));
        s.add(String.format("\t\tCALLED_REF\tCALLED_VAR_HET\tCALLED_VAR_HOM\tNO_CALL\t\t\tTOTALS"));
        for (int i=0; i < 4; i++) {
            StringBuffer sb = new StringBuffer();
            sb.append(TRUTH_NAMES[i] + "\t");
            for (int j=0; j < 4; j++)
                sb.append(table[i][j] +" (" + cellPercent(table[i][j], truth_totals[i]) + ")\t\t");
            sb.append(truth_totals[i]);
            s.add(sb.toString());
        }
        s.add("\n");
        s.add(String.format("\t\tCALLED_REF\tCALLED_VAR_HET\tCALLED_VAR_HOM\tNO_CALL"));
        for (int i=0; i < 4; i++) {
            StringBuffer sb = new StringBuffer();
            sb.append(TRUTH_NAMES[i] + "\t");
            for (int j=0; j < 4; j++)
                sb.append(table[i][j] + " (" + cellPercent(table[i][j], calls_totals[j]) + ")\t\t");
            s.add(sb.toString());
        }
        s.add(String.format("TOTALS\t%d\t\t%d\t\t%d\t\t%d", calls_totals[CALL_REF], calls_totals[CALL_VAR_HET], calls_totals[CALL_VAR_HOM], calls_totals[NO_CALL]));
        s.add("\n");
        for (int i=0; i < 4; i++) {
            for (int j=0; j < 4; j++) {
                s.add(TRUTH_NAMES[i]+"_"+CALL_NAMES[j]+"_COUNT "+table[i][j]);
                s.add(TRUTH_NAMES[i]+"_"+CALL_NAMES[j]+"_PERCENT_OF_TRUTH "+cellPercent(table[i][j], truth_totals[i]));
                s.add(TRUTH_NAMES[i]+"_"+CALL_NAMES[j]+"_PERCENT_OF_CALLS "+cellPercent(table[i][j], calls_totals[j]));
            }
        }
        return s;
    }

    private static String cellPercent(int count, int total) {
        StringBuffer sb = new StringBuffer();
        if ( total == 0 )
            sb.append(0);
        else
            sb.append(100*count/total);
        sb.append("%");
        return sb.toString();
    }

    private static boolean isHet(AllelicVariant var) {
        if ( var instanceof Genotype )
            return ((Genotype)var).isHet();

        List<String> genotype = var.getGenotype();
        if ( genotype.size() < 1 )
            return false;

        return genotype.get(0).charAt(0) != genotype.get(0).charAt(1);
    }
}