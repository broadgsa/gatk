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

    private static final int CALL_REF = 0;
    private static final int CALL_VAR_HET = 1;
    private static final int CALL_VAR_HOM = 2;
    private static final int CALL_NO_CONF = 3;
    private static final int UNCALLABLE = 4;

    private int[][] table = new int[4][5];

    public GenotypeConcordance(final String name) {
        super("genotype_concordance");
        dbName = name;
        for ( int i = 0; i < 4; i++ ) {
            for ( int j = 0; j < 5; j++ ) {
                table[i][j] = 0;
            }
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
            callIndex = UNCALLABLE;
        else if ( eval.isReference() && Utils.countOccurrences(ref, eval.getGenotype().get(0)) == eval.getGenotype().get(0).length() )
            callIndex = CALL_REF;
        else if ( isHet(eval) )
            callIndex = CALL_VAR_HET;
        else
            callIndex = CALL_VAR_HOM;

        if ( chip != null || eval != null ) {
            System.out.printf("TESTING ME: %d/%d %s vs. %s%n", truthIndex, callIndex, chip, eval);

            table[truthIndex][callIndex]++;
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
        s.add(String.format("\t\tCALLED_REF\tCALLED_VAR_HET\tCALLED_VAR_HOM\tNO_CONF\tUNCALLABLE"));
        s.add(String.format("IS_REF\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d", table[TRUTH_REF][CALL_REF], table[TRUTH_REF][CALL_VAR_HET], table[TRUTH_REF][CALL_VAR_HOM], table[TRUTH_REF][CALL_NO_CONF], table[TRUTH_REF][UNCALLABLE]));
        s.add(String.format("IS_VAR_HET\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d", table[TRUTH_VAR_HET][CALL_REF], table[TRUTH_VAR_HET][CALL_VAR_HET], table[TRUTH_VAR_HET][CALL_VAR_HOM], table[TRUTH_VAR_HET][CALL_NO_CONF], table[TRUTH_VAR_HET][UNCALLABLE]));
        s.add(String.format("IS_VAR_HOM\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d", table[TRUTH_VAR_HOM][CALL_REF], table[TRUTH_VAR_HOM][CALL_VAR_HET], table[TRUTH_VAR_HOM][CALL_VAR_HOM], table[TRUTH_VAR_HOM][CALL_NO_CONF], table[TRUTH_VAR_HOM][UNCALLABLE]));
        s.add(String.format("UNKNOWN\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d", table[TRUTH_UNKNOWN][CALL_REF], table[TRUTH_UNKNOWN][CALL_VAR_HET], table[TRUTH_UNKNOWN][CALL_VAR_HOM], table[TRUTH_UNKNOWN][CALL_NO_CONF], table[TRUTH_UNKNOWN][UNCALLABLE]));
        return s;
    }

    private static boolean isHet(AllelicVariant var) {
        if ( var instanceof Genotype )
            return ((Genotype)var).isHet();

        List<String> genotype = var.getGenotype();
        if ( genotype.size() < 1 )
            return false;


        boolean het = genotype.get(0).charAt(0) != genotype.get(0).charAt(1);
        System.out.printf("********* Genotype %s is het = %b%n", genotype, het);
        return het;        
    }
}