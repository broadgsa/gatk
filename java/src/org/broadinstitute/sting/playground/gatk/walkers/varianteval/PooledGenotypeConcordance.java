package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 9, 2009
 * Time: 4:45:21 PM
 * To change this template use File | Settings | File Templates.
 */

//todo -- onTraversalDone

public class PooledGenotypeConcordance extends BasicVariantAnalysis implements GenotypeAnalysis {
    private String[] individuals = null;

    private static final boolean DEBUG = true;

    private static final int REF = 0;
    private static final int VAR_MATCH = 1;
    private static final int VAR_HET = 1;
    private static final int VAR_MISMATCH = 2;
    private static final int VAR_HOM = 2;
    private static final int UNKNOWN = 3;
    private static final int NO_CALL = 3;   // synonym
    private static final String[] TRUTH_NAMES = {"IS_REF", "IS_SINGLE_VARIANT", "IS_MULTIPLE_VARIANT", "UNKNOWN"};
    private static final String[] CALL_NAMES = {"CALLED_REF", "CALLED_SINGLE_VARIANT", "CALLED_MULTIPLE_VARIANT", "NO_CALL"};

    private int[][][] tableByIndividual;
    private int[][] truthTotalsByIndividual;
    private int[][] callTotalsByIndividual;

    public PooledGenotypeConcordance(String pathToHapmapPoolFile) {
        super("genotype_concordance");
        if(pathToHapmapPoolFile == null) {
            // do nothing
        } else {
            BufferedReader fileReader = null;
            try {
                fileReader = new BufferedReader(new FileReader(pathToHapmapPoolFile));
            } catch (IOException e) {
                String errorMsg = "Could not open any file at " + pathToHapmapPoolFile + " Please check to see if the path is accurate.";
                throw new StingException(errorMsg, e);
            }

            generateNameTableFromFile(fileReader);

            if(DEBUG) {
                printAllToCheck(individuals);
            }
            
            truthTotalsByIndividual = new int[individuals.length][4];
            callTotalsByIndividual = new int[individuals.length][4];
            tableByIndividual = new int[individuals.length][4][4];
            for(int individual = 0; individual < individuals.length; individual ++) {
                for(int call1 = 0; call1 < 4; call1++) {
                    for(int call2 = 0; call2 < 4; call2++)
                    {
                        tableByIndividual[individual][call1][call2] = 0;
                    }
                    callTotalsByIndividual[individual][call1] = 0;
                    truthTotalsByIndividual[individual][call1] = 0;
                }
            }
        }
    }
    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context){
        if(DEBUG && context != null) {
            System.out.println( "Update Called at location: " + context.getLocation().toString() );
        } else if (DEBUG) {
            System.out.println("Update Called, but no alignment context was given.");
        }

        for( int nameOffset = 0; nameOffset < individuals.length; nameOffset++ ) {
            inc(eval, tracker, ref, context, nameOffset);
        }
        return null;
    }


    public void inc(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context, int nameOffset) {
        Variation chip = (Variation) tracker.lookup(individuals[nameOffset],null);
        if ((chip != null && !(chip instanceof VariantBackedByGenotype) || (eval != null && !(eval instanceof VariantBackedByGenotype)))) {
            String errMsg = "Trying to compare non-pooled data or non-genotype data.";
            throw new StingException(errMsg);
        }

        DiploidGenotype g = DiploidGenotype.valueOf(Utils.dupString(ref,2));

        if(DEBUG) {
            System.out.print(incDebugString(nameOffset, chip, eval, ref));
        }

        // todo -- is this how want to do pooled concordance checking?
        int truthIndex, callIndex;
        if( chip == null ) {
            truthIndex = UNKNOWN;
        } else if ( chip.getAlternateBases().equals(g.toString()) ) {
            truthIndex = REF;
        } else if ( chip.getAlternateBases().charAt(0) != chip.getAlternateBases().charAt(1)) {
            truthIndex = VAR_HET;
        } else {
            truthIndex = VAR_HOM;
        }

        if( eval == null ) {
            callIndex = UNKNOWN;
        } else if ( eval.getAlternateBases().equals(g.toString()) ) {
            callIndex = REF;
        } else if ( callWrongBase(eval,chip,ref) ) {
            callIndex = VAR_MISMATCH;
        } else {
            callIndex = VAR_MATCH;
        }

        if( chip != null || eval != null ) {
            tableByIndividual[nameOffset][truthIndex][callIndex]++;
            truthTotalsByIndividual[nameOffset][truthIndex]++;
            if ( callIndex != NO_CALL ) {
                callTotalsByIndividual[nameOffset][callIndex]++;
            }
        }

        if(DEBUG) {
            System.out.printf("TruthIndex: %d  CallIndex: %d%n", truthIndex, callIndex);
        }



    }

    public boolean isCorrectVariantType(Variation eval) {
        // todo -- this. Check if eval is a TypeOf some ROD class that's the right pooled call output that we
        // todo -- want to deal with.

        return true;
    }

    public boolean callWrongBase(Variation eval, Variation chip, char ref) {
        boolean wrongCall;
        if ( chip == null ) {
            wrongCall = true;
        } else if ( chip.getAlternateBases().charAt(0) == chip.getAlternateBases().charAt(1) ) { // homozygous nonref
            if( eval.getAlternateBases().charAt(0) == chip.getAlternateBases().charAt(0) ||
                eval.getAlternateBases().charAt(1) == chip.getAlternateBases().charAt(0) ) {
                wrongCall = false;
            } else {
                wrongCall = true;
            }
        } else if ( chip.getAlternateBases().charAt(0) == ref || chip.getAlternateBases().charAt(1) == ref ) { // het nonref
            wrongCall = ( getNonrefBase(eval,ref) != getNonrefBase(chip,ref) );
        } else { // todo -- what do we really want to do if ref is C, but chip is AG ??
            wrongCall = true;
        }
        return wrongCall;
    }

    public char getNonrefBase(Variation var, char ref) {
        char nonRefBase;
        if( var.getAlternateBases().charAt(0) == ref ) {
            nonRefBase =  var.getAlternateBases().charAt(1);
        } else {
            nonRefBase = var.getAlternateBases().charAt(0);
        }

        return nonRefBase;
    }



     //
     // private methods for reading names into individualsByPool from an external file
     //

    private void generateNameTableFromFile(BufferedReader reader) {
        LinkedList<String> nameList = new LinkedList<String>();

        while(continueReading(reader)) {
            String line = readLine(reader);
            nameList.add(line);
        }

        individuals = nameList.toArray(new String[nameList.size()]);
    }

    private boolean continueReading(BufferedReader reader) {
        boolean continueReading = false;
        try {
            continueReading = reader.ready();
        } catch(IOException e) {
            continueReading = false;
        }
        return continueReading;
    }

    private String readLine(BufferedReader reader) {
        String line;
        try {
            line = reader.readLine();
        } catch( IOException e) {
            String errMsg = "BufferedReader pointing to "+reader.toString()+" was declared ready but no line could be read from it.";
            throw new StingException(errMsg,e);
        }
        return line;
    }

    // code strictly for debug

    private void printAllToCheck(String[] blah) {
        System.out.println("Checking:");
        for( String s : blah) {
            System.out.println(s);
        }
    }

    private String incDebugString(int nameOffset, Variation chip, Variation eval, char ref) {
        String truth;
        if(chip == null) {
            truth = "NoChip";
        } else {
            truth = chip.getAlternateBases();
        }

        String call;
        if(eval == null) {
            call = "NoCall";
        } else {
            call = eval.getAlternateBases();
        }
        return String.format("Person: %s   Ref: %s   Truth: %s   Call: %s%n",
                             individuals[nameOffset], Character.toString(ref), truth,
                             call);
    }

}
