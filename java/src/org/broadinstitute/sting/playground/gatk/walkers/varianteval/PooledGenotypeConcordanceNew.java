package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.playground.utils.SQuad;

import java.io.*;
import java.util.LinkedList;
import java.util.List;

class HapmapConcordanceTableOld {

    private int[][][] truthTableByName;
    private int[][] poolCallConcordance;
    private int[][] callsByName;
    private int[][] truthByName;
    private int absoluteFalsePositives;
    private int sitesWithFullDataAndNoSNPs;

    private final int REF = 0;
    private final int NOVARIANT = 0;
    private final int HET = 1;
    private final int HOM = 2;
    private final int UNKNOWN = 3;
    private final int MATCH = 1;
    private final int MISMATCH = 2;
    private final int NO_CALL = 3;
    private final int DE_NOVO_CALL = 3; //synonym
    private final int ONE_SNP = 0;
    private final int MULTI_SNP = 1;
    private final int DE_NOVO = 2;
    private final int TRUE_INT = 1;
    private final int FALSE_INT = 0;

    public HapmapConcordanceTableOld( int numberOfPeople ) {
        absoluteFalsePositives = 0;
        sitesWithFullDataAndNoSNPs = 0;
        truthTableByName = new int[numberOfPeople][4][4];
        callsByName = new int[numberOfPeople][4];
        truthByName = new int[numberOfPeople][4];
        poolCallConcordance = new int[3][4];

        for( int i = 0; i < numberOfPeople; i ++ ) {
            for( int j = 0; j < 4; j ++) {
                callsByName[i][j] = 0;
                truthByName[i][j] = 0;
                for( int k = 0; k < 4; k ++ ) {
                    truthTableByName[i][j][k] = 0;
                }
            }
        }

        for( int i = 0; i < 3; i ++) {
            for( int j = 0; j < 4; j ++ ) {
                poolCallConcordance[i][j] = 0;
            }
        }
    }

    public void updateTable(RefMetaDataTracker tracker, char ref, String[] names, Variation eval, GenomeLoc loc, DebugWriter debug) {
        if( ref != 'N' && ref != 'n') {
            int nHapmapVariants = 0;
            int nCorrectCalls = 0;
            int hapmapCoverage = 0;
            for ( int name = 0; name < names.length; name ++ ) {
                Variation chip = (Variation) tracker.lookup(names[name], null);
                SQuad<Integer> variantAndCall = update(eval, name, chip, Character.toUpperCase(ref), debug);
                nHapmapVariants += variantAndCall.getFirst();
                nCorrectCalls += variantAndCall.getSecond();
                hapmapCoverage += variantAndCall.getThird();
                updateFalsePositives(hapmapCoverage,nHapmapVariants,variantAndCall.getFourth());
                debug.printLocusInfo(eval, names[name], ref, chip, loc, variantAndCall.getFourth());
            }
            updatePoolCallConcordance(nHapmapVariants,nCorrectCalls);
            debug.printVerbose(String.format("%s   Variants In Pool: %d     Correctly Called: %d%n", loc.toString(), nHapmapVariants, nCorrectCalls));
        } else {
            // todo -- what do we want to do if ref is N ?? Nothing at the moment . . .
        }
    }

    public void updatePoolCallConcordance(int nHapmapVariants, int nCorrectCalls) {
        int hapmapAllelesIndex, callConcordanceIndex;
        if( nHapmapVariants == 0) {
            hapmapAllelesIndex = DE_NOVO;
            if ( nCorrectCalls == 0 ) {
                callConcordanceIndex = MATCH;
            } else {
                callConcordanceIndex = MISMATCH;
            }
        } else if ( nHapmapVariants == 1 ) {
            hapmapAllelesIndex = ONE_SNP;
            if ( nCorrectCalls > 0 ) {
                callConcordanceIndex = MATCH;
            } else {
                callConcordanceIndex = MISMATCH;
            }
        } else {
            hapmapAllelesIndex = MULTI_SNP;
            if ( nCorrectCalls > 0 ) { // get any of them it's good. We assume variants all have the same nonref base.
                callConcordanceIndex = MATCH;
            } else {
                callConcordanceIndex = MISMATCH;
            }
        }

        poolCallConcordance[hapmapAllelesIndex][callConcordanceIndex]++;
    }

    public void updateFalsePositives(int coverage, int nHapmapVariants, int callConcordanceIndex) {
        if( coverage == truthByName.length ) {
            if( nHapmapVariants == 0  ) {
                if( callConcordanceIndex != NO_CALL ) {
                    if (callConcordanceIndex != NOVARIANT ) {
                        absoluteFalsePositives++;
                        sitesWithFullDataAndNoSNPs++;
                    } else {
                        sitesWithFullDataAndNoSNPs++;
                    }
                }
            }
        }
    }

    public SQuad<Integer> update(Variation eval, int name, Variation chip, char ref, DebugWriter debug) {
        int truthIndex, callIndex, isVariant, isCorrectCallAndThereIsVariant, hapmapCoverage;
        if( chip == null ) {
            truthIndex = UNKNOWN;
            isVariant = FALSE_INT;
            hapmapCoverage = FALSE_INT;
        } else if ( chip.getAlternateBases().charAt(0) == ref && chip.getAlternateBases().charAt(1) == ref) {
            truthIndex = REF;
            isVariant = FALSE_INT;
            hapmapCoverage = TRUE_INT;
        } else if (chip.getAlternateBases().charAt(0) != ref && chip.getAlternateBases().charAt(1) != ref){
            truthIndex = HOM;
            isVariant = TRUE_INT;
            hapmapCoverage = TRUE_INT;
        } else {
            truthIndex = HET;
            isVariant = TRUE_INT;
            hapmapCoverage = TRUE_INT;
        }

        if( eval == null ) {
            callIndex = NO_CALL;
            isCorrectCallAndThereIsVariant = FALSE_INT;
        } else if ( chip == null ) {
            callIndex = DE_NOVO_CALL;
            isCorrectCallAndThereIsVariant = FALSE_INT;
        } else if ( eval.getAlternateBases().charAt(0) == ref && eval.getAlternateBases().charAt(1) == ref ) {
            callIndex = NOVARIANT;
            isCorrectCallAndThereIsVariant = FALSE_INT;
        } else if ( callWrongBase(eval, chip, ref, debug) ) {
            callIndex = MISMATCH;
            isCorrectCallAndThereIsVariant = FALSE_INT;
        } else {
            callIndex = MATCH;
            isCorrectCallAndThereIsVariant = TRUE_INT;
        }

        truthByName[name][truthIndex] ++;
        callsByName[name][callIndex] ++;
        truthTableByName[name][truthIndex][callIndex] ++;

        debug.printDebug(String.format("Person: %d    Truth Index:   %d    Call Index: %d%n", name, truthIndex, callIndex));

        return new SQuad<Integer>(isVariant, isCorrectCallAndThereIsVariant, hapmapCoverage, callIndex);
    }

    public boolean callWrongBase(Variation eval, Variation chip, char ref, DebugWriter debug) {
        // eval and chip guaranteed to be non-null
        char evalRef;
        char evalSNP;

        if ( eval.getAlternateBases().charAt(0) == ref ) {
            evalRef = eval.getAlternateBases().charAt(0);
            evalSNP = eval.getAlternateBases().charAt(1);
        } else {
            evalRef = eval.getAlternateBases().charAt(1);
            evalSNP = eval.getAlternateBases().charAt(0);
        }

        debug.printDebug(String.format("CallWrongBase Check: ref is %s      evalRef is %s      evalSNP is %s%n", ref, evalRef, evalSNP));

        return (evalSNP != chip.getAlternateBases().charAt(0) && evalSNP != chip.getAlternateBases().charAt(1));
    }

    public List<String> standardOutput() {
        LinkedList<String> outLines = new LinkedList<String>();
        int numSingleSNPSites = poolCallConcordance[ONE_SNP][MATCH] + poolCallConcordance[ONE_SNP][MISMATCH];
        int numMultiSNPSites = poolCallConcordance[MULTI_SNP][MATCH] + poolCallConcordance[MULTI_SNP][MISMATCH];
        int numSNPSites = numSingleSNPSites + numMultiSNPSites;
        int numCorrect = poolCallConcordance[ONE_SNP][MATCH] + poolCallConcordance[MULTI_SNP][MATCH];
        int numIncorrect = poolCallConcordance[ONE_SNP][MISMATCH] + poolCallConcordance[MULTI_SNP][MISMATCH];
        outLines.add(String.format("Number of Hapmap SNP Sites:  %d", numSNPSites));
        outLines.add(String.format("\tNumber correctly called:   %d  (%f)", numCorrect, ((double) numCorrect)/numSNPSites));
        outLines.add(String.format("\tNumber incorrectly called: %d  (%f)", numIncorrect, ((double) numIncorrect)/numSNPSites));
        outLines.add(String.format("\tSites with variants in One hapmap individual: %d", numSingleSNPSites));
        outLines.add(String.format("\t\tNumber called correctly:   %d   (%f)",poolCallConcordance[ONE_SNP][MATCH], ((double)poolCallConcordance[ONE_SNP][MATCH]/numSingleSNPSites)));
        outLines.add(String.format("\t\tNumber called incorrectly: %d   (%f)", poolCallConcordance[ONE_SNP][MISMATCH], ((double)poolCallConcordance[ONE_SNP][MISMATCH])/numSingleSNPSites));
        outLines.add(String.format("\tSites with variants in multiple hapmap individuals: %d", numMultiSNPSites));
        outLines.add(String.format("\t\tNumber called correctly:   %d   (%f)", poolCallConcordance[MULTI_SNP][MATCH], ((double)poolCallConcordance[MULTI_SNP][MATCH])/numMultiSNPSites));
        outLines.add(String.format("\t\tNumber called incorrectly: %d   (%f)", poolCallConcordance[MULTI_SNP][MISMATCH], ((double)poolCallConcordance[MULTI_SNP][MISMATCH])/numMultiSNPSites));
        outLines.add(String.format("\tNumber of called sites with homozygous reference chip data on all individuals: %d", sitesWithFullDataAndNoSNPs));
        outLines.add(String.format("\t\tNumber incorrectly called SNPs: %d   (%f)", absoluteFalsePositives, ((double)absoluteFalsePositives)/sitesWithFullDataAndNoSNPs));
        return outLines;
    }

    public List<String> verboseOutput(String[] names) {
        LinkedList<String> outLines = (LinkedList<String>) this.standardOutput();
        for( int nameOffset = 0; nameOffset < names.length; nameOffset ++ ) {
            outLines.add(String.format("Concordance for Hapmap individual %s:", names[nameOffset]));
            int numSNPSites = truthByName[nameOffset][HET] + truthByName[nameOffset][HOM];
            int numCallSites = callsByName[nameOffset][MATCH] + callsByName[nameOffset][MISMATCH];
            outLines.add(String.format("hi %s", "hello"));
        }
        
        return outLines;
    }

}

class DebugWriter {

    private final int SUPPRESS = -1;
    private final int STANDARD = 0;
    private final int INFO = 1;
    private final int VERBOSE = 2;
    private final int DEBUG = 3;
    private final int ALL = 4;

    private PrintWriter writer;
    private int verbosity;

    public DebugWriter() {
        writer = null;
        verbosity = SUPPRESS;
    }

    public DebugWriter(String file, int verbosity) {
        try {
            writer = new PrintWriter(file);
        } catch ( FileNotFoundException e) {
            String errMsg = "Debug file at "+file+" was not found and could not be created. Likely a directory.";
            throw new StingException(errMsg, e);
        }

        this.verbosity = verbosity;
    }

    public void print(String s) {
        if( verbosity >= STANDARD ) {
            writer.print(s);
        }
    }

    public void printInfo(String s) {
        if ( verbosity == INFO || verbosity == ALL) {
            this.print(s);
        }
    }

    public void printVerbose(String s) {
        if ( verbosity == VERBOSE || verbosity == ALL ) {
            this.print(s);
        }
    }

    public void printDebug(String s) {
        if ( verbosity == DEBUG || verbosity == ALL ) {
            this.print(s);
        }
    }

    public boolean isDebug() {
        return verbosity == DEBUG;
    }

    public boolean isVerbose() {
        return verbosity == VERBOSE;
    }

    public boolean isInfo() {
        return verbosity == INFO;
    }

    public int verbosity() {
        return verbosity;
    }

    public void printLocusInfo(Variation eval, String name, char ref, Variation chip, GenomeLoc loc, int callIndex) {
        if(verbosity == INFO || verbosity == ALL) {
            String callStr = (eval == null) ? "NoCall" : eval.getAlternateBases();
            String varStr = (chip == null) ? "NoChip" : chip.getAlternateBases();
            this.print(String.format("%s %s   Call: %s     Chip: %s     Ref: %s   CallIndex: %d%n", name, loc, callStr, varStr, ref, callIndex));
        }
    }

}