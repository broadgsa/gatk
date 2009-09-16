package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 9, 2009
 * Time: 4:45:21 PM
 * To change this template use File | Settings | File Templates.
 */
    
public class PooledGenotypeConcordance extends BasicVariantAnalysis implements GenotypeAnalysis {

    private PooledConcordanceTable table;
    private String[] hapmapNames;

    public PooledGenotypeConcordance( String pathToPoolFile ) {
        super("Pooled Genotype Concordance");
        if( pathToPoolFile == null ) {
            table = null;
            hapmapNames = null;
        } else {
            generateNameTableFromFile( pathToPoolFile );
            table = new PooledConcordanceTable( hapmapNames.length );
        }
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        table.updateTable(tracker, ref, eval, hapmapNames);
        return null;
    }

    public List<String> done() {
        return table.standardOutput();
    }

    // private methods for reading in names from a file

    private void generateNameTableFromFile(String file) {
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(file));
        } catch( FileNotFoundException e) {
            String errMsg = "Hapmap pool file at "+file+" was not found. Please check filepath.";
            throw new StingException(errMsg, e);
        }

        LinkedList<String> nameList = new LinkedList<String>();

        while(continueReading(reader)) {
            String line = readLine(reader);
            nameList.add(line);
        }

        hapmapNames = nameList.toArray(new String[nameList.size()]);
    }

    private boolean continueReading(BufferedReader reader) {
        boolean continueReading;
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

}

class PooledConcordanceTable {
    private final int CALL_INDECES = 5;
    private final int TRUTH_INDECES = 4;
    private final int NO_CALL = 0;
    private final int REF_CALL = 0; // synonym
    private final int VARIANT_CALL_NONHAPMAP = 1;
    private final int VARIANT_CALL_MATCH = 2;
    private final int VARIANT_CALL = 2; // re-using index
    private final int VARIANT_CALL_MISMATCH = 3;
    private final int VARIANT_CALL_UNKNOWN = 4;
    private final int NO_TRUTH_DATA = 0;
    private final int TRUTH_REF = 1;
    private final int TRUTH_VAR = 2;
    private final int TRUTH_UNKNOWN = 3;

    private int[][] table;
    private int[][][] tableByHMFrequency;
    private int variantCallsAtRefN;
    private int poolSize;

    public PooledConcordanceTable(int poolSize) {
        this.poolSize = poolSize;
        table = new int[TRUTH_INDECES][CALL_INDECES];
        tableByHMFrequency = new int[calculateNumFrequencyIndeces(poolSize)][TRUTH_INDECES][CALL_INDECES];
        variantCallsAtRefN = 0;

        for ( int j = 0; j < TRUTH_INDECES; j ++ ) {
            for ( int k = 0; k < CALL_INDECES; k ++ ) {
                table[j][k] = 0;
                for( int i = 0; i < calculateNumFrequencyIndeces(poolSize); i ++ ) {
                    tableByHMFrequency[i][j][k] = 0;
                }
            }
        }
    }

    public void updateTable(RefMetaDataTracker tracker, char ref, Variation eval, String[] names) {
        int truthIndex, callIndex, frequencyIndex;
        List<Variation> chips = getChips(names, tracker);
        if ( ref == 'N' || ref == 'n' ) {
            variantCallsAtRefN += ( eval == null ) ? 0 : 1;
            truthIndex = NO_TRUTH_DATA; // don't want calls on Ns to factor in to calculation
            callIndex = NO_CALL; // don't want calls on Ns to factor in to calculation
        } else if( chips.isEmpty() ) {
            truthIndex = NO_TRUTH_DATA;
            if ( eval == null ) {
                callIndex = NO_CALL;
            } else {
                callIndex = ( pooledCallIsRef( eval, ref ) ) ? REF_CALL : VARIANT_CALL_NONHAPMAP;
            }
        } else {
            frequencyIndex = calcVariantFrequencyIndex(ref, chips);

            if ( freqIndexToFrequency( frequencyIndex ) != 0 ) {
                truthIndex = TRUTH_VAR;
            } else if ( chips.size() == poolSize ) {
                truthIndex = TRUTH_REF;
            } else {
                truthIndex = TRUTH_UNKNOWN;
            }

            if ( eval == null ) {
                callIndex = NO_CALL;
            } else {
                if ( pooledCallIsRef( eval, ref ) ) {
                    callIndex = REF_CALL;
                } else if (truthIndex == TRUTH_REF ) {
                    callIndex = VARIANT_CALL;
                } else if (truthIndex == TRUTH_UNKNOWN) {
                    callIndex = VARIANT_CALL_UNKNOWN;
                } else if ( mismatchingCalls( eval, chips, ref ) ) {
                    callIndex = VARIANT_CALL_MISMATCH;
                } else {
                    callIndex = VARIANT_CALL_MATCH;
                }
            }
            
            tableByHMFrequency[frequencyIndex][truthIndex][callIndex] ++; // note that this updates only on HAPMAP sites. This is good.
        }

        table[truthIndex][callIndex] ++;
    }

    public List<Variation> getChips(String[] names, RefMetaDataTracker tracker) {
        LinkedList<Variation> chips = new LinkedList<Variation>();
        for ( String name : names ) {
            Variation chip = (Variation) tracker.lookup(name,null);
            if ( chip != null ) {
                chips.add(chip);
            }
        }
        
        return chips;
    }

    public boolean pooledCallIsRef(Variation eval, char ref) {
        // code broken out for easy alteration when we start using pool-specific variations
        return eval.getAlternateBases().equalsIgnoreCase((Utils.dupString(ref,2)));
    }

    public int calculateNumFrequencyIndeces(int poolSize) {
        // code broken out for easy alteration when we start using pool-specific variations
        return 2*(poolSize+1);
    }

    public int calcVariantFrequencyIndex( char ref, List<Variation> evals ) {
        return frequencyToFrequencyIndex( calcVariantFrequency( evals, ref ) );
    }

    public int frequencyToFrequencyIndex( double frequency ) {
        // code broken out for easy alteration when we start using pool-specific variations
        return (int) frequency;
    }

    public double calcVariantFrequency( List<Variation> evals, char ref ) {
        // code broken out for easy alteration when we start using pool-specific variations
        Variation firstEval = evals.get(0);
        double alternateFrequency = 0.0;
        for ( Variation eval : evals ) {
            if ( mismatchingCalls(firstEval, eval, ref) ) {
                // todo -- make this not a StingException but go to the log
                throw new StingException("Tri-Allelic Position "+eval.getAlternateBases()+"/"+firstEval.getAlternateBases() + " Ref: "+ ref + " not supported");
            } else {
                alternateFrequency += calledVariantFrequency(eval,ref);
            }
        }
        return alternateFrequency;
    }

    public boolean mismatchingCalls(Variation var1, List<Variation> vars, char ref) {
        ListIterator varIter = vars.listIterator();
        try {
            while( pooledCallIsRef( (Variation) varIter.next(), ref ) ) {
                // don't do squat!
            }
        } catch (NoSuchElementException e) {
            String errMsg = "Comparison data given to mismatchingCalls(Variation, List<Variation>, char) was entirely reference. This is a bug and should never happen.";
            throw new StingException(errMsg, e);
        }

        return mismatchingCalls( var1, (Variation) varIter.previous(), ref);
    }

    public boolean mismatchingCalls(Variation eval, Variation chip, char ref) {
        // eval and chip guaranteed to be non-null
        char chipF = chip.getAlternateBases().charAt(0);
        char chipS = chip.getAlternateBases().charAt(1);
        char evalF = chip.getAlternateBases().charAt(0);
        char evalS = chip.getAlternateBases().charAt(1);
        boolean mismatch;
        if (chipF == ref) {
            if ( chipS == ref ) {
                mismatch = false; // no mismatch against hom ref
            } else {
                mismatch = ( evalF != chipS && evalS != chipS ); // test against het nonref
            }
        } else if ( chipS == ref ) {
            mismatch = ( evalF != chipF && evalS != chipF ); // test against het nonref
        } else {
            if( evalF == ref ) {
                mismatch = ( evalS != chipF && evalS != chipF ); // test against hom nonref
            } else if ( evalS == ref ) {
                mismatch = ( evalF != chipF && evalF != chipF ); // test against hom nonref
            } else {
                // both are hom nonref
                mismatch = ( evalF != chipF );
            }
        }
        return mismatch;
    }

    public double calledVariantFrequency( Variation var, char ref ) {
        // code broken out for easy alteration when we start using pool-specific variations
        String varStr = var.getAlternateBases();
        double freq;
        if ( varStr.charAt(0) != ref && varStr.charAt(1) != ref ) {
            freq = (double) 2;
        } else if ( varStr.charAt(0) == ref && varStr.charAt(1) == ref ) {
            freq = (double) 0;
        } else {
            freq = (double) 1;
        }

        return freq;
    }

    public double freqIndexToFrequency( int freqIndex ) {
        // code broken out for easy alteration when we start using pool-specific variations
        return (double) freqIndex;
    }

    public List<String> standardOutput() {
        LinkedList<String> out = new LinkedList<String>();
        int nSNPCallSites = 0;
        int nHapmapSites = 0;
        int nHapmapSNPSites = 0;
        int nFullHapmapRefSites = 0;
        int nUnknownHapmapSites = 0;

        for( int truthIndex = 0; truthIndex < TRUTH_INDECES; truthIndex ++ ) {
            nSNPCallSites += table[truthIndex][VARIANT_CALL_MATCH] + table[truthIndex][VARIANT_CALL_MISMATCH] + table[truthIndex][VARIANT_CALL_NONHAPMAP];
            nSNPCallSites += table[truthIndex][VARIANT_CALL_UNKNOWN];
        }

        for ( int callIndex = 0; callIndex < CALL_INDECES; callIndex ++ ) {
            nHapmapSites += table[TRUTH_REF][callIndex] + table[TRUTH_VAR][callIndex] + table[TRUTH_UNKNOWN][callIndex];
            nUnknownHapmapSites += table[TRUTH_UNKNOWN][callIndex];
            nHapmapSNPSites += table[TRUTH_VAR][callIndex];
            nFullHapmapRefSites += table[TRUTH_REF][callIndex];
        }

        int nRefsCalledCorrectly = table[TRUTH_REF][NO_CALL];
        int nRefsCalledAsSNP = table[TRUTH_REF][VARIANT_CALL];
        int nSNPsCalledCorrectly = table[TRUTH_VAR][VARIANT_CALL_MATCH];
        int nSNPsCalledIncorrectly = table[TRUTH_VAR][VARIANT_CALL_MISMATCH];
        int nSNPsCalledAsRef = table[TRUTH_VAR][REF_CALL];
        int nSNPsOnHapmapSNP = table[TRUTH_VAR][VARIANT_CALL_MATCH] + table[TRUTH_VAR][VARIANT_CALL_MISMATCH];
        int nSNPsAtNonHapmap = table[NO_TRUTH_DATA][VARIANT_CALL_NONHAPMAP];
        int nSNPsAtHapmapWithRefDataOnSubsetOfIndividuals = table[TRUTH_UNKNOWN][VARIANT_CALL_UNKNOWN];

        out.add(String.format("Total Number of SNP Calls:                 %d", nSNPCallSites));
        out.add(String.format("Total SNP calls on non-Hapmap sites        %d", nSNPsAtNonHapmap));
        out.add(String.format("Number of Hapmap Sites:                    %d", nHapmapSites));
        out.add("Data on Hapmap Reference Sites");
        out.add(String.format("\tSites where all Hapmap chips were ref:     %d", nFullHapmapRefSites));
        out.add(String.format("\t\tReference sites correctly called:    %d   (%f)", nRefsCalledCorrectly, ((double)nRefsCalledCorrectly/nFullHapmapRefSites)));
        out.add(String.format("\t\tReference sites called as variant:   %d   (%f)", nRefsCalledAsSNP, ((double)nRefsCalledAsSNP/nFullHapmapRefSites)));
        out.add(String.format("\tSites where all seen Hapmap chips were ref, but not all Hapmap chips available: %d", nUnknownHapmapSites));
        out.add(String.format("\t\tPutative reference sites called ref: %d   (%f)", table[TRUTH_UNKNOWN][REF_CALL], ((double)table[TRUTH_UNKNOWN][REF_CALL]/nUnknownHapmapSites)));
        out.add(String.format("\t\tPutative reference sites called SNP: %d   (%f)", nSNPsAtHapmapWithRefDataOnSubsetOfIndividuals, ((double)nSNPsAtHapmapWithRefDataOnSubsetOfIndividuals/nUnknownHapmapSites)));
        out.add("Data on Hapmap Variant Sites");
        out.add(String.format("\tNumber of Hapmap SNP Sites:                %d", nHapmapSNPSites));
        out.add(String.format("\t\tSNP sites incorrectly called ref:    %d   (%f)", nSNPsCalledAsRef, ((double)nSNPsCalledAsRef/nHapmapSNPSites)));
        out.add(String.format("\tSNP calls on Hapmap SNP Sites:             %d", nSNPsOnHapmapSNP));
        out.add(String.format("\t\tSNP sites correctly called SNP:      %d   (%f)", nSNPsCalledCorrectly, ((double)nSNPsCalledCorrectly/nSNPsOnHapmapSNP)));
        out.add(String.format("\t\tSNP sites called a different base:   %d   (%f)", nSNPsCalledIncorrectly, ((double)nSNPsCalledIncorrectly/nSNPsOnHapmapSNP)));
        out.add(String.format("Calls on reference N:                       %d", variantCallsAtRefN));

        return out;
    }
}
