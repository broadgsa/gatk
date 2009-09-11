package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Sep 9, 2009
 * Time: 4:45:21 PM
 * To change this template use File | Settings | File Templates.
 *
public class PooledGenotypeConcordance extends BasicVariantAnalysis implements GenotypeAnalysis {
    private String[][] individualsByPool = null;
    private int nPools;

    private static final int REF = 0;
    private static final int VAR_SINGLE = 1;
    private static final int VAR_MULTI = 2;
    private static final int UNKNOWN = 3;
    private static final int NO_CALL = 3;   // synonym
    private static final String[] TRUTH_NAMES = {"IS_REF", "IS_SINGLE_VARIANT", "IS_MULTIPLE_VARIANT", "UNKNOWN"};
    private static final String[] CALL_NAMES = {"CALLED_REF", "CALLED_SINGLE_VARIANT", "CALLED_MULTIPLE_VARIANT", "NO_CALL"};
    // todo -- consider the resolution: single vs multi enough, or should there be deeper resolution?

    private int[][][] table;
    private int[][] truthTotalsByPool;
    private int[][] callTotalsByPool;

    public PooledGenotypeConcordance(String pathToHapmapPoolFile) {
        super("genotype_concordance");
        BufferedReader fileReader = null;
        try {
        fileReader = new BufferedReader(new FileReader(pathToHapmapPoolFile));
        } catch (IOException e) {
            String errorMsg = "Could not open any file at " + pathToHapmapPoolFile + " Please check to see if the path is accurate.";
            throw new StingException(errorMsg, e);
        }

        generateNameTableFromFile(fileReader);

        truthTotalsByPool = new int[nPools][4];
        callTotalsByPool = new int[nPools][4];
        table = new int[nPools][4][4];
        for(int pool = 0; pool < nPools; pool++) {
            for(int call1 = 0; call1 < 4; call1++) {
                for(int call2 = 0; call2 < 4; call2++)
                {
                    table[pool][call1][call2] = 0;
                }
                callTotalsByPool[pool][call1] = 0;
                truthTotalsByPool[pool][call1] = 0;
            }
        }

    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, AlignmentContext context){
        for(int pool = 0; pool < nPools; pool ++) {
            incorporatePoolCall(eval, tracker, ref, context, pool);
        }

        return null;
    }

    public void incorporatePoolCall(AllelicVariant eval, RefMetaDataTracker tracker, char ref, AlignmentContext context, int pool) {
        for( int nameOffset = 0; nameOffset < individualsByPool[pool].length; nameOffset++) {
            inc(eval, tracker, ref, context, pool, nameOffset);
        }
    }

    public void inc(AllelicVariant eval, RefMetaDataTracker tracker, char ref, AlignmentContext context, int pool, int nameOffset) {
        AllelicVariant chip = (AllelicVariant) tracker.lookup(individualsByPool[pool][nameOffset],null);
        if( (chip != null && !chip.isGenotype()) || ! isCorrectVariantType(eval)) {
            String errMsg = "Trying to compare non-pooled data or non-genotype data."
            throw new StingException(errMsg);
        }

        int truthIndex, callIndex;
        if(chip == null) {
            truthIndex = UNKNOWN;
        } else if(chip.isReference()) // todo -- how do we want to do pooled concordance checking?

    }

    public boolean isCorrectVariantType(AllelicVariant eval) {
        // todo -- this. Check if eval is a TypeOf some ROD class that's the right pooled call output that we
        // todo -- want to deal with. For now this will work

        return eval.isPooled();
    }


     //
     // private methods for reading names into individualsByPool from an external file
     //

    private void generateNameTableFromFile(BufferedReader reader) {
        // first line is a special case
        String firstline;
        if(continueReading(reader)) {
            firstline = readLine(reader);
        } else {
            String errMsg = "First line of the HapmapPoolFile " + reader.toString() + " could not be read. Please check that the path and file properly formatted.";
            throw new StingException(errMsg);
        }

        StringTokenizer tokFirstLine = new StringTokenizer(firstline);
        nPools = tokFirstLine.countTokens();
        LinkedList<String>[] namesByPool = new LinkedList[nPools];

        for( int firstLineIterator = 0; firstLineIterator < nPools; firstLineIterator ++) {
            namesByPool[firstLineIterator] = new LinkedList();
            namesByPool[firstLineIterator].add(tokFirstLine.nextToken());
        }

        while(continueReading(reader)) {
            String line = readLine(reader);
            StringTokenizer tokLine = new StringTokenizer(line);
            int newNames = tokLine.countTokens();
            for(int lineIt = 0; lineIt < newNames; lineIt ++ ) {
                namesByPool[lineIt].add(tokLine.nextToken());
            }
        }

        convertListOfNamesToMatrix(namesByPool);
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

    private void convertListOfNamesToMatrix(LinkedList<String>[] names) {
        // initialize matrix
        for( int pool = 0; pool < nPools; pool ++) {
            individualsByPool[pool] = new String[names[pool].size()];
            individualsByPool[pool] = names[pool].toArray(individualsByPool[pool]);
        }

    }
}
*/