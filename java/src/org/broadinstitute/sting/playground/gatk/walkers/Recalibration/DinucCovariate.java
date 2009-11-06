package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 3, 2009
 */
public class DinucCovariate implements Covariate {

    public static ArrayList<String> BASES;

    public DinucCovariate() { // empty constructor is required by CovariateCounterWalker
        BASES = new ArrayList<String>();
        BASES.add("A");
        BASES.add("G");
        BASES.add("C");
        BASES.add("T");
        BASES.add("a");
        BASES.add("g");
        BASES.add("c");
        BASES.add("t");
    }
    
    public Comparable<?> getValue(SAMRecord read, int offset, char[] refBases) {
        byte[] bases = read.getReadBases();
        char base = (char)bases[offset];
        char prevBase = (char)bases[offset - 1];
        if( read.getReadNegativeStrandFlag() ) {
            base = BaseUtils.simpleComplement(base);
            prevBase = BaseUtils.simpleComplement( (char)bases[offset + 1] );
        }
        // Check if bad base, probably an 'N'
        if( !BASES.contains( String.format( "%c", prevBase ) ) || !BASES.contains( String.format( "%c", base) ) ) {
            return null; // CovariateCounterWalker and TableRecalibrationWalker will recognize that null means skip this particular location in the read
        } else {
            return String.format("%c%c", prevBase, base);
        }
    }
    
    public Comparable<?> getValue(String str) {
    	return str;
    }

    public String toString() {
        return "Dinucleotide";
    }
}