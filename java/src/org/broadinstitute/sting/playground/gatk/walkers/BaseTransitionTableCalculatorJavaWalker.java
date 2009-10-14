package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;

import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;
import java.io.PrintStream;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 12, 2009
 * Time: 2:43:06 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
public class BaseTransitionTableCalculatorJavaWalker extends LocusWalker<ReferenceContextWindow,BaseTransitionTable>{
@Argument(fullName="usePreviousBases", doc="Use previous bases as part of the calculation, uses the specified number, defaults to 0", required=false)
    int nPreviousBases = 0;
@Argument(fullName="useSecondaryBase",doc="Use the secondary base of a read as part of the calculation", required=false)
    boolean useSecondaryBase = false;
@Argument(fullName="confidentRefThreshold",doc="Set the lod score that defines confidence in ref, defaults to 4", required=false)
    int confidentRefThreshold = 4;
@Argument(fullName="pileupMismatchThreshold",doc="Set the maximum number of mismatches at a locus before choosing not to use it in calculation. Defaults to 1.", required=false)
    int pileupMismatchThreshold = 1;

    private UnifiedGenotyper ug;
    private ReferenceContextWindow refWindow;

    public void initialize() {
        ug = new UnifiedGenotyper();
        ug.initialize();
        refWindow = new ReferenceContextWindow(nPreviousBases);
    }

    public BaseTransitionTable reduceInit() {
        return new BaseTransitionTable(nPreviousBases,useSecondaryBase);
    }

    public ReferenceContextWindow map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref.getBase(),context);
        refWindow.update(ref,pileup,baseIsUsable(tracker,ref,pileup,context));

        return refWindow;
    }

    public BaseTransitionTable reduce( ReferenceContextWindow map, BaseTransitionTable confusionCounts ) {
                                            // todo -- refactor
        if ( map.isValidWindow() ) {
            List<SAMRecord> reads = map.getPileup().getReads();
            List<Integer> offsets = map.getPileup().getOffsets();
            ReferenceContext ref = map.getMiddleReferenceContext();
            String forwardContext = map.getForwardRefString();
            String reverseContext = map.getReverseRefString();
            for ( int r = 0; r < reads.size(); r ++ ) {
                if ( includeRead(reads.get(r) ,offsets.get(r)) ) {
                    confusionCounts.update(createConfusionContext(reads.get(r),offsets.get(r),ref,forwardContext,reverseContext));
                }
            }
        }

        return confusionCounts;
    }

    public void onTraversalDone( BaseTransitionTable table ) {
        out.printf("%s%n", makeHeader() );
        table.print(out);
    }


    public boolean baseIsConfidentRef( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        List<GenotypeCall> calls = ug.map(tracker,ref,context);
        return ( ! calls.get(0).isVariant(ref.getBase()) && calls.get(0).getNegLog10PError() >= confidentRefThreshold );

    }

    public boolean baseIsUsable ( RefMetaDataTracker tracker, ReferenceContext ref, ReadBackedPileup pileup, AlignmentContext context ) {
        return pileupBelowMismatchThreshold(ref,pileup) && baseIsConfidentRef(tracker,ref,context) && refIsNotN(ref);
    }

    public boolean pileupBelowMismatchThreshold( ReferenceContext ref, ReadBackedPileup pileup ) {
        char[] bases = pileup.getBases().toCharArray();
        int mismatches = 0;
        for ( char b : bases ) {
            if ( Character.toUpperCase(b) == ref.getBase() ) {
                mismatches++;
            }
        }

        return mismatches > pileupMismatchThreshold;
    }

    public boolean refIsNotN(ReferenceContext ref) {
        return ref.getBase() != 'N';
    }

    public boolean includeRead ( SAMRecord read, int offset) {
        // todo -- do we want to filter out individual reads?

        return true;
    }

    public boolean readWindowContainsNonBaseCharacters( SAMRecord read, int offset, int posNeg ) {
        byte[] bases = read.getReadBases();
        if ( posNeg > 0 ) {
            for ( int i = offset; i < offset + nPreviousBases; i ++ ) {
                char base = Character.toUpperCase(convertIUPACByteToChar(bases[i]));
                if ( ! ( base == 'A' || base == 'G' || base == 'C' || base == 'T') ) {
                    return true;
                }
            }
            return false;
        } else {
            for ( int i = offset; i > offset - nPreviousBases; i -- ) {
                char base = Character.toUpperCase(convertIUPACByteToChar(bases[i]));
                if ( ! ( base == 'A' || base == 'G' || base == 'C' || base == 'T') ) {
                    return true;
                }
            }
            return false;
        }
    }

    public Pair<String,char[]> createConfusionContext(SAMRecord read, int offset, ReferenceContext ref, String refFor, String refRev) {
        char[] baseChars;
        if ( useSecondaryBase ) {
            baseChars = new char[3];
        } else {
            baseChars = new char[2];
        }

        baseChars[0] = Character.toUpperCase(ref.getBase());
        baseChars[1] = Character.toUpperCase(convertIUPACByteToChar(read.getReadBases()[offset]));

        if ( useSecondaryBase ) {
            baseChars[2]  = Character.toUpperCase(BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex(((byte[]) read.getAttribute("SQ"))[offset])));
        }

        Pair<String,char[]> confusionContext;

        if ( read.getReadNegativeStrandFlag() ) {
            confusionContext = new Pair<String,char[]>(refRev, baseChars);
        } else {
            confusionContext = new Pair<String,char[]>(refFor, baseChars);
        }

        return confusionContext;

    }

    private char convertIUPACByteToChar(byte iupacBase) {
        char compBase;

        switch (iupacBase) {
            case 'A':
            case 'a': compBase = 'A'; break;
            case 'C':
            case 'c': compBase = 'C'; break;
            case 'G':
            case 'g': compBase = 'G'; break;
            case 'T':
            case 't': compBase = 'T'; break;
            default:  compBase = '.'; break;
        }

        return compBase;
    }

    public String makeHeader() {

        String output = "Reference_Base\tObserved_Base";

        if ( useSecondaryBase ) {
            output = output + "\tSecondary_Base";
        }

        for ( int i = 0; i < nPreviousBases; i ++ ) {
            output = String.format("%s\t%s", output, "Previous_Base_"+Integer.toString(i));
        }

        output = String.format("%s\t%s\t%s\t%s", output, "Hash", "N_Observations", "As_Proportion");

        return output;
    }

}


class BaseTransitionTable {

    final int N_BASES = 4;
    final int A_OFFSET = 0;
    final int C_OFFSET = 1;
    final int G_OFFSET = 2;
    final int T_OFFSET = 3;

    private BaseStringHash strHash;
    protected int[][] confusionTable;
    protected boolean useSecondaryBase;

    public BaseTransitionTable( int prevBaseStringLength, boolean useSecondaryBase ) {
        this.useSecondaryBase = useSecondaryBase;

        if ( useSecondaryBase ) {
            strHash = new BaseStringHash(prevBaseStringLength+2);
        } else {
            strHash = new BaseStringHash(prevBaseStringLength+1);
        }

        confusionTable = new int[strHash.hashSize()][N_BASES];
        for ( int i = 0; i < strHash.hashSize(); i ++ ) {
            for ( int j = 0; j < N_BASES; j ++ ) {
                confusionTable[i][j] = 0;
            }
        }
    }

    public void update( Pair<String, char[]> confusionContext ) {
        //TODO -- THIS [done]
        String context = confusionContext.getFirst();
        char[] bases = confusionContext.getSecond();
        // secondary base is a part of the hash
        if ( useSecondaryBase ) {
            context = context + bases[2] + bases[1];
        } else {
            context = context + bases[1];
        }
        confusionTable[strHash.hash(context)][strHash.hash(bases[0])] ++;
    }

    public void print( PrintStream out ) {
        long totalCounts = countOverTable();

        for ( int hash = 0; hash < strHash.maxHash(); hash ++ ) {
            for ( int ref = 0; ref < N_BASES; ref ++ ) {
                char[] contextBases = strHash.inverse(hash).toCharArray();
                String output = Character.toString(strHash.invert(ref));
                for ( int b = contextBases.length-1; b > 0; b ++ ) {
                    output = output + "\t" + Character.toString(contextBases[b]);
                }

                output = String.format("%s\t%f\t%d\t%f",output,hash+(double)ref/4,confusionTable[hash][ref], ((double)confusionTable[hash][ref])/totalCounts);
                out.printf("%s%n",output);
            }
        }
    }

    public long countOverTable() {
        long count = 0l;
        for (int hash = 0; hash < strHash.maxHash(); hash++ ) {
            for ( int ref = 0; ref < N_BASES; ref ++ ) {
                count += confusionTable[hash][ref];
            }
        }

        return count;
    }

}

class BaseStringHash {

    // character-level mappings and inverses for
    // recursive hash

    final int A_HASH = 0;
    final int C_HASH = 1;
    final int G_HASH = 2;
    final int T_HASH = 3;
    final char INV_0 = 'A';
    final char INV_1 = 'C';
    final char INV_2 = 'G';
    final char INV_3 = 'T';

    int stringLength;

    public BaseStringHash( int stringLength ) {
        this.stringLength = stringLength;
    }

    public int maxHash() {
        return (int) Math.round(Math.pow(4,stringLength));
    }

    public int hashSize() {
        return (int) Math.round(Math.pow(4,stringLength+1));
    }

    public int hash( char b ) {
        switch(b) {
            case 'A':
                return A_HASH;
            case 'C':
                return C_HASH;
            case 'G':
                return G_HASH;
            case 'T':
                return T_HASH;
            default:
                throw new StingException("Incorrect base type sent to base hashing function, was: "+b);
        }
    }

    public int hash ( String s ) {
        return recursiveHash(s,0);
    }

    public int recursiveHash( String s, int offset ) {
        if ( offset == s.length() ) {
            return 0;
        } else {
            return (int) Math.round(hash(s.charAt(offset))*Math.pow(4,s.length()-offset)) + recursiveHash(s, offset+1);
        }
    }

    public String inverse( int h ) {
        return recursiveInverse(h, 0);
    }

    public char invert( int h ) {
        switch(h) {
            case 0:
                return INV_0;
            case 1:
                return INV_1;
            case 2:
                return INV_2;
            case 3:
                return INV_3;
            default:
                throw new StingException("Non-primitive base-string hash code sent to invert. Expected [0,1,2,3]; received "+ Integer.toString(h));
        }
    }

    public String recursiveInverse(int h, int k ) {
        if ( h == 0 ) {
            return "";
        } else {
            double quaternary = Math.pow(4,stringLength-k);
            int coef = (int) Math.floor( h/quaternary );
            return Character.toString(invert(coef)) + recursiveInverse(h - (int) Math.floor(quaternary), k+1);
        }
    }
}


class ReferenceContextWindow {

    protected int windowSize;
    protected int nPrevBases;
    protected LinkedList<ReadBackedPileup> prevAlignments;
    protected LinkedList<ReferenceContext> prevRefs;
    protected LinkedList<Boolean> usePrevious;
    protected boolean initialized;

    public ReferenceContextWindow( int nPrevBases ) {
        windowSize = 2*nPrevBases + 1;
        this.nPrevBases = nPrevBases;
        prevAlignments = new LinkedList<ReadBackedPileup>();
        prevRefs = new LinkedList<ReferenceContext>();
        usePrevious = new LinkedList<Boolean>();
        initialized = false;
    }

    public void update( ReferenceContext ref, ReadBackedPileup pileup, boolean useLocus ) {
        if ( ! initialized ) {
            prevAlignments.add(pileup);
            prevRefs.add(ref);
            usePrevious.add(useLocus);
            if ( prevAlignments.size() == windowSize ) {
                initialized = true;
            }
        } else {
            prevAlignments.removeFirst();
            prevRefs.removeFirst();
            usePrevious.removeFirst();
            prevAlignments.add(pileup);
            prevRefs.add(ref);
            usePrevious.add(useLocus);
        }
    }

    public String getReferenceString() {
        String ref = "";
        for ( ReferenceContext c : prevRefs ) {
            ref = ref + c.getBase();
        }

        return ref;
    }

    public String getForwardRefString() {
        String ref = "";
        for ( ReferenceContext c : prevRefs.subList(0,nPrevBases+1) ) {
            ref = ref + c.getBase();
        }

        return ref;
    }

    public String getReverseRefString() { // todo -- make sure we want to flip this done (yes we do)
        String ref = "";
        for ( int base = prevRefs.size()-1; base >= nPrevBases; base -- ) {
            ref = ref + prevRefs.get(base).getBase();
        }

        return ref;
    }

    public ReadBackedPileup getPileup() {
        // because lists are 0-indexed, this returns the alignments
        // to the middle base in the window.
        return prevAlignments.get(nPrevBases);
    }

    public ReferenceContext getMiddleReferenceContext() {
        return prevRefs.get(nPrevBases);
    }

    public boolean isValidWindow() {
        boolean valid;
        if ( ! initialized ) {
            valid = false;
        } else {
            valid = true;
            // check if everything is confident ref
            for ( Boolean b : usePrevious ) {
                if ( !b ) {
                    valid = false;
                    break;
                }
            }
            // if still valid, check distances
            if ( valid ) {
                ListIterator<ReferenceContext> iter = prevRefs.listIterator();
                ReferenceContext prev = iter.next();
                while ( iter.hasNext() ) {
                    ReferenceContext cur = iter.next();

                    if ( cur.getLocus().distance(prev.getLocus()) > 1 ) {
                        valid = false;
                        break;
                    }

                    prev = cur;
                }
            }
        }

        return valid;
    }

    public int getWindowSize() {
        return windowSize;
    }

}