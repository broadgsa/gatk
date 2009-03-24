package org.broadinstitute.sting.playground.indels;

import org.broadinstitute.sting.utils.Pair;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 21, 2009
 * Time: 4:25:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class ConsensusSequence {
    private List< int[] > coverage; // counts of observations of every nucleotide at each position
    private long referencePos;// (arbitrary) reference position; when adding sequences, their offsets should be wrt this position
    private int startOffset; // offset of the leftmost base of this consensus sequence wrt referencePos; negative=left, positive=right
    private final static int NUMBINS = 4;
    private final static char[] BASES = { 'A','C','G','T' };

    public ConsensusSequence(int refPos) {
        coverage = new ArrayList< int[] >();
        referencePos = refPos;
        startOffset = 0;
    }

    public ConsensusSequence() {
        this(0);
    }

    /** Adds sequence <code>se</code> to the consensus, in the sense that all bases of seq are counted
     * into observed coverage kept by the consensus. The position of the sequence is specified by the
     * <code>offset</code> with respect to the fixed reference position of the consensus (the latter does not
     * have to be consensus start), and if the sequence extends beyound the consensus on either end, the
     * consensus will be extended appropriately to accomodate the full sequence.
     * @param seq nucleotide sequence ('ACGT...')
     * @param offset position of the start of the sequence relative to the fixed reference position of the consensus
     */
    public void addSequence(String seq, int offset) {
        // if sequence starts before than the currently held consensus oes, extend consensus to the left
        if ( offset < startOffset ) {
            coverage.addAll(0,instantiateCoverageList(startOffset-offset));
            startOffset = offset;
        }
        // if the sequence ends beyound the currently held consensus, extend consensus to the right
        if ( offset + seq.length() > startOffset + coverage.size() ) {
            coverage.addAll( instantiateCoverageList(offset+seq.length() - startOffset - coverage.size()) );
        }

        // count bases from the sequence into the coverage
        int posOnConsensus = offset - startOffset;
        for ( int i = 0 ; i < seq.length() ; i++, posOnConsensus++ ) {
            char base = Character.toUpperCase(seq.charAt(i));
            if ( base == 'N') continue;
            coverage.get(posOnConsensus)[baseToInt(base)]++;
        }
    }

    /** Removes sequence <code>seq</code> from the consensus. More exactly, 1 will be subtracted from current
     * observation counts kept by the consensus for each observed base at every position of the sequence. The
     * position of the sequence is specified by the <code>offset</code> with respect to the reference position
     * of the consensus. NOTE: this method is unchecked and does not verify that the sequence being subtracted
     * was indeed previously added to the consensus and/or that the consenus does accomodate full length of
     * the sequence. If it is not the case, the results can be unpredictable or assert failure may occur.
     *
     * @param seq nucleotide sequence ('ACGT...')
     * @param offset position of the start of the sequence relative to the fixed reference position of the consensus
     */
    public void removeSequence(String seq, int offset) {
        assert offset >= startOffset :
                "Attempt to remove from consensus a sequence that starts prior to consenus start";
        assert (offset+seq.length() < startOffset + coverage.size()) :
                "Attempt to remove from consensus a sequence that extends beyond consensus end";
        // subtract sequence bases from the coverage
        int posOnConsensus = offset - startOffset;
        for ( int i = 0 ; i < seq.length() ; i++, posOnConsensus++ ) {
            char base = Character.toUpperCase(seq.charAt(i));
            if ( base == 'N') continue;
            coverage.get(posOnConsensus)[ baseToInt(base) ]--;
        }
    }

    /** Returns offset of the start of consensus sequence with respect to the reference position the
     * consensus is pinned to.
     * @return
     */
    public int getStartOffset() { return startOffset; }

    /** Returns the length (number of bases) of the consensus sequence.
     *
     * @return
     */
    public int length() { return coverage.size(); }

    /** Returns the "distance" (score measuring the agreement) from the currently held consensus sequence to
     * the specified sequence <code>seq</code> starting at position <code>offset</code> wrt consenus reference position.
     * @param seq
     * @param offset
     * @return
     */
    public double distance(String seq, int offset) {
        int posOnConsensus; //  index into the currently held consensus sequence
        int i ; // index into the passed sequence argument
        if ( offset < startOffset ) {
            posOnConsensus = 0;
            i = startOffset - offset;
        } else {
            i = 0 ;
            posOnConsensus = offset - startOffset;
        }
        // stop position on the passed sequence (can be less than sequence length if consensus stops prematurely)
        int stop = Math.min(offset+seq.length(), startOffset+coverage.size() ) - offset;

        for (  ; i < stop ; i++, posOnConsensus++ ) {
            int base = baseToInt(Character.toUpperCase(seq.charAt(posOnConsensus)));
            int [] cov = coverage.get(posOnConsensus);
            int totalcov = cov[0]+cov[1]+cov[2]+cov[3];

        }
        return 0.0;
    }

    /** Returns consensus base at the specified offset wrt the consesus sequence's reference position.
     * Specified offset must be within the span of currently held consensus sequence. Consensus base is the
     * one with the maximum count of observations. If two different nucleotides were observed exactly the
     * same number of times (and that number is greater than the number of observations for othe nucleotides),
     * the "lesser" one, (order being ACGT) will be returned. If coverage at specified position is zero, 'N' will
     * be returned.
     * @param offset
     * @return
     */
    public char baseAt(int offset) {
        assert offset >= startOffset && offset < startOffset + coverage.size() : "Offset out of bounds";
        int [] cov = coverage.get(offset+(int)referencePos);
        int total_cov = cov[0] + cov[1] + cov[2] + cov[3];
        int bmax = 0;
        char base = 'N';
        for ( int z = 0; z < 4 ; z++ ) {
            if ( cov[z] > bmax ) {
                bmax = cov[z];
                base = BASES[z];
            }
        }
        return base;        
    }

    /** Returns consensus base at the specified offset together with its observation count.
     *
     * @param offset
     * @return
     * @see #baseAt(int)
     */
    public Pair<Character,Integer> baseWithCountAt(int offset) {
        assert offset >= startOffset && offset < startOffset + coverage.size() : "Offset out of bounds";
        int [] cov = coverage.get(offset+(int)referencePos);
        int total_cov = cov[0] + cov[1] + cov[2] + cov[3];
        int bmax = 0;
        char base = 'N';
        for ( int z = 0; z < 4 ; z++ ) {
            if ( cov[z] > bmax ) {
                bmax = cov[z];
                base = BASES[z];
            }
        }
        return new Pair<Character,Integer>(base,bmax);
    }

    /** Returns total coverage (all observations regardless of what base what observed) at position
     * specified by offset with respect to the conensus' reference position. offset does not have to be within
     * the bounds of the currently kept consensus sequence, if it falls outside, a 0 will be silently returned.
     * @param offset
     * @return
     */
    public int coverageAt(int offset) {
        if ( offset < startOffset || offset >= startOffset + coverage.size() ) return 0;
        int [] cov = coverage.get(offset+(int)referencePos);
        return cov[0]+cov[1]+cov[2]+cov[3];
    }

    /** Returns consesus sequence as a astring of bases (ACGTN); N will be returned for positions with zero
     * coverage.
     * @return
     */
    public String getSequence() {
        char [] b = new char[coverage.size()];
        for ( int i = 0 ; i < b.length ; i++ ) {
            b[i] = baseAt(i+startOffset);
        }
        return new String(b);
    }

    private List<int[]> instantiateCoverageList(int n) {
        List< int[] > subseq = new ArrayList<int[] >(n);
        for ( int i = 0 ; i < n ; i++ ) subseq.add(new int[NUMBINS]);
        return subseq;
    }

    private int baseToInt(char c) {
        int base;
        switch( Character.toUpperCase(c) ) {
            case 'A': base = 0; break;
            case 'C': base = 1; break;
            case 'G': base = 2; break;
            case 'T': base = 3; break;
            case 'N': base = -1; break;
            default : throw new IllegalArgumentException("Sequence can contain only ACGTN symbols");
        }
        return base;
    }
}
