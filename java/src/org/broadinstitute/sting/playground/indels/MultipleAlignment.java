package org.broadinstitute.sting.playground.indels;

import java.util.*;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;


public class MultipleAlignment implements Iterable<Integer>  {
	private static final int IMPOSSIBLE = 1000000000;
	private Map<Integer,Integer> index;    // maps external id of the sequence onto its index in the pile
	private List<String> seqs;             // sequences, in order they were added
	private List<Integer> ext_ids;         // external ids of the sequences, in order they were added to the pile
	private List<Integer>  alignment_offsets; // offset of seqs[i] w/respect to seqs[0] (i.e. in order the seqs were added)
	private int best_mm; // mismatch count
	private int next_mm; // next-best mismatch count
	private ConsensusSequence consensus;

	public MultipleAlignment() {
		index = new HashMap<Integer,Integer>();
		seqs = new ArrayList<String>();
		alignment_offsets = new ArrayList<Integer>();
		ext_ids = new ArrayList<Integer>();
        consensus = new ConsensusSequence(); // we use reference position 0, e.g. we hook onto the first read in the pile
	}

	public void clear() {
		seqs.clear();
		index.clear();
		alignment_offsets.clear();
		ext_ids.clear();
	}

	/** Adds  single sequence with id set to i. Pile must be empty, or IllegalStateException will be thrown
	 * 
	 * @param seq sequence to add
	 * @param i id of the sequence (can be use later to query the pile)
     * @see #add(String,int,int)
	 */
	public void add( String seq, int i ) throws IllegalStateException {
		if ( size() != 0 ) throw new IllegalStateException("Single sequence can be added to an empty pile only");
        add(seq,i,0);
	}
	
    /** Adds  single sequence with id set to i and places it at the specified offset wrt the first sequence
     * in this pile (i.e. wrt reference position 0).
     *
     * @param seq sequence to add
     * @param i id of the sequence (can be use later to query the pile)
     * @see #add(String,int)
     */
    public void add( String seq, int i, int offset ) throws IllegalStateException {
        index.put(i,index.size());
        ext_ids.add(i);
        seqs.add(seq);
        alignment_offsets.add(offset);
        consensus.addSequence(seq,offset);
    }

	public void add( PairwiseAlignment a) {
		if ( a.id1() == -1 || a.id2() == -1 ) throw new IllegalArgumentException("Attempt to add pairwise alignemnt with sequence ids not properly set");
		add(a,a.id1(),a.id2());
	}
	
	/** Adds pair of aligned sequences to the pile, with the external ids of the first and second sequences being i and j,
	 * respectively. Pairwise alignment can be always added to an empty pile. If the pile is non-empty, exactly
     * one of the sequences held by the pair-wise alignment should be already in the pile; this sequence (and the
     * pairwise alignment itself) will be used to stitch the other sequence to the pile. If either both or
	 * none of the specified ids are already in the pile, an IllegalStateException will be thrown.
	 * @param a
	 * @param i
	 * @param j
	 */
	public void add( PairwiseAlignment a, int i, int j ) throws IllegalStateException {
		if ( seqs.size() == 0 ) {
            add(a.getSequence1(),i,0);
            add(a.getSequence2(),j,a.getBestOffset2wrt1());
			return;
		}
		
		Integer first = index.get(i);
		Integer second = index.get(j);
		
		if ( first != null && second != null ) {
			throw new IllegalStateException("Attempt to add pairwise alignment for two sequences that are already in the pile");
		}

		if ( first == null && second == null ) {
			throw new IllegalStateException("Attempt to add pairwise alignment for two sequences none of which is already in the pile");
		}
		
		if ( second == null ) add(a.getSequence2(),j, a.getBestOffset2wrt1() + alignment_offsets.get( first ) );
		else add(a.getSequence1(),i, -a.getBestOffset2wrt1() + alignment_offsets.get( second ) );
	}

    /** Adds another pile of aligned sequences to this pile, stitching them together using specified pairwise alignment
     * p of the sequences with external ids i and j. One of the indices i, j must be in this pile, and the other in
     * the pile being added, otherwise an IllegalArgumentException is thrown. Sequence id's i and j MUST be the ids
     * of the first and second sequences in the pairwise alignment, in that order. Specified ids override
     * ids, if any, set for the sequences in the pairwise alignment; it is not checked whether the specified and
     * stored ids match. The piles can not overlap.
     */
    public void add(MultipleAlignment a, PairwiseAlignment p, int i, int j) {
        int off2; // offset of the first sequence in pile 'a' wrt the first sequence in this pile
        if ( this.contains(i) ) {
            if ( ! a.contains(j)) throw new IllegalArgumentException("Sequence is not in the pile");
            off2 = getOffsetById(i)+p.getBestOffset2wrt1()-a.getOffsetById(j);
        } else {
            if ( this.contains(j)) {
                if ( ! a.contains(i)) throw new IllegalArgumentException("Sequence is not in the pile");
                off2 = getOffsetById(j)-p.getBestOffset2wrt1()-a.getOffsetById(i);
            } else throw new IllegalArgumentException("Sequence is not in the pile");
        }
        // stitch sequences from a into this pile:
        for ( Integer id : a ) {
            if ( this.contains(id) ) throw new IllegalArgumentException("Attempt to add a pile that shares sequences with the current one");
            add(a.getSequenceById(id),id,off2+a.getOffsetById(id));
        }
    }


    /** Adds another pile of aligned sequences (a) to this pile, stitching them together using specified
     * pairwise alignment p. Sequence ids must be set in the pairwise alignment, and one of those ids
     * must be in this pile, and the other in the pile 'a' being added, otherwise an IllegalArgumentException
     * is thrown. If pairwise alignment does not have sequence ids set, IllegalArgumentException is thrown.
     * The piles can not overlap.
     */
    public void add(MultipleAlignment a, PairwiseAlignment p) {
        if ( p.id1() == -1 || p.id2() == -1 ) throw new IllegalArgumentException("Attempt to add MSA based on pairwise alignemnt with sequence ids not properly set");
        add(a,p,p.id1(),p.id2());
    }

	/** Returns sequence associated with the specified external id, or null if sequence with this external id is
     * not found in the pile
	 * 
	 * @param id query id
	 * @return sequence for specified id or null
	 */
	public String getSequenceById(int id) {
		if ( ! contains(id)) return null;
		return seqs.get(index.get(id));
	}
	
	/** Returns offset relative to the first sequence in the pile for sequence associated with the specified
     * external id. If sequence with specified id is not found in the pile, RuntimeException is thrown.
	 * 
	 * @param id query id
	 * @return offset for sequence with specified id
	 */
	public int getOffsetById(int id) {
		if ( ! contains(id) ) throw new RuntimeException("Specified id is not in the pile");
		return alignment_offsets.get(index.get(id));
	}

    /** Returns external id of the read the offsets of this multiple alignment are based upon (i.e. all the offsets
     * are specified wrt the base read).
     * @return
     */
    public int getBaseReadId() { return ext_ids.get(0); }

    /** Returns offset of the read specified by its external id wrt the start of the consensus sequence in this
     * multiple alignment (consenus sequence is a major vote union of all the reads in this alignment).
     * @param id
     * @return
     */
    public int getOffsetWrtConsensus(int id) {
        return getOffsetById (id)- consensus.getStartOffset();
    }

	/** Returns true if the alignment already contains sequence with the specified id.
	 * 
	 * @param id
	 * @return
	 */
	public boolean contains(int id) {
		return index.containsKey(id);
	}

    /** Returns number of mismatches between sequences i and j (external ids) in the currently held multiple alignment.
     * Will return 0 if sequences do not overlap. Will throw RuntimeException if any of the specified ids is not
     * found in the current pile. 
     * @param i id of the first sequence
     * @param j id of the second sequence
     * @return mismatch count
     *
     * */
	public int countMismatches(int i, int j) {
		return PairwiseAlignment.countMismatches(getSequenceById(i), getSequenceById(j), getOffsetById(j)-getOffsetById(i));
	}

	/** Returns the length of the overlapping region of the two sequences specified by their external ids i and j.
	 * 
	 * @return overlap size
	 */
	public int getOverlap(int i, int j) {
		if ( ! contains(i) || ! contains(j)  ) throw new RuntimeException("Sequence with specified id is not in MSA pile");
		int off = getOffsetById(j) - getOffsetById(i);
        int L;
		if ( off >= 0 ) L = Math.min(getSequenceById(i).length()-off, getSequenceById(j).length());
		else L = Math.min(getSequenceById(j).length()+off, getSequenceById(i).length());
		return ( L < 0 ? 0 : L );
	}
	
	/** Given the two sequence ids, one of which has to be already in the pile, returns the one that is not in the pile.
	 * 
	 * @param i sequence id
	 * @param j sequence id
	 * @return one of the input arguments that is not found in the pile
	 * @throws IllegalArgumentException when either both or none of the specified indices are in the pile
	 */
	public int selectExternal(int i, int j) {
		if ( contains(i) ) {
			if ( contains(j) ) throw new IllegalArgumentException("Can not select external when both indices are in the pile");
			return j;
		} else {
			if ( ! contains(j) ) throw new IllegalArgumentException("Attempt to select external when both indices are not in the pile");
			return i;
		}
	}

    /** Returns a string consisting of n spaces.
     *
     * @param n
     * @return
     */
    private String skipN(int n) {
        StringBuilder b=new StringBuilder();
        for ( int k = 0 ; k < n ; k++ ) b.append(' ');
        return b.toString();
    }

    /** Prints n spaces directly into the specified string builder.
     *
     * @param n
     * @param b
     */
    private void skipN(int n, StringBuilder b) {
        for ( int k = 0 ; k < n ; k++ ) b.append(' ');
    }

	/** Returns a (multiline) string that represents the alignment visually: the sequences are appropriately
	 *  shifted and ready for printout;  
	 */
	public String toString(boolean inorder) {

		StringBuilder b = new StringBuilder();
		java.util.Formatter frmt = new java.util.Formatter(b);
		
		if ( seqs.size() == 0 ) return b.toString();
		
        final int first_offset = -consensus.getStartOffset();

        final int msa_length = consensus.length();
        char[][] consensusString = new char[4][msa_length];

        for ( int i = 0 ; i < msa_length ; i++ ) {

            Pair<Character,Integer> base = consensus.baseWithCountAt(i-first_offset);
            consensusString[3][i] = base.first;
            int mm = consensus.coverageAt(i-first_offset) - base.second;
            if ( mm > 0 ) {
                consensusString[2][i] = '*';
                if ( mm > 9 ) consensusString[0][i] = Character.forDigit(mm/10,10);
                else consensusString[0][i] = ' ';
                consensusString[1][i] = Character.forDigit(mm%10,10);
            } else {
                consensusString[0][i] = consensusString[1][i] = consensusString[2][i] = ' ';
            }
        }

        b.append("    "); b.append(consensusString[0]); b.append('\n');
        b.append("    "); b.append(consensusString[1]); b.append('\n');
        b.append("    "); b.append(consensusString[2]); b.append('\n');
        b.append("    "); b.append(consensusString[3]); b.append('\n');

        Integer[] perm = null;
        if ( inorder ) perm = Utils.SortPermutation(alignment_offsets);
		
		for ( int i = 0 ; i < seqs.size() ; i++ ) {
            int index = (inorder ? perm[i] : i);
			frmt.format("%3d:", ext_ids.get(index));
			skipN(alignment_offsets.get(index)+ first_offset,b);
			b.append(seqs.get(index));
			b.append('\n');
		}
//		b.append(best_mm+" mismatches, "+ next_mm + " next best, " + getOverlap() + " overlapping bases, distance=" + distance() + "\n");
		return b.toString();
	}

    public String getConsensus() {
        return consensus.getSequence();
    }

    public String toString() { return toString(true); }

	public int size() { return seqs.size(); }
	
	/** Returns an iterator over the id's of the sequences currently stored in the pile 
	 * 
	 * @return
	 */
	public Iterator<Integer> sequenceIdIterator() { return index.keySet().iterator(); }

    /** Returns an iterator over external seuqnce ids of the sequences stored in the pile, presenting them in
     * the order of ascending alignment offsets.
     * @return
     */
    public Iterator<Integer> sequenceIdByOffsetIterator() {
        final Integer[] perm = Utils.SortPermutation(alignment_offsets);
        return new Iterator<Integer>() {
            private int i = 0;
            public boolean hasNext() {
                return i < perm.length;
            }
            public Integer next() {
                return ext_ids.get(perm[i++]);
            }
            public void remove() {
                throw new UnsupportedOperationException("remove not supported");
            }
        }   ;

    }

	public Iterator<Integer> iterator() {
		return sequenceIdIterator();
	}


}
