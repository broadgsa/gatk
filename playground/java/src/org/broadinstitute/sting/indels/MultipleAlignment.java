package org.broadinstitute.sting.indels;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;


public class MultipleAlignment implements Iterable<Integer>  {
	private static final int IMPOSSIBLE = 1000000000;
	private Map<Integer,Integer> index; 
	private List<String> seqs;
	private List<Integer> ext_ids;
	private List<Integer>  alignment_offsets; // offset of seqs[i] w/respect to seqs[0]
	private int best_mm; // mismatch count
	private int next_mm; // next-best mismatch count
		
	public MultipleAlignment() {
		index = new HashMap<Integer,Integer>();
		seqs = new ArrayList<String>();
		alignment_offsets = new ArrayList<Integer>();
		ext_ids = new ArrayList<Integer>();
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
	 */
	public void add( String seq, int i ) throws IllegalStateException {
		if ( size() != 0 ) throw new IllegalStateException("Single sequence can be added to an empty pile only");
		index.put(i,0);
		ext_ids.add(i);
		seqs.add(seq);
		alignment_offsets.add(0);
	}
	
	public void add( PairwiseAlignment a) {
		if ( a.id1() == -1 || a.id2() == -1 ) throw new IllegalArgumentException("Attempt to add pairwise alignemnt with sequence ids not properly set");
		add(a,a.id1(),a.id2());
	}
	
	/** Adds pair of aligned sequences to the pile, with the external ids of the first and second sequences being i and j,
	 * respectively. Pairwise alignment can be always added to an empty pile. If the pile is non-empty and either both or 
	 * none of the specified ids are already in the pile, an IllegalStateException will be thrown
	 * @param a
	 * @param i
	 * @param j
	 */
	public void add( PairwiseAlignment a, int i, int j ) throws IllegalStateException {
		if ( seqs.size() == 0 ) {
			index.put(i,0);
			ext_ids.add(i);
			seqs.add(a.getSequence1());
			index.put(j,1);
			ext_ids.add(j);
			seqs.add(a.getSequence2());
			alignment_offsets.add(0);
			alignment_offsets.add(a.getBestOffset2wrt1());
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
		
		if ( second == null ) {
			index.put(j,seqs.size());
			seqs.add(a.getSequence2());
            ext_ids.add(j);
			alignment_offsets.add( a.getBestOffset2wrt1() + alignment_offsets.get( first ));
		} else {
			// first = null
			index.put(i,seqs.size());
			seqs.add(a.getSequence1());
            ext_ids.add(i);
			alignment_offsets.add( -a.getBestOffset2wrt1() + alignment_offsets.get( second ));
		}
	}
	
	/** Returns sequence associated with the specified external id, or null if sequence with this external id is not found in the pile
	 * 
	 * @param id query id
	 * @return sequence for specified id or null
	 */
	public String getSequenceById(int id) {
		if ( ! contains(id)) return null;
		return seqs.get(index.get(id));
	}
	
	/** Returns offset relative to the first sequence in the pile for sequence associated with the specified external id
	 * 
	 * @param id query id
	 * @return offset for sequence with specified id
	 */
	public int getOffsetById(int id) {
		//TODO: do something meaningful when id is not in the pile (exception?)
		if ( ! contains(id)) return 0;
		return alignment_offsets.get(index.get(id));
	}

	/** Adds pile of aligned sequences to this pile, stitching them together using specified pairwise alignment
	 * p of the sequences with external ids i and j. One of the indices i, j must be in this pile, and the other in the pile being added,
	 * otherwise an IllegalArgumentException is thrown. The piles can not overlap.
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
			index.put(id,seqs.size());
			ext_ids.add(id);
			seqs.add(a.getSequenceById(id));
			alignment_offsets.add(off2+a.getOffsetById(id));
		}
	}
	
	/** Returns true if the alignment already contains sequence with the specified id 
	 * 
	 * @param id
	 * @return
	 */
	public boolean contains(int id) {
		return index.containsKey(id);
	}
	
	public int countMismatches(int i, int j) {
		return PairwiseAlignment.countMismatches(getSequenceById(i), getSequenceById(j), getOffsetById(j)-getOffsetById(i));
	}

	/** Returns the length of the overlapping region of sequences si and sj  .
	 * 
	 * @return overlap size
	 */
	public int getOverlap(int i, int j) {
		if ( ! contains(i) || ! contains(j)  ) return -1;
		int off = getOffsetById(j) - getOffsetById(i);
		if ( off >= 0 ) {
			return Math.min(getSequenceById(i).length()-off, getSequenceById(j).length());
		} else {
			return Math.min(getSequenceById(j).length()+off, getSequenceById(i).length());
		}
	}
	
	/** Given the two indices, one of which has to be already in the pile, returns the one that is not in the pile.
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
	
	/** Returns a (multiline) string that represents the alignment visually: the sequences are appropriately
	 *  shifted and ready for printout;  
	 */
	public String toString() {
		StringBuilder b = new StringBuilder();
		java.util.Formatter frmt = new java.util.Formatter(b);
		
		if ( seqs.size() == 0 ) return b.toString();
		
		int skip_first = 0;
		for ( int i = 0 ; i < seqs.size() ; i++ ) {
			if ( -alignment_offsets.get(i) > skip_first ) skip_first = -alignment_offsets.get(i);
		}
		
		frmt.format("%3d:", ext_ids.get(0));
		for ( int k = 0 ; k < skip_first ; k++ ) b.append(' ');
		b.append(seqs.get(0));
		b.append('\n');
		
		for ( int i = 1 ; i < seqs.size() ; i++ ) {
			frmt.format("%3d:", ext_ids.get(i));
			int skip = alignment_offsets.get(i) + skip_first ;
			for ( int k = 0 ; k < skip ; k++ ) b.append(' ');
			b.append(seqs.get(i));
			b.append('\n');
		}
//		b.append(best_mm+" mismatches, "+ next_mm + " next best, " + getOverlap() + " overlapping bases, distance=" + distance() + "\n");
		return b.toString();
	}
	
	public int size() { return seqs.size(); }
	
	/** Returns an iterator over the id's of the sequences currently stored in the pile 
	 * 
	 * @return
	 */
	public Iterator<Integer> sequenceIdIterator() { return index.keySet().iterator(); }

	@Override
	public Iterator<Integer> iterator() {
		return sequenceIdIterator();
	}
}
