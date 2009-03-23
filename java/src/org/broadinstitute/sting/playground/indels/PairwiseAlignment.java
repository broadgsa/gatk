package org.broadinstitute.sting.playground.indels;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class PairwiseAlignment {
		    private static final int IMPOSSIBLE = 1000000000;
			private String s1;
			private String s2;
			private int i1; // (external) id of the first sequence
			private int i2; // (external) id of the second sequence
			private int alignment_offset; // offset of s2 w/respect to s1
			private int best_mm; // mismatch count
			private int next_mm; // next-best mismatch count
			
			/** Initializes the alignment with pair of sequences (that will be immediately aligned) and
			 * stores their specified external ids id1, id2.
			 * @param is1 first nucleotide sequence (pre-indexed)
			 * @param is2 second nucleotide sequence (pre-indexed)
			 * @param id1 external id of the first sequence
			 * @param id2 external id of the second sequence
			 */
			public PairwiseAlignment(IndexedSequence is1, IndexedSequence is2, int id1, int id2 ) {
				s1 = new String(is1.getSequence());
				s2 = new String(is2.getSequence());
				i1 = id1;
				i2 = id2;
				best_mm = IMPOSSIBLE;
				next_mm = IMPOSSIBLE; 
				align(is1,is2);
			}
			
			/** Initializes the alignment with pair of sequences (that will be immediately aligned) and
			 * sets their external ids to -1. Such un-annotated pairwise alignment can not be added to MultipleAlignment.
			 *
			 */
			public PairwiseAlignment(IndexedSequence is1, IndexedSequence is2) {
				this(is1,is2,-1,-1);
			}
			
			/**
			 * Returns offset of sequence 2 with respect to sequence 1 in the best alignment
			 * @return positive offset if s2 is shifted right (starts later) wrt s1, or negative offset
			 *                 if s2 is shifted left (starts earlier) wrt s1
			 */
			public int getBestOffset2wrt1() { return alignment_offset; }

            /** Returns offset of the sequence j wrt sequence i in the best pairwise alignment found.
             *
             * @param i extrenal id of a sequence, must be one of the sequences kept by this alignment
             * @param j extrenal id of a sequence, must be one of the sequences kept by this alignment
             * @return offset of 2nd arg (j) wrt to the first arg (i)
             */
            public int getBestOffset2wrt1(int i, int j ) {
               if ( i == i1 && j == i2 ) return alignment_offset;
               else if ( i == i2 && j == i1 ) return -alignment_offset;
               throw new RuntimeException("Specified sequence id not found in the alignment");
            }

			public String getSequence1() { return s1; }
			public String getSequence2() { return s2; }
            public String getSequenceById(int i) {
                if ( i == i1 ) return s1;
                else if ( i == i2 ) return s2;
                throw new RuntimeException("Specified sequence id not found in the alignment");
            }
			public int id1() { return i1;}
			public int id2() { return i2;}
			
			/** Returns mismatch count in the best alignment found.
			 * 
			 * @return count of mismatches or impossibly large number of no mismatches were found
			 */
			public int getBestMMCount() { return best_mm; }
			
			/** Returns the number of mismatches in the next-best alignment found
			 * 
			 * @return next-best count of mismatches or impossibly large number if at most one alignment
			 * was ever found (that one would make the best then)
			 */
			public int getNextBestMMCount() { return next_mm; }
			
			/** Returns the length of the overlapping region of sequences s1 and s2 in the best alignment found, or -1 if
			 *   sequences do not align.
			 * 
			 * @return overlap size; can not be smaller than the size of the kmer used in IndexedSequence arguments the
			 * alignment was built from
			 */
			public int getOverlap() {
				if ( ! alignmentExists() ) return -1;
				if ( alignment_offset >= 0 ) {
					return Math.min(s1.length()-alignment_offset, s2.length());
				} else {
					return Math.min(s2.length()+alignment_offset, s1.length());
				}
			}

            public static int getOverlap(String seq1, String seq2, int offset2wrt1) {
                int L ;
                if ( offset2wrt1 >= 0 ) {
                    L = Math.min(seq1.length()-offset2wrt1, seq2.length());
                } else {
                    L = Math.min(seq2.length()+offset2wrt1, seq1.length());
                }
                return ( L < 0 ? 0 : L );
            }
			
			/** Returns true if at least one alignment, no matter how bad, was found between the two sequences
			 * (i.e. the sequences have at least one kmer in common).
			 */
			public boolean alignmentExists() { return best_mm < IMPOSSIBLE; }
			
			public void align(IndexedSequence is1, IndexedSequence is2) {
				
				Set<Integer> offsets = new HashSet<Integer>() ; // possible offsets of s2 wrt s1 as suggested by matching kmers
				for ( Map.Entry<Short,List<Integer>> e : is1 ) { // for each kmer in s1
					List<Integer> kmer_offsets_2 = is2.getOffsets(e.getKey());
					if ( kmer_offsets_2 == null ) continue; // uh-oh, kmer is not found in the other sequence
					for ( Integer i1 : e.getValue() ) {
						for ( Integer i2 : kmer_offsets_2 ) {
							offsets.add(i1-i2); // offset of seq 2 wrt seq1 as suggested by the currently inspected  occurences of the same kmer e.getKey() in both sequences 
						}
					}
				}
				// we have now a collection of distinct s1-s2 offsets seeded by matching kmers.
				// lets extend these kmer matches and count mismatches:
				
				for ( Integer trial_offset : offsets ) {
					int mm_cnt = countMismatches(is1.getSequence(), is2.getSequence(), trial_offset,next_mm); 
					if ( mm_cnt < best_mm ) {
						alignment_offset = trial_offset;
						next_mm = best_mm;
						best_mm = mm_cnt;
					} else {
						if ( mm_cnt < next_mm ) next_mm = mm_cnt;
					}
				}
			}
			
			public static int countMismatches(String seq1, String seq2, int offset2wrt1) {
				int pos1 = ( offset2wrt1 >= 0 ? offset2wrt1 : 0 );
				int pos2 = ( offset2wrt1 >= 0 ? 0 : -offset2wrt1 );  
				int cnt = 0;
				while ( pos1 < seq1.length() && pos2 < seq2.length() ) {
					if ( Character.toUpperCase(seq1.charAt(pos1++)) ==
                            Character.toUpperCase(seq2.charAt(pos2++)) ) continue;
					cnt++; // found mismatch						
				}
				return cnt;
			}

			public static int countMismatches(String seq1, String seq2, int offset2wrt1, int maxerr) {
				int pos1 = ( offset2wrt1 >= 0 ? offset2wrt1 : 0 );
				int pos2 = ( offset2wrt1 >= 0 ? 0 : -offset2wrt1 );  
				int cnt = 0;
				while ( pos1 < seq1.length() && pos2 < seq2.length() && cnt < maxerr ) {
					if ( Character.toUpperCase(seq1.charAt(pos1++)) ==
                            Character.toUpperCase(seq2.charAt(pos2++)) ) continue;
					cnt++; // found mismatch						
				}
				return cnt;
			}
			
			/** Returns a (multiline) string that represents the alignment visually: the sequences are appropriately
			 *  shifted and ready for printout; the pairwise alignment is followed by a stats line 
			 */
			public String toString() {
				StringBuffer b = new StringBuffer();
				int skip1 = ( alignment_offset >= 0 ? 0 : -alignment_offset );
				int skip2 = ( alignment_offset >=0 ? alignment_offset : 0 );
				for ( int k = 0 ; k < skip1 ; k++ ) b.append(' ');
				b.append(s1);
				b.append('\n');
				for ( int k = 0 ; k < skip2 ; k++ ) b.append(' ');
				b.append(s2);
				b.append('\n');
				b.append(best_mm+" mismatches, "+ next_mm + " next best, " + getOverlap() + " overlapping bases, distance=" + distance() + "\n");
				return b.toString();
			}

			public double distance() {
				int L = getOverlap();
				if ( L <=0 ) return 1e100;
				double l = ( best_mm==0? 1.0 : (double)best_mm + Math.sqrt((double)best_mm) );
				return ( l / (double)L );	
			}

            public static double distance(String seq1, String seq2, int offset2wrt1) {
                int L = getOverlap(seq1,seq2,offset2wrt1);
                if ( L <= 0 ) return 1e100;
                int mm = countMismatches(seq1,seq2,offset2wrt1);
                double l = ( mm == 0 ? 1.0 : (double)mm + Math.sqrt((double)mm) );
                return ( l / (double) L );
            }

}

	
