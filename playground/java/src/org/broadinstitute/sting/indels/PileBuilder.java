package org.broadinstitute.sting.indels;

import net.sf.samtools.SAMRecord;

import java.util.Collection;


public class PileBuilder implements RecordPileReceiver {
	private StrictlyUpperTriangularMatrix distances ;
	private Matrix<PairwiseAlignment> alignments ;
    private static final int KmerSize = 8;
    private MultipleAlignment alignments1;
    private MultipleAlignment alignments2;

    private static class SelectedPair {
		private int i_;
		private int j_;
		private double d_;
		
		private SelectedPair(int i, int j, double d) {
			set(i,j,d);
		}
		
		private SelectedPair() {
			set(-1,-1,1e100);
		}
		
		private double d() { return d_; }
		private int i() { return i_; }
		private int j() { return j_; }

		private void set(int i, int j, double d) {
			i_ = i;
			j_ = j;
			d_ = d;
		}

		/** Returns true if any of the two indices kept by this pair is equal to i.
		 * 
		 * @param i
		 * @return
		 */
		private boolean contains(int i) {
			return (  ( i_ == i ) || ( j_ == i ) );
		}
	
	}
	
	public PileBuilder() {}

    public void receive(Collection<SAMRecord> c) {
           IndexedSequence[] seqs = new IndexedSequence[c.size()];
           int i = 0;
           for ( SAMRecord r : c ) {
                seqs[i++] = new IndexedSequence(r.getReadString(),KmerSize);
           }
           doMultipleAlignment(seqs);
           System.out.print("Distance between final piles: "+distance(alignments1, alignments2));
           System.out.println("; diameter of PILE1: "+ diameter(alignments1));
           System.out.println("; diameter of PILE2: "+ diameter(alignments2));

		   System.out.println("PILE 1: \n"+alignments1.toString());
		   System.out.println("PILE 2: \n"+alignments2.toString());

    }

	public void initPairwiseAlignments( IndexedSequence [] seqs ) {
		 distances = new StrictlyUpperTriangularMatrix( seqs.length );
		 alignments = new Matrix<PairwiseAlignment>( seqs.length );
		 for( int i = 0; i < seqs.length ; i++ ) {
			 for ( int j = i+1 ; j < seqs.length ; j++ ) {
					PairwiseAlignment a = new PairwiseAlignment(seqs[i],seqs[j],i,j); // compute pairwise alignment
					alignments.set(i, j, a); // save it
					distances.set(i,j,a.distance());
			 }
		 }
	}	
	
	/** Finds the best pairwise alignment across all available ones. The object must be initialized first,
	 * so that the alignments are pre-computed.
	 * @return id's of the two sequences and the distance between them in a compound object.
	 */
	public SelectedPair findClosestPair() {
		
		SelectedPair p = new SelectedPair(-1,-1,1e100);
		
		for( int i = 0; i < distances.size() ; i++ ) {
			for ( int j = i+1 ; j < distances.size() ; j++ ) {
				double d = distances.get(i,j); 
				if ( d < p.d() ) p.set(i,j,d);
			}
		}
		return p;
	}

	/** Finds the worst pairwise alignment across all available ones. The object must be initialized first,
	 * so that the alignments are pre-computed.
	 * @return id's of the two sequences and the distance between them in a compound object.
	 */
	public SelectedPair findWorst() {
		
		SelectedPair p = new SelectedPair(-1,-1,-1.0);
		
		for( int i = 0; i < distances.size() ; i++ ) {
			for ( int j = i+1 ; j < distances.size() ; j++ ) {
				double d = distances.get(i,j); 
				if ( d > p.d() ) p.set(i,j,d);
			}
		}
		return p;
	}
	
	
	/** Finds the best pairwise alignment across all available ones, subject to the constraint that neither
	 *  of the two sequences found can be listed (by its id) in the supplied SelectedPair object. If the best pair is passed
	 *  as an argument, this method will find the next best pair. 
	 * 
	 * @param pexclude neither of the two sequences in the returned pair can have its id listed in pexclude pair. 
	 * @return Best pairwise alignment excluding alignments between pairs involving at least one sequence from pexclude
	 */
	public SelectedPair findNextClosestPairAfter(SelectedPair pexclude) {
		
		SelectedPair p = new SelectedPair(-1,-1,1e100);
		
		for( int i = 0; i < distances.size() ; i++ ) {
			if ( pexclude.contains(i) ) continue;
			for ( int j = i+1 ; j < distances.size() ; j++ ) {
				if ( pexclude.contains(j)) continue;
				double d = distances.get(i,j); 
				if ( d < p.d() ) p.set(i,j,d);
			}
		}
		return p;
	}
	
	/** Finds the closest sequence to the specified pile among all sequences, which are not yet in that pile. Being
	 *  the 'closest' is defined in terms of minimum distance.
	 * 
	 * @param a alignment pile to find the closest sequence for
	 * @return a compound SelectedPair object that contains the index of the closest sequence found (is guaranteed to
	 * be not in the pile), the index of the sequence in the pile it is closest to, and the actual distance between the two. 
	 */
	public SelectedPair findClosestToPile(MultipleAlignment a) {

		SelectedPair p = new SelectedPair(-1,-1,1e100);
		
		for ( Integer id : a ) {
			for (int i = 0; i < id; i++) {
					if (a.contains(i))	continue; // a already contains both sequences (i,id)
					double d = distances.get(i, id);
					if (d < p.d() ) p.set(i,id,d);
			}
			for (int j = id + 1; j < distances.size() ; j++) {
					if (a.contains(j))	continue; // a already contains both sequences (id, j)
					double d = distances.get(id, j);
					if (d < p.d()) p.set(id,j,d);
			}
		}
		return p;
	}

	
	/** Finds, among all sequences, the one farthest from the specified pile. Being
	 *  the 'farthest' is defined as having the largest lower bound of the distances to all sequences in the pile.
	 * 
	 * @param a alignment pile to find the closest sequence for
	 * @return index of the farthest sequence 
	 */
	public int findFarthestFromPile(MultipleAlignment a) {

		double dmaxmin = 0;
		int i_out = -1;
		
		for ( int i = 0 ; i < distances.size() ; i++ ) {

			if ( a.contains(i) ) continue;
			double d_min = 1e100; // smallest distance from sequence i to the pile

			for ( Integer id : a ) {
				double d = ( i < id ? distances.get(i, id) : distances.get(id,i) );
				if (d < d_min ) d_min = d;
			}
			// d_min is the smallest distance from sequence i to pile a
			if ( d_min > dmaxmin ) {
				// sequence i is so far the farthest...
				dmaxmin = d_min;
				i_out = i;
			}
		}
		return i_out;
	}
	
	public double distance(MultipleAlignment a1, MultipleAlignment a2) {
		double d = 1e100;
		for ( Integer id1 : a1 ) {
			for ( Integer id2 : a2 ) {
				if ( id1 < id2 && distances.get(id1,id2) < d ) d = distances.get(id1,id2);
				if ( id1 > id2 && distances.get(id2,id1) < d ) d = distances.get(id2,id1);
			}
		}
		return d;
	}

    /** Computes the distances from each sequence in the pile to its closest
     * neighbor (within the same pile), and returns the greatest among these distances.
     * In other words, no sequence in the pile is farther than diameter() from its closest neighbor.
     * @param a alignment pile to compute diameter for
     * @return the greatest distance from a sequence to its closest neighbor within the pile
     */
	public double diameter(MultipleAlignment a) {
		double dmaxmin = 0.0;
		for ( Integer id1 : a ) {
			double d = 1e100;
			for ( Integer id2 : a ) {
				if ( id2 <= id1 ) continue;
				if ( distances.get(id1,id2) < d ) d = distances.get(id1,id2);
			}
            // d = distance from id1 to its closest neighbor within the pile
			if ( d < 1e99 ) System.out.printf("%8.4g",d);
			if ( d < 1e99 && d > dmaxmin ) dmaxmin = d;
		}
        // dmaxmin = the largest distance from a sequence in this pile to its closest neighbor
		System.out.println();
		return dmaxmin;
	}
	
	public static void main(String argv[]) {
		
		int K=8;
//		IndexedSequence [] seqs = testSet1(K); // initialize test set data
//		IndexedSequence [] seqs = testSet2(K); // initialize test set data
		IndexedSequence [] seqs = testSet3(K); // initialize test set data
		
		PileBuilder pb = new PileBuilder();

        pb.doMultipleAlignment(seqs);
    }

    public void doMultipleAlignment(IndexedSequence[] seqs) {
		// two piles we are going to grow until all sequences are assigned to one of them.
		// we intend to keep the piles disjoint, e.g. no sequence should be placed in both
		MultipleAlignment pile1 = new MultipleAlignment();
		MultipleAlignment pile2 = new MultipleAlignment();
	
		initPairwiseAlignments(seqs);
		
		
		// all the pairwise alignments are computed and disjoint best and next-best pairs are found

		//System.out.println( distances.format("%8.4g "));

	/*	
		SelectedPair pworst = pb.findWorst();

		pile1.add(seqs[pworst.i()].getSequence(), pworst.i());
		pile2.add(seqs[pworst.j()].getSequence(), pworst.j());
*/
				
		// initialize piles with best and next-best pairs
		SelectedPair p_best = findClosestPair();
		SelectedPair p_nextbest = findNextClosestPairAfter(p_best);
		pile1.add( alignments.get(p_best.i(), p_best.j()),p_best.i(),p_best.j());
		pile2.add( alignments.get(p_nextbest.i(), p_nextbest.j()),p_nextbest.i(),p_nextbest.j());
/*		
		System.out.println("Best pair ("+p_best.i() + "," + p_best.j()+", d="+p_best.d()+"):");
		System.out.println(pile1.toString());
		System.out.println("Next best pair ("+p_nextbest.i() + "," + p_nextbest.j()+", d="+p_nextbest.d()+ "):");
		System.out.println(pile2.toString());
*/
		SelectedPair p1 = null;
		SelectedPair p2 = null;
		
		// grow piles hierarchical clustering-style
	   while ( pile1.size() + pile2.size() < seqs.length ) {
		   // find candidate sequences closest to each of the two piles
		   p1 = findClosestToPile(pile1);
		   p2 = findClosestToPile(pile2);
		   int id1_cand = pile1.selectExternal(p1.i(), p1.j()); // id of the sequence closest to the pile 1
		   int id2_cand = pile2.selectExternal(p2.i(), p2.j()); // id of the sequence closest to the pile 2
		   if ( pile2.contains(id1_cand) && pile1.contains(id2_cand)) { 
			   // pile1 and pile 2 are mutually the closest, so we need to merge them.
			   // if piles are mutually the closest, then p1 and p2 are the same pair (id1, id2), 
			   // so we just merge on one of the (redundant) instances:
			   pile1.add(pile2, alignments.get( p1.i(), p1.j()),p1.i(), p1.j());
			   pile2.clear(); // need to reset pile 2 to something else
			   int z = findFarthestFromPile(pile1); // get sequence farthest from merged pile 1
			   pile2.add(seqs[z].getSequence(), z); // and reinitialize pile 2 with that sequence
		   } else {
				   if ( p1.d() < p2.d() ) {
					   if ( pile2.contains(id1_cand) ) {
						   pile1.add(pile2, alignments.get( p1.i(), p1.j()),p1.i(), p1.j());
						   pile2.clear(); // need to reset pile 2 to something else
						   int z = findFarthestFromPile(pile1); // get sequence farthest from merged pile 1
						   pile2.add(seqs[z].getSequence(), z); // and reinitialize pile 2 with that sequence
					   } else pile1.add( alignments.get(p1.i(), p1.j()) );
				   } else {
					   if ( pile1.contains(id2_cand) ) {
						   pile2.add(pile1, alignments.get( p2.i(), p2.j()),p2.i(), p2.j());
						   pile1.clear(); // need to reset pile 2 to something else
						   int z = findFarthestFromPile(pile2); // get sequence farthest from merged pile 1
						   pile1.add(seqs[z].getSequence(), z); // and reinitialize pile 2 with that sequence
					   } else pile2.add( alignments.get(p2.i(), p2.j()) );
			   } 
		   }
	   } // end WHILE

       alignments1 = pile1;
       alignments2 = pile2;
/*
 * System.out.println("Closest distance to the pile: " + best_d
					+ "(adding: " + best_i + "," + best_j + "):");
			System.out.println(pile.toString());
		}
*/
	}
	
	public static IndexedSequence[] testSet1(int K) {
		IndexedSequence [] seqs = new IndexedSequence[9];
		seqs[0] = new IndexedSequence("CAAAAAAAGCAAAACTCTGAAGAAAGAGAGAGAGAGGGAGAGAGGGAGAGAGAAAGGGAGAGACGATGAGAGACAG",K);                                                                                
		seqs[1] = new IndexedSequence("GCAAAACTCTGAAGAAAGAGAGAGAGAGGGAGAGAGGGAGAGAGAAAGGAAGAGACGAT",K);                                                                                         
		seqs[2] = new IndexedSequence("AACTCTGAAGAAAGAGAGAGAGAGGGAGAGAGGGAGAGAGAAAGGAAGAGACGATGAGA",K);                                                                                     
		seqs[3] = new IndexedSequence("GAGAGGGAGAGAGAAAGGAAGAGACGATGAGAGACAGAGAAGGAGAGAGAAAGTACAAAAGAACGAATGAACGAAC",K);                                       
		seqs[4] = new IndexedSequence("ACGATGAGAGACAGAGAAGGAGAGAGAAAGTACAAAAGAACGAATGAACGAACAAACTAGAAATCGAGCAGGAAAA",K);                
		seqs[5] = new IndexedSequence("GAGAGACAGAGAAGGAGAGAGAAAGTACAAAAGAACGAATGAACGAACAAACTAGAAATCGAGCAGGAACCTTGGA",K);           
		seqs[6] = new IndexedSequence("TGAGACAGAGAAGGAGAGAGAAAGTACAAAAGAACGAATGAACGAACAAACTAGAAATC",K);                           
		seqs[7] = new IndexedSequence("AGACAGAGAAGGAGAGAGAAAGTACAAAAGAACGAATGAACGAACAAACTAGAAATCGAGCAGGAACCTTGGAGGA",K);        
		seqs[8] = new IndexedSequence("AGACAGAGAAGGAGAGAGAAAGTACAAAAGAACGAATGAACGAACAAACTAGAAATCGAGCAGGAACCTTGGAGGA",K);
		return seqs;
	}
	
	public static IndexedSequence[] testSet2(int K) {
		IndexedSequence [] seqs = new IndexedSequence[11];
		seqs[0] = new IndexedSequence("TGCAATGAGATGAGATCGTGCCTCTGCACTCCAGCCTGGGCGACAGAGTGAGAGACCCTGTCTCAAAAACACAAAA",K);
		seqs[1] = new IndexedSequence("AATGAGATGAGATCGTGCCTCTGCACTCCAGCCTGGGCGACAGAGTGAGAGACCCTGTCTCAAAAACACAAAAACA",K);                                                                            
		seqs[2] = new IndexedSequence("CCTCTGCACTCCAGCCTGGGCGACAGAGTGAGAGACCCTGTCTCAAAAACACAAAAACAACAACAACAAAAAAACA",K);                                                        
		seqs[3] = new IndexedSequence("CAGAGTGAGAGACCCTGTCTCAAAAACACAAAAACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCT",K);                                       
		seqs[4] = new IndexedSequence("CAGAGTGAGAGACCCTGTCTCAAAAACACAAAAACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCT",K);                                       
		seqs[5] = new IndexedSequence("GAGACCCTGTCTCAAAAACACAAAAACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCTAAACCGAG",K);                               
		seqs[6] = new IndexedSequence("CCCTGTCTCAAAAACACAAAAACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCTAAACCGAGTGTT",K);                           
		seqs[7] = new IndexedSequence("CCAAAAACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCTAAACCGAGTG",K);                             
		seqs[8] = new IndexedSequence("CAAAAACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCTAAACCGAGTGTTGTTATCTCTGGGGAGT",K);           
		seqs[9] = new IndexedSequence("AACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCTAAACCGAGTGTTGTTATCTCTGGGTAGTTTGG",K);       
		seqs[10] = new IndexedSequence("ACAACAACAACAAAAAAACACCAATCTGAGCAAATACTGCCCTAAACCGAGTGTTGTTATCTCTGGGTAGCTTGGA",K);      
		return seqs;
	}
	
	public static IndexedSequence[] testSet3(int K) {
		IndexedSequence [] seqs = new IndexedSequence[11];
		seqs[0] = new IndexedSequence("TGGAAATTTATTTCTCAGAGTACTGGAAGCTGGGAATCCAAGATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGG",K);                                                                                   
		seqs[1] = new IndexedSequence("TGGAAATTTATTTCTCAAAGTACTGGAAGCTGGGAATCCAAGATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGG",K);                                                                                   
		seqs[2] = new IndexedSequence("GGAAATTTATTTCTCAGAGTACTGGAAGCTGGGAATCCAAGATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGGG",K);                                                                                  
		seqs[3] = new IndexedSequence("GGAAATTTATTTCACAGAGTAATGGAAGCTGGGAATCCAAGATCAAAATGCCAGCAGCTTCTAAGTCTGCTGAGGG",K);                                                                                  
		seqs[4] = new IndexedSequence("ATTTCTCAGAGTACTGGAAGCTGGGAATCCAAGATCGAAATGCCAGCAGATTCTAAGTC",K);                                                                                                
		seqs[5] = new IndexedSequence("ATTTCTCAGAGTACTGGAAGCTGGGACTCCAAGATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGGGTAGGGTGC",K);                                                                     
		seqs[6] = new IndexedSequence("GTACTGGAAGCTGGGAATCCAAGATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGGGTAGGGTGCACTCTCTGCT",K);                                                                
		seqs[7] = new IndexedSequence("AATCCAAGATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGGGTAGGGTGCACTCTCTGCTTCATAAATGGGTCTC",K);                                                 
		seqs[8] = new IndexedSequence("CAAGATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGGGTAGGGCGCACTCTCTGCTTCATAAATGGGTCTCTTGC",K);                                             
		seqs[9] = new IndexedSequence("ATCAAAATGCCAGCAGATTCTAAGTCTGGTGAGGGTAGGGTGCACTCTCTGCTTCATAAATGGGTCTCTTGCCGCA",K);                                         
		seqs[10] = new IndexedSequence("GTCTGGTGAGGGTAGGGTGCACTCTCTGCTTCATAAATGGGTCTCTTGCCGCAAAAAAATCTGTTTGCTCCTCCAG",K); 		
		return seqs;
	}
}
