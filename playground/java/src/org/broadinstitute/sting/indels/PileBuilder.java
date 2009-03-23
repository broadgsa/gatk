package org.broadinstitute.sting.indels;

import net.sf.samtools.SAMRecord;

import java.util.*;
import org.broadinstitute.sting.utils.PrimitivePair;


public class PileBuilder implements RecordPileReceiver {
	private SymmetricMatrix distances ;
	private Matrix<PairwiseAlignment> alignments ;
    private static final int KmerSize = 8;
    private MultipleAlignment alignments1;
    private MultipleAlignment alignments2;
    private String referenceSequence;
    private int reference_start;


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

    public class SelectedSequence {
        private int id_;
        private double d_;

        private SelectedSequence(int i, double d) {
            set(i,d);
        }

        private SelectedSequence() { this(-1,1e100) ; }
        private void set(int i, double d) { id_ = i; d_ = d; }

        public double d() { return d_;}
        public int i() { return id_; }

    }
	
	public PileBuilder() {
        referenceSequence = null;
        reference_start = -1;
    }

    public void setReferenceSequence(String seq, int start) {
        referenceSequence = seq;
        reference_start = start;
    }

    public void setReferenceSequence(String seq) {
        referenceSequence = seq;
        reference_start = -1;
    }

    public void receive(Collection<SAMRecord> c) {
           IndexedSequence[] seqs = new IndexedSequence[c.size()];
           int i = 0;
           int left = 1000000000;
           int right = 0;
           for ( SAMRecord r : c ) {
                seqs[i++] = new IndexedSequence(r.getReadString(),KmerSize);
                left = Math.min(left, r.getAlignmentStart() );
                right = Math.max(right,r.getAlignmentEnd());
           }

           SequencePile originalAligns = new SequencePile(referenceSequence.substring(left-1,right));
           for ( SAMRecord r : c ) {
               originalAligns.addAlignedSequence(r.getReadString(), r.getReadNegativeStrandFlag(),
                       r.getCigar(), r.getAlignmentStart() - left );
           }
           System.out.println("\n#############################################################################");
           System.out.println("ORIGINAL ALIGNMENT: \n");
           originalAligns.colorprint(true);
           System.out.println("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++") ;

           List<MultipleAlignment> piles = doMultipleAlignment2(seqs);
//           System.out.print("Distance between final piles: "+distance(alignments1, alignments2));
//           System.out.print("; diameter of PILE1: "+ diameter(alignments1));
//           System.out.println("; diameter of PILE2: "+ diameter(alignments2));

//		   System.out.println("PILE 1: \n"+alignments1.toString());
//		   System.out.println("PILE 2: \n"+alignments2.toString());

        SymmetricMatrix d = new SymmetricMatrix(piles.size());
        for ( int n = 0 ; n < piles.size() ; n++ ) {
            d.set(n,n,diameter(piles.get(n)));
            for ( int m = n+1 ; m < piles.size() ; m++ ) {
                d.set(n,m,distance(piles.get(n), piles.get(m)));
            }
        }
        System.out.println(d.format("%8.4g"));
        for ( int n = 0 ; n < piles.size() ; n++ ) {
            System.out.println("PILE " + n +"\n" +piles.get(n).toString());
        }
    }

	public void initPairwiseAlignments( IndexedSequence [] seqs ) {
		 distances = new SymmetricMatrix( seqs.length );
		 alignments = new Matrix<PairwiseAlignment>( seqs.length );
		 for( int i = 0; i < seqs.length ; i++ ) {
			 for ( int j = i+1 ; j < seqs.length ; j++ ) {
					PairwiseAlignment a = new PairwiseAlignment(seqs[i],seqs[j],i,j); // compute pairwise alignment
					alignments.set(i, j, a); // save it
                    alignments.set(j, i, a); // save it
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
			for (int i = 0; i < distances.size(); i++) {
					if (a.contains(i))	continue; // a already contains both sequences (i,id)
					double d = distances.get(i, id);
					if (d < p.d() ) p.set(i,id,d);
			}
		}
		return p;
	}

    public SelectedPair findClosestToPileAverage(MultipleAlignment a) {

        SelectedPair p = new SelectedPair(-1,-1,1e100);

        // currently, we compute the average distance from each sequence to the pile, but if the average
        // distance is small enough, we will try to stitch that sequence to the pile based on the *best*
        // available pairwise alignment, best_id will keep the id of that sequence from the pile that
        // has the best alignment with the sequence that is the closest on average
        int best_id=-1;

        Set<Integer> offsets = new HashSet<Integer>(); // all putative offsets suggested by different p-wise alignments
        for ( int i = 0 ; i < distances.size() ; i++ ) { // for each sequence i

            if ( a.contains(i) ) continue; // sequence i is already in the pile, ignore it

            offsets.clear();

            for ( Integer id : a ) { // for all sequences from the msa pile
                PairwiseAlignment pa = alignments.get(i,id);
                if ( pa.getOverlap() <= 0 ) continue; // at this step we do not take into account sequences with no overlap
                // alignment pa suggests this offset of i wrt the first sequence in the msa
                offsets.add( pa.getBestOffset2wrt1(id,i)+a.getOffsetById(id));
            }
            // we got all suggested offsets; now lets recompute distances:

            for( Integer off : offsets ) {
                SelectedPair spo = averageDistanceForOffset(a,i,off);
                if ( spo.d() < p.d() ) p.set(spo.i(),spo.j(),spo.d());
            }
        }
        return p;
    }

    public Matrix<SelectedPair> averageClosestDistanceMatrix(List<MultipleAlignment> la, int n) {

        Matrix<SelectedPair> mp = new Matrix<SelectedPair>(n);

        for ( int i = 0 ; i < n ; i++ ) {
            for ( int j = i + 1 ; j < n ; j++  ) {
                mp.set(i,j, findBestAlignment(la.get(i),la.get(j)) );
                mp.set(j,i, mp.get(i,j) );
            }
        }
        return mp;
    }

    public SelectedPair findBestAlignment(MultipleAlignment a1, MultipleAlignment a2) {

        Map<Integer, PrimitivePair.Int > all_offsets = new HashMap<Integer, PrimitivePair.Int >();
        SelectedPair p = new SelectedPair(-1,-1,1e100);

        for ( Integer id1 : a1 ) {
            for ( Integer id2 : a2 ) {
                PairwiseAlignment pa = alignments.get(id1,id2);
                if ( pa.getOverlap() <= 0 ) continue; // id1 and id2 do not overlap and/or we don't have p-wise alignment

                // record suggested offset of a2 wrt a1 (defined by their first sequences), and remember the
                // pairwise alignment that suggested it
                int suggested_offset = a1.getOffsetById(id1) + pa.getBestOffset2wrt1(id1,id2) - a2.getOffsetById(id2);
                if ( ! all_offsets.containsKey(suggested_offset) ) {
                    all_offsets.put( suggested_offset , new PrimitivePair.Int(id1,id2)) ;
                }
            }
        }
        for ( Map.Entry<Integer,PrimitivePair.Int> offset_record : all_offsets.entrySet() ) {

            double d = averageDistanceForOffset(a1,a2,offset_record.getKey());
            if ( d < p.d() ) p.set(offset_record.getValue().first,offset_record.getValue().second,d);
        }
        return p;
    }

    public double averageDistanceForOffset(MultipleAlignment a1, MultipleAlignment a2, int offset) {
        SelectedPair p = new SelectedPair();

        double d_av = 0;
        int nseq = 0;
        int i1 = -1;
        int i2 = -1;

        for ( Integer id2 : a2 ) {
            SelectedPair spo = averageDistanceForOffset(a1,id2,offset+a2.getOffsetById(id2));
            if ( spo.d() > 1e99 ) continue;
            nseq++;
            d_av += spo.d();
        }
        if ( nseq == 0 ) return 1e100;
        d_av /= nseq;
        return d_av;
    }

    /** Computes average distance from sequence i to multiple alignment a for the specified offset of i wrt a
     *  and returns that distance and pair of sequence indices, on which the specified offset is realized
     * @param a
     * @param i
     * @param offset
     * @return
     */
    public SelectedPair averageDistanceForOffset(MultipleAlignment a, int i, int offset) {

        SelectedPair p = new SelectedPair(-1,-1,1e100);

        double d = 0; // will hold average distance
        double dmin = 1e100; // used to find the nearest individual sequence in the pile
        int nseq = 0; // number of sequences in the pile that have distance to sequence i defined
        int best_id = -1;

        for ( Integer id : a ) { // for all sequences from the msa pile
            PairwiseAlignment pa = alignments.get(i,id);
            int new_off = offset - a.getOffsetById(id); // offset of i wrt id as suggested by <offset>
            double dist_for_off; // distance between i and id for the given offset off

            // check if p-wise alignment has data for the specified offset:
            boolean canuse = false;
            if ( pa.alignmentExists() && pa.getBestOffset2wrt1(id,i) == new_off ) {
                dist_for_off = distances.get(i,id);
                canuse = true; // can use this alignment to stitch i to a
            }
            else {
                // offset is different from what the pwise alignment suggests; recompute!
                dist_for_off = PairwiseAlignment.distance(pa.getSequenceById(id),pa.getSequenceById(i),new_off);
            }
            if ( dist_for_off > 1e99 ) continue; // at this offset, i and id do not overlap, go check next id
            d += dist_for_off;
            nseq++;
            if ( dist_for_off < dmin && canuse ) {
                dmin = dist_for_off;
                best_id = id;
            }

        }
        if ( nseq == 0 ) return p;
        d /= (double)nseq;
        p.set(i,best_id,d);
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
				double d = distances.get(i, id) ;
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
				if ( distances.get(id1,id2) < d ) d = distances.get(id1,id2);
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
        System.out.print("\n[");
        Iterator<Integer> ids1 = a.sequenceIdByOffsetIterator();
		while ( ids1.hasNext() ) {
            Integer id1 = ids1.next();
			double d = 1e100; // will hold distance from id1 to its closest neighbor
            for ( Integer id2 : a ) {
				if ( id2 == id1 ) continue;
                double dpair = distances.get(id1,id2) ;
				d = Math.min(d,dpair);
			}
            // d = distance from id1 to its closest neighbor within the pile
			if ( d < 1e99 ) System.out.printf("%8.4g",d);
			if ( d < 1e99 && d > dmaxmin ) dmaxmin = d;
		}
        System.out.println(" ]");
        // dmaxmin = the largest distance from a sequence in this pile to its closest neighbor
//		System.out.println();
		return dmaxmin;
	}
	
	public static void main(String argv[]) {
		
		int K=8;
//		IndexedSequence [] seqs = testSet1(K); // initialize test set data
//		IndexedSequence [] seqs = testSet2(K); // initialize test set data
//        IndexedSequence [] seqs = testSet3(K); // initialize test set data
        IndexedSequence [] seqs = testSet4(K); // initialize test set data

		PileBuilder pb = new PileBuilder();

        //pb.doMultipleAlignment(seqs);
        pb.doMultipleAlignment2(seqs);
        System.out.print("Distance between final piles: "+pb.distance(pb.alignments1, pb.alignments2));
        System.out.print("; diameter of PILE1: "+ pb.diameter(pb.alignments1));
        System.out.println("; diameter of PILE2: "+ pb.diameter(pb.alignments2));

        System.out.println("PILE 1: \n"+pb.alignments1.toString());
        System.out.println("PILE 2: \n"+pb.alignments2.toString());
    }

    public void doMultipleAlignment(IndexedSequence[] seqs) {
		// two piles we are going to grow until all sequences are assigned to one of them.
		// we intend to keep the piles disjoint, e.g. no sequence should be placed in both


		MultipleAlignment pile1 = new MultipleAlignment();
		MultipleAlignment pile2 = new MultipleAlignment();
	
		initPairwiseAlignments(seqs);
		
		
		// all the pairwise alignments are computed and disjoint best and next-best pairs are found

//		System.out.println( distances.format("%8.4g "));


		SelectedPair pworst = findWorst();

		pile1.add(seqs[pworst.i()].getSequence(), pworst.i());
		pile2.add(seqs[pworst.j()].getSequence(), pworst.j());


		// initialize piles with best and next-best pairs
/*
		SelectedPair p_best = findClosestPair();
		SelectedPair p_nextbest = findNextClosestPairAfter(p_best);
		pile1.add( alignments.get(p_best.i(), p_best.j()));
		pile2.add( alignments.get(p_nextbest.i(), p_nextbest.j()));
*/
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

//		   p1 = findClosestToPileAverage(pile1); // findClosestToPile(pile1);
//		   p2 = findClosestToPileAverage(pile2); //findClosestToPile(pile2);
           p1 = findClosestToPile(pile1); // findClosestToPile(pile1);
           p2 = findClosestToPile(pile2); //findClosestToPile(pile2);
		   int id1_cand = pile1.selectExternal(p1.i(), p1.j()); // id of the sequence closest to the pile 1
		   int id2_cand = pile2.selectExternal(p2.i(), p2.j()); // id of the sequence closest to the pile 2
		   if ( pile2.contains(id1_cand) && pile1.contains(id2_cand)) { 
			   // pile1 and pile 2 are mutually the closest, so we need to merge them.
			   // if piles are mutually the closest, then p1 and p2 are the same pair (id1, id2), 
			   // so we just merge on one of the (redundant) instances:
			   pile1.add(pile2, alignments.get( p1.i(), p1.j()));
			   pile2.clear(); // need to reset pile 2 to something else
			   int z = findFarthestFromPile(pile1); // get sequence farthest from merged pile 1
			   pile2.add(seqs[z].getSequence(), z); // and reinitialize pile 2 with that sequence
		   } else {
				   if ( p1.d() < p2.d() ) {
					   if ( pile2.contains(id1_cand) ) {
						   pile1.add(pile2, alignments.get( p1.i(), p1.j()));
						   pile2.clear(); // need to reset pile 2 to something else
						   int z = findFarthestFromPile(pile1); // get sequence farthest from merged pile 1
						   pile2.add(seqs[z].getSequence(), z); // and reinitialize pile 2 with that sequence
					   } else pile1.add( alignments.get(p1.i(), p1.j()) );
				   } else {
					   if ( pile1.contains(id2_cand) ) {
						   pile2.add(pile1, alignments.get( p2.i(), p2.j()));
						   pile1.clear(); // need to reset pile 2 to something else
						   int z = findFarthestFromPile(pile2); // get sequence farthest from merged pile 1
						   pile1.add(seqs[z].getSequence(), z); // and reinitialize pile 2 with that sequence
					   } else pile2.add( alignments.get(p2.i(), p2.j()) );
			   } 
		   }
           System.out.println("PILE 1: \n"+pile1.toString());
           System.out.println("PILE 2: \n"+pile2.toString());
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

public List<MultipleAlignment> doMultipleAlignment2(IndexedSequence[] seqs) {

    initPairwiseAlignments(seqs);

    List<MultipleAlignment> piles = new LinkedList<MultipleAlignment>();

    int npiles = seqs.length;

    for ( int i = 0 ; i < seqs.length ; i++ ) {
        MultipleAlignment m = new MultipleAlignment();
        m.add(seqs[i].getSequence(),i);
        piles.add(m);
    }

    while ( npiles > 2 ) {
        Matrix<SelectedPair> dist = averageClosestDistanceMatrix(piles,npiles);
        int best_i = -1;
        int best_j = -1;
        int pile_i = -1;
        int pile_j = -1;
        double d = 1e100;
        for ( int i = 0 ; i < npiles ; i++ ) {
            for ( int j = i+1 ; j < npiles ; j++ ) {
                SelectedPair p = dist.get(i,j);
                if ( p.d() < d ) {
                    d = p.d();
                    pile_i = i;
                    pile_j = j;
                    best_i = p.i();
                    best_j = p.j();
                }
            }
        }
//        System.out.println("joining pile "+pile_i +" and pile " + pile_j +" on seqs " + best_i +" and " + best_j );
        // got the closest pair
        piles.get(pile_i).add(piles.get(pile_j),alignments.get(best_i,best_j));
//        System.out.println("JOINED PILE: \n"+piles.get(pile_i).toString());
        piles.remove(pile_j);
        npiles--;
    }

    alignments1 = piles.get(0);
    alignments2 = piles.get(1);


//    System.out.println("PILE 1: \n"+piles.get(0).toString());
//    System.out.println("PILE 2: \n"+piles.get(1).toString());
    return piles;
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

    public static IndexedSequence[] testSet4(int K) {
        IndexedSequence [] seqs = new IndexedSequence[19];
        seqs[0]  = new IndexedSequence("CGTGTGTGTGTGTGCAGTGCGTGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTTTGTGAGATC",K);
        seqs[1]  = new IndexedSequence("ATGTGTGTGTGTGCAGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCAT",K);
        seqs[2]  = new IndexedSequence("GTGTGTGTGTGCAGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGC",K);
        seqs[3]  = new IndexedSequence("TGTGTGTGTGCAGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCA",K);
        seqs[4]  = new IndexedSequence("GTGTGTGTGCAGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCAT",K);
        seqs[5]  = new IndexedSequence("GTGTGTGTGCCGTGCTTTGTGCTGTGAGATCTGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCTGCAT",K);
        seqs[6]  = new IndexedSequence("GTGTGTGCAGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGT",K);
        seqs[7]  = new IndexedSequence("GTGCAGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGT",K);
        seqs[8]  = new IndexedSequence("TGCAGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGTG",K);
        seqs[9]  = new IndexedSequence("AGTGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGTGTGT",K);
        seqs[10] = new IndexedSequence("TGGGCATGGTGCTGTGAGATCAGCGTGTGTGTGTGCAGCGCATGGTGCTGTGTGAGATCAGCGTGTGTGTGTGCAG",K);
        seqs[11]  = new IndexedSequence("GCTGTGAGATCAGCGTGTGTGTGTGAGCAGTGCATGGGGATGTGTGAGATCAGCATGTGTGTGTGTGTGCAGCGCG",K);
        seqs[12]  = new IndexedSequence("GCTGTGAGATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGTGTGTGTGCAGTGCA",K);
        seqs[13]  = new IndexedSequence("AGATCAGCATGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGTGTGTGTGCAGTGCATGGTGC",K);
        seqs[14] = new IndexedSequence("AGATCAGCGTGTGTGTGTGCAGCGCATGGCGCTGTGTGAGATCAGCATGTGTGTGTGTGTGCGGCGCATGGGGGTG",K);
        seqs[15]  = new IndexedSequence("GATCAGCGTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGAATGTGTGTGTGTGTGCAGTGCATGGTGCT",K);
        seqs[16]  = new IndexedSequence("ATCAGCATGGGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGGGTGTGTGGGGTGGGTGGTGGTG",K);
        seqs[17]  = new IndexedSequence("ATCAGCATGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGTGTGTGTGCAGTGCATGGGGCTG",K);
        seqs[18]  = new IndexedSequence("GTGTGTGTGTGTGCAGTGCATGGTGCTGTGTGAGATCAGCATGTGTGTGTGTGTGCAGTGCATGGTGCTGAGTGTG",K);
        return seqs;
    }
}
