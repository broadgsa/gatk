/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.utils;

import org.broadinstitute.sting.utils.collections.PrimitivePair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import net.sf.samtools.util.StringUtil;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Aug 3, 2010
 * Time: 2:20:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class Assembly {
        private byte[] consensus;
        private short[] coverage;
        private short[] mismatches;
        private short [][] base_counts;

        private boolean debug = false;
        private List<String> seq_ids;
        private List<byte []> seqs;
        private List<Integer> seq_offsets;

        private KmerIndex lookup; // assembled consensus sequence is indexed here

        private int hookedAt = -1; // if set, specifies start on the ref of the assembled consensus sequence

        private static List<PrimitivePair.Int> EMPTY_KMER_LIST = new ArrayList<PrimitivePair.Int>(0);

        private int K = 15;

        private AlignmentStrategy strategy = null;
    
    /** Creates new assembly seeded with the specified sequence; default key length (15) is used.
     *
     * @param seq
     * @param id
     */
        public Assembly(final byte[] seq, String id) {
            this(15,seq,id);
        }

    /** Creates new assembly seeded with the specified sequence and sets kmer (key) length K for the internally maintained
     * lookup index tables.
     * @param K
     * @param seq
     * @param id
     */
        public Assembly(int K, final byte[] seq, String id) {
            this.K = K;
            seq_ids = new ArrayList<String>();
            seq_offsets = new ArrayList<Integer>();
            seqs = new ArrayList<byte[]>();
            seq_ids.add(id);
            seq_offsets.add(0);
            seqs.add(seq);
            consensus = Arrays.copyOf(seq,seq.length);
            coverage = new short[seq.length];
            Arrays.fill(coverage,(short)1);
            mismatches = new short[seq.length]; // filled with 0's
            base_counts = new short[4][seq.length];
            for ( int i = 0 ; i < seq.length ; i++ ) {
                int j = BaseUtils.simpleBaseToBaseIndex(seq[i]);
                if ( j != -1) base_counts[j][i] = 1;
            }
            lookup = new KmerIndex(K,seq);
            strategy = new DefaultAlignmentStrategy();
        }

    /** Creates new assembly seeded with the specified sequence; default key length (15) is used and the position on the
     * reference of the entire assembly is set (as assemblly grows, position on the ref will be updated properly).
     *
     * @param seq
     * @param id
     */
        public Assembly(final byte[] seq, String id, int posOnRef) {
            this(seq,id);
            hookedAt = posOnRef;
        }

    /** Creates new assembly seeded the specified sequence and sets kmer (key) length K for the internally maintained
     * lookup index tables. Parameter <code>posOnRef</code> specifies the (initial) position of the entire assembly on the
     * ref; as the assembly grows, the position on ref will be updated properly.
     * @param K
     * @param seq
     * @param id
     */
        public Assembly(int K, final byte[] seq, String id, int posOnRef) {
            this(K,seq,id);
            hookedAt = posOnRef;
        }

    /** Returns total number of sequences currently held by this assembly.
     * 
     * @return
     */
        public int getNumSeqs() { return seqs.size() ; }

        /** Attempts to align <code>seq</code> to this assembly's consensus. Does NOT add
         * the sequence to the consensus even if it aligns! This methods returns a list of alternative
         * best alignments found (according to the strategy used) in a newly allocated AlignmentList object.
         * @param seq sequence to align to this consensus
         * @param tryRC if true, will try aligning both seq and its reverse complement; otherwise
         * only forward alignment will be attempted (i.e. best placement of the seq, as it is provided,
         * along the assembled consensus sequence)
         * @return a newly allocated alignment list; returned list can be empty if no alignments are found
         */
        public AlignmentList align(final byte[] seq, boolean tryRC) {
            return align(seq,tryRC,null);
        }

        /** Attempts to align <code>seq</code> to this assembly's consensus. Does NOT add
         * the sequence to the consensus even if it aligns! This method uses existing list of alignments
         * (which can contain alignments to a different assembly) and updates it as necessary if a better alignment
         * (or multiple better alignments) than the one(s) already held in the list is found. Reference to the
         * <i>same</i> alignment list object is returned: this method modifies it's argument. If alignment list argument
         * is <code>null</code>, new alignment list object will be allocated and returned by this method.
         *
         * @param seq sequence to align to this consensus
         * @param tryRC if true, will try aligning both seq and its reverse complement; otherwise
         * only forward alignment will be attempted (i.e. best placement of the seq, as it is provided,
         * along the assembled consensus sequence)
         * @return a newly allocated alignment list; returned list can be empty if no alignments are found
       */
        public AlignmentList align(final byte[] seq, boolean tryRC, AlignmentList a) {
            if ( debug ) System.out.println("Assembly:: aligning sequence of length "+seq.length+"; tryRC="+tryRC+"; K="+K);

            List<PrimitivePair.Int> fw_kmers = KmerIndex.toKeyOffsetList(K,seq);

            if ( debug ) {
                for( PrimitivePair.Int kmer: fw_kmers) {
                    System.out.println("id="+kmer.getFirst()+" seq="+new String(KmerIndex.idToSeq(K,kmer.getFirst()))+" offset on seq="+kmer.getSecond());
                }
            }

            byte [] rc_seq = (tryRC ? BaseUtils.simpleReverseComplement(seq) : null );
            List<PrimitivePair.Int> rc_kmers = (tryRC ? KmerIndex.toKeyOffsetList(K,rc_seq) : EMPTY_KMER_LIST );

            if ( a == null ) a = new AlignmentList(strategy);

            // i is the position on the sequence seq or on its reverse complement
            for(PrimitivePair.Int kmer : fw_kmers ) {

                List<Integer> offsets = lookup.getOffsets(kmer.first);
                if ( offsets != null ) {
                    // kmer present in consensus sequence
                    for ( int s : offsets ) {      // s=offset of the current kmer on the assembled consensus
                        int trial_offset = s - kmer.second; // offset of the seq on the assembled consensus suggested by current kmer/offset
                        int trial_mm = countMismatches(seq,trial_offset,a.getNextBestMMCount());
                        a.tryAdd(new AlignmentInfo(trial_mm,trial_offset,false,overlap(trial_offset,seq.length),this));
                    }
                }
            }

            for ( PrimitivePair.Int kmer : rc_kmers ) {

                List<Integer> offsets = lookup.getOffsets(kmer.first);
                if ( offsets != null ) {
                    // kmer present in consensus sequence
                    for ( int s : offsets ) {
                        int trial_offset = s - kmer.second;
                        int trial_mm = countMismatches(rc_seq,trial_offset,a.getNextBestMMCount());
                        a.tryAdd(new AlignmentInfo(trial_mm,trial_offset,true,overlap(trial_offset,seq.length),this));
                    }
                }
            }
            return a;
        }

        public void setDebug(boolean d) { this.debug = d; lookup.setDebug(d);}

        public int numSequences() { return seq_ids.size(); }

        private int overlap(int offset, int seq_length ) {
            return Math.min(consensus.length,offset+seq_length)-Math.max(0,offset);
        }

        private int countMismatches(final byte seq[], int offset, int cutoff) {
            int mm = 0;

            int i ;
            if ( offset >= 0 ) i = 0;
            else { i = (-offset); offset = 0; }
            for ( ; i < seq.length && offset < consensus.length ; i++ , offset++ ) {
                if ( seq[i] != consensus[offset] ) {
                    mm++;
                    if ( mm > cutoff ) break;
                }
            }

            return mm;
        }

        public byte[] getConsensus() { return consensus; }

        public int getPosOnRef() { return hookedAt; }

        public int getConsensusLength() { return consensus.length; }

        public List<Integer> getOffsets() { return seq_offsets; }
        public int getOffset(int i) { return seq_offsets.get(i); }

        public List<String> getSeqIds() { return Collections.unmodifiableList(seq_ids); }

    /** Adds specified sequence to this assembly according to the provided
     * alignment information. Will properly update consensus sequence of this assembly
     * and all associated information (mismatches, base counts etc)
     * @param seq
     * @param id
     * @param a
     */
        public void add(final byte[] seq, String id, AlignmentInfo a) {

            if ( ! a.isAligned() ) throw new StingException("Can not add sequence to the assembly: provided alignment is empty");

            seq_ids.add(id);

            int offset = a.getOffset();
            int oldConsensusLength = consensus.length;

            byte [] seq_to_add = ( a.isNegativeStrand() ? BaseUtils.simpleReverseComplement(seq) : seq);

            seqs.add(seq_to_add);

            int pos_on_seq = 0;
            int pos_on_cons = 0;

            int leftExtension = 0; // how many bases we added to the consensus on the left
            int rightExtension = 0; // how many bases we added to the consensus on the right

            if ( offset < 0 ) {
                // if sequence sticks out to the left of the current consensus:

                leftExtension = -offset;
                for(int i = 0 ; i < seq_offsets.size() ; i++ ) {
                    // we are going to extend consensus to the left, so we need to update all current offsets:
                    seq_offsets.set(i,seq_offsets.get(i)+leftExtension);
                }

                if ( hookedAt > 0 ) hookedAt -= leftExtension;
                // extend consensus and associated arrays to the left :

                consensus = Utils.extend(consensus,offset,(byte)0); // remember, offset is negative here, extending to the left
                coverage = Utils.extend(coverage,offset,(short)1) ;
                mismatches = Utils.extend(mismatches,offset,(short)0);
                for ( int i = 0 ; i < 4 ; i++ ) base_counts[i] = Utils.extend(base_counts[i],offset,(short)0);

                for ( int j = 0 ; j < -offset ; j++ ) {
                    consensus[j] = seq_to_add[j];
                    int b = BaseUtils.simpleBaseToBaseIndex(seq_to_add[j]);
                    if ( b != -1 ) base_counts[b][j]=1;
                }

                pos_on_seq = pos_on_cons = -offset;

                offset = 0;
            }
            if ( offset > 0 ) pos_on_cons = offset;

            seq_offsets.add(offset);

            boolean consensus_changed = false;

            for ( ; pos_on_seq < seq_to_add.length && pos_on_cons < consensus.length ; pos_on_seq++, pos_on_cons++ ) {
                coverage[pos_on_cons]++;
                final byte base = seq_to_add[pos_on_seq];
                final int b = BaseUtils.simpleBaseToBaseIndex(base);
                if ( b != -1 ) {
                    // if base on seq is not a regular base, there is nothing to do;
                    // otherwise count mismatches and optionally update consensus if current base tips the balance
                    base_counts[b][pos_on_cons]++;
                    int maxcount = 0;
                    int maxb = -1;
                    for ( int j = 0 ; j < 4 ; j++ ) {
                        if ( base_counts[j][pos_on_cons] > maxcount ) {
                            maxcount = base_counts[j][pos_on_cons];
                            maxb = j;
                        }
                    }
                    // we are guaranteed here that maxb != -1 since we just added one regular base (the current one)
                    // few lines above...
                    byte newbase = BaseUtils.baseIndexToSimpleBase(maxb);
                    if ( newbase != consensus[pos_on_cons] ) { // need to change the consensus base (will recompute mismatches)
                        consensus[pos_on_cons] = newbase;
                        consensus_changed = true;
                        mismatches[pos_on_cons] = 0;
                        for ( int i = 0 ; i < 4 ; i++ ) {
                             if ( i == maxb ) continue;
                             mismatches[pos_on_cons] += base_counts[i][pos_on_cons];
                        }
                    } else { // consensus base did not change; just increment mismatches if current sequence's base differs from consensus
                        if ( base != consensus[pos_on_cons]) mismatches[pos_on_cons]++;
                    }
                }

            }

            // Last step: if sequence sticks out of current consensus on the right, we need to extend the latter:

            if ( pos_on_seq < seq_to_add.length ) {
                // sequence sticks out of consensus to the right
                rightExtension = seq_to_add.length - pos_on_seq;
                consensus = Utils.extend(consensus,rightExtension,(byte)0);
                coverage = Utils.extend(coverage,rightExtension,(short)1);
                mismatches = Utils.extend(mismatches,rightExtension,(short)0);
                for ( int i = 0 ; i < 4 ; i++ ) base_counts[i] = Utils.extend(base_counts[i],rightExtension,(short)0);
                for ( ; pos_on_seq < seq_to_add.length ; pos_on_seq++, pos_on_cons++ ) {
                    byte base = seq_to_add[pos_on_seq];
                    consensus[pos_on_cons] = base;
                    int b = BaseUtils.simpleBaseToBaseIndex(base);
                    if ( b != -1 ) base_counts[b][pos_on_cons] = base;
                }
            }

            // finally, the new sequence we just added could have mismatches that tip some consensus bases into new values;
            // let's catch those cases:


            for ( int i = 0 ; i < consensus.length ; i++ ) {
                byte cons_base = consensus[i];
                int b = BaseUtils.simpleBaseToBaseIndex(cons_base);
            }

            // there is probably a better way, but for now we just recompute the whole lookup table when consensus
            // changes somewhere in the middle (if we want to be samrt we need to identify just the kmers that changed
            // and find/change them in the lookup table).
            if ( consensus_changed ) {
                lookup.clear();
                lookup.index(consensus);
            } else {
                if ( leftExtension > 0 || rightExtension > 0 ) lookup.updateIndex(consensus,leftExtension,oldConsensusLength);
            }


        }

        public String toAlignmentString(boolean mismatchesOnly, boolean printNames) {

            int maxNameLength = 0;
            int spacing=3;

            if ( printNames ) {
                for ( String n : seq_ids ) if ( n.length() > maxNameLength ) maxNameLength++;
            }

            StringBuilder b = new StringBuilder();
            if ( printNames ) b.append(Utils.dupString(' ',maxNameLength+spacing));
            b.append(new String(consensus));
            b.append('\n');

            for ( int j = 0; j < seqs.size() ; j++ ) {
                int offset = seq_offsets.get(j);
                byte [] seq = seqs.get(j);

                if ( printNames ) {
                    b.append(seq_ids.get(j));
                    b.append(Utils.dupString(' ',maxNameLength-seq_ids.get(j).length()+spacing));
                }

                for ( int i = 0 ; i < offset ; i++ ) b.append(' ');

                for ( int i = 0 ; i < seq.length ; i++ ) {

                    byte base = seq[i];
                    if ( mismatchesOnly && base == consensus[i+offset] ) {
                        b.append('.');
                    } else b.append((char)base);
                }
                b.append('\n');
            }
            return b.toString();
        }


        public static void testMe(String [] argv ) {
            byte [] seq1 =      "ACGTTGCGTGGTTCACTGCAGTAACTGACTGATGCA".getBytes();
            byte [] seq2 =           "GCGTGGTTTACTGCAGTAACTGACTGATGCAACGTGTTTG".getBytes();
            byte [] seq3 = "GGNTGACGTTGCGTGGTTTACTGCAGTAACTGACT".getBytes();
            byte [] seq4 =      "NNNTTNCGTGGTTTACTGCAGTAACTGACTGATGCA".getBytes();

            Assembly a = new Assembly(seq1,"1");

            AlignmentList al = a.align(seq2,false);
            if ( al.isAligned() ) System.out.println("seq 2 aligned");
            else System.out.println("seq 2 did NOT align");

            if ( al.size() == 1 ) a.add(seq2,"2",al.getAlignments().get(0));
            else System.out.println("Multiple alignments found for seq 2");

            al = a.align(seq3,false);
            if ( al.isAligned() ) System.out.println("seq 3 aligned");
            else System.out.println("seq 3 did NOT align");

            if ( al.size() == 1 ) a.add(seq3,"3",al.getAlignments().get(0));
            else System.out.println("Multiple alignments found for seq 3");

            al = a.align(seq4,false);
            if ( al.isAligned() ) System.out.println("seq 4 aligned");
            else System.out.println("seq 4 did NOT align");

            if ( al.size() == 1 ) a.add(seq4,"4",al.getAlignments().get(0));
            else System.out.println("Multiple alignments found for seq 4");

            System.out.println(a.toAlignmentString(true, true));

        }

}
