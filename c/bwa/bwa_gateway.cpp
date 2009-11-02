#include <cstdio>
#include <cstring>

#include "bwase.h"
#include "bwa_gateway.h"

BWA::BWA(const char* ann_filename,
	 const char* amb_filename,
	 const char* pac_filename,
	 const char* forward_bwt_filename, 
	 const char* forward_sa_filename, 
	 const char* reverse_bwt_filename, 
	 const char* reverse_sa_filename) 
{
  bns = bns_restore_core(ann_filename,amb_filename,pac_filename);
  bwts[0] = bwt_restore_bwt(forward_bwt_filename);
  bwt_restore_sa(forward_sa_filename, bwts[0]);
  bwts[1] = bwt_restore_bwt(reverse_bwt_filename);
  bwt_restore_sa(reverse_sa_filename, bwts[1]);
  load_default_options();
}

BWA::~BWA() {
  bwt_destroy(bwts[0]);
  bwt_destroy(bwts[1]);
}

void BWA::align(const char* bases, const unsigned read_length, Alignment*& alignments, unsigned& num_alignments) 
{
  bwa_seq_t* sequence = create_sequence();
  copy_bases_into_sequence(sequence, bases, read_length);

  // Calculate the suffix array interval for each sequence, storing the result in sequence.aln (and sequence.n_aln).
  // This method will destroy the contents of seq and rseq.
  bwa_cal_sa_reg_gap(0,bwts,1,sequence,&options);

  // Translate suffix array indices into exactly how many alignments have been found.
  bwa_aln2seq(sequence->n_aln,sequence->aln,sequence);

  // Calculate and refine the position for each alignment.  This position may be inaccurate 
  // if the read contains indels, etc.  Refinement requires the original sequences in the proper order.
  copy_bases_into_sequence(sequence, bases, read_length);
  create_alignments(sequence, alignments, num_alignments);

  bwa_free_read_seq(1,sequence);
}

void BWA::load_default_options() 
{
  options.s_mm = 3; 
  options.s_gapo = 11; 
  options.s_gape = 4; 
  options.mode = 3; 
  options.indel_end_skip = 5; 
  options.max_del_occ = 10; 
  options.max_entries = 2000000; 
  options.fnr = 0.04; 
  options.max_diff = -1; 
  options.max_gapo = 1; 
  options.max_gape = 6; 
  options.max_seed_diff = 2; 
  options.seed_len = 2147483647; 
  options.n_threads = 1; 
  options.max_top2 = 30; 
  options.trim_qual = 0;
}

/**
 * Create a sequence with a set of reasonable initial defaults.  
 * Will leave seq and rseq empty.
 */
bwa_seq_t* BWA::create_sequence() 
{
  bwa_seq_t* sequence = new bwa_seq_t;

  sequence->tid = -1;

  sequence->name = 0;

  sequence->seq = NULL;
  sequence->rseq = NULL;
  sequence->qual = 0;
  sequence->aln = 0;
  sequence->md = 0;

  sequence->cigar = NULL;
  sequence->n_cigar = 0;

  return sequence;
}

void BWA::copy_bases_into_sequence(bwa_seq_t* sequence, const char* bases, const unsigned read_length) 
{
  // seq, rseq will ultimately be freed by bwa_cal_sa_reg_gap
  sequence->seq = new ubyte_t[read_length];
  sequence->rseq = new ubyte_t[read_length];
  for(unsigned i = 0; i < read_length; i++) sequence->seq[i] = nst_nt4_table[(unsigned)bases[i]];
  memcpy(sequence->rseq,sequence->seq,read_length);

  // BWA expects the read bases to arrive reversed.
  seq_reverse(read_length,sequence->seq,0);
  seq_reverse(read_length,sequence->rseq,1);

  sequence->full_len = sequence->len = read_length;
}

void BWA::create_alignments(bwa_seq_t* sequence, Alignment*& alignments, unsigned& num_alignments) {
  bool debug = false;

  num_alignments = 0;
  for(unsigned i = 0; i < (unsigned)sequence->n_aln; i++)
    num_alignments += (sequence->aln + i)->l - (sequence->aln + i)->k + 1;

  alignments = new Alignment[num_alignments];
  unsigned alignment_idx = 0;

  // backup existing alignment blocks.
  bwt_aln1_t* alignment_blocks = sequence->aln;
  int num_alignment_blocks = sequence->n_aln;

  for(unsigned alignment_block_idx = 0; alignment_block_idx < (unsigned)num_alignment_blocks; alignment_block_idx++) {
    // Stub in a 'working' alignment block, so that only the desired alignment is local-aligned.
    const bwt_aln1_t* alignment_block = alignment_blocks + alignment_block_idx;
    bwt_aln1_t working_alignment_block = *alignment_block;

    // Loop through all alignments, aligning each one individually.
    for(unsigned sa_idx = alignment_block->k; sa_idx <= alignment_block->l; sa_idx++) {
      working_alignment_block.k = working_alignment_block.l = sa_idx;
      sequence->aln = &working_alignment_block;
      sequence->n_aln = 1;

      sequence->sa = sa_idx;
      sequence->strand = alignment_block->a;
      sequence->score = alignment_block->score;

      // Each time through bwa_refine_gapped, seq gets reversed.  Revert the reverse.
      // TODO: Fix the interface to bwa_refine_gapped so its easier to work with.
      if(alignment_idx > 0)
	seq_reverse(sequence->len, sequence->seq, 0);

      // Calculate the local alignment.
      bwa_cal_pac_pos_core(bwts[0],bwts[1],sequence,options.max_diff,options.fnr);
      if(debug) {
	printf("alignment_idx = %d, k = %d, l = %d, sa_idx = %d\n", alignment_idx, alignment_block->k, alignment_block->l, sa_idx);
	printf("sequence->pos = %d\n",sequence->pos);
      }

      bwa_refine_gapped(bns, 1, sequence, 0, NULL);

      // Copy the local alignment data into the alignment object.
      Alignment& alignment = *(alignments + alignment_idx);
      alignment.type = sequence->type;
      bns_coor_pac2real(bns, sequence->pos, pos_end(sequence) - sequence->pos, &alignment.contig);
      alignment.pos = sequence->pos - bns->anns[alignment.contig].offset + 1;
      alignment.negative_strand = sequence->strand;
      alignment.mapQ = sequence->mapQ;

      alignment.cigar = NULL;
      if(sequence->cigar) {
	alignment.cigar = new uint16_t[sequence->n_cigar];
	memcpy(alignment.cigar,sequence->cigar,sequence->n_cigar*sizeof(uint16_t));
      }
      alignment.n_cigar = sequence->n_cigar;

      alignment_idx++;
    }
  }

  // Restore original alignment blocks.
  sequence->aln = alignment_blocks;
  sequence->n_aln = num_alignment_blocks;
}
