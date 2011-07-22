#include <cstdio>
#include <cstring>
#include <cstdlib>

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
  // Load the bns (?) and reference
  bns = bns_restore_core(ann_filename,amb_filename,pac_filename);
  reference = new ubyte_t[bns->l_pac/4+1];
  rewind(bns->fp_pac);
  fread(reference, 1, bns->l_pac/4+1, bns->fp_pac);
  fclose(bns->fp_pac);
  bns->fp_pac = NULL;

  // Load the BWTs (both directions) and suffix arrays (both directions)
  bwts[0] = bwt_restore_bwt(forward_bwt_filename);
  bwt_restore_sa(forward_sa_filename, bwts[0]);
  bwts[1] = bwt_restore_bwt(reverse_bwt_filename);
  bwt_restore_sa(reverse_sa_filename, bwts[1]);
  load_default_options();

  // Always reinitialize the random seed whenever a new set of files are loaded.
  initialize_random_seed();

  // initialize the bwase subsystem
  bwase_initialize();
}

BWA::~BWA() {
  delete[] reference;
  bns_destroy(bns);
  bwt_destroy(bwts[0]);
  bwt_destroy(bwts[1]);
}

void BWA::find_paths(const char* bases, const unsigned read_length, bwt_aln1_t*& paths, unsigned& num_paths, unsigned& best_path_count, unsigned& second_best_path_count) 
{
  bwa_seq_t* sequence = create_sequence(bases, read_length);

  // Calculate the suffix array interval for each sequence, storing the result in sequence->aln (and sequence->n_aln).
  // This method will destroy the contents of seq and rseq.
  bwa_cal_sa_reg_gap(0,bwts,1,sequence,&options);

  paths = new bwt_aln1_t[sequence->n_aln];
  memcpy(paths,sequence->aln,sequence->n_aln*sizeof(bwt_aln1_t));
  num_paths = sequence->n_aln;

  // Call aln2seq to initialize the type of match present.
  bwa_aln2seq(sequence->n_aln,sequence->aln,sequence);
  best_path_count = sequence->c1;
  second_best_path_count = sequence->c2;

  bwa_free_read_seq(1,sequence);
}

Alignment* BWA::generate_single_alignment(const char* bases, const unsigned read_length) {
  bwa_seq_t* sequence = create_sequence(bases,read_length);

  // Calculate paths.
  bwa_cal_sa_reg_gap(0,bwts,1,sequence,&options);

  // Check for no alignments found and return null.
  if(sequence->n_aln == 0) {
    bwa_free_read_seq(1,sequence);
    return NULL;
  }

  // bwa_cal_sa_reg_gap destroys the bases / read length.  Copy them back in.
  copy_bases_into_sequence(sequence,bases,read_length);

  // Pick best alignment and propagate its information into the sequence.
  bwa_aln2seq(sequence->n_aln,sequence->aln,sequence);

  // Generate the best alignment from the sequence.
  Alignment* alignment = new Alignment;
  *alignment = generate_final_alignment_from_sequence(sequence);

  bwa_free_read_seq(1,sequence);

  return alignment;
}

void BWA::generate_alignments_from_paths(const char* bases, 
                                         const unsigned read_length, 
                                         bwt_aln1_t* paths, 
                                         const unsigned num_paths, 
                                         const unsigned best_count,
                                         const unsigned second_best_count,
                                         Alignment*& alignments, 
                                         unsigned& num_alignments) 
{
  bwa_seq_t* sequence = create_sequence(bases,read_length);

  sequence->aln = paths;
  sequence->n_aln = num_paths;

  // (Ab)use bwa_aln2seq to propagate values stored in the path out into the sequence itself.
  bwa_aln2seq(sequence->n_aln,sequence->aln,sequence);

  // But overwrite key parts of the sequence in case the user passed back only a smaller subset
  // of the paths.
  sequence->c1 = best_count;
  sequence->c2 = second_best_count;
  sequence->type = sequence->c1 > 1 ? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;

  num_alignments = 0;
  for(unsigned i = 0; i < (unsigned)sequence->n_aln; i++)
    num_alignments += (sequence->aln + i)->l - (sequence->aln + i)->k + 1;

  alignments = new Alignment[num_alignments];
  unsigned alignment_idx = 0;

  for(unsigned path_idx = 0; path_idx < (unsigned)num_paths; path_idx++) {
    // Stub in a 'working' path, so that only the desired alignment is local-aligned.
    const bwt_aln1_t* path = paths + path_idx;
    bwt_aln1_t working_path = *path;

    // Loop through all alignments, aligning each one individually.
    for(unsigned sa_idx = path->k; sa_idx <= path->l; sa_idx++) {
      working_path.k = working_path.l = sa_idx;
      sequence->aln = &working_path;
      sequence->n_aln = 1;

      sequence->sa = sa_idx;
      sequence->strand = path->a;
      sequence->score = path->score;

      // Each time through bwa_refine_gapped, seq gets reversed.  Revert the reverse.
      // TODO: Fix the interface to bwa_refine_gapped so its easier to work with.
      if(alignment_idx > 0)
        seq_reverse(sequence->len, sequence->seq, 0);

      // Copy the local alignment data into the alignment object.
      *(alignments + alignment_idx) = generate_final_alignment_from_sequence(sequence);

      alignment_idx++;
    }
  }

  sequence->aln = NULL;
  sequence->n_aln = 0;

  bwa_free_read_seq(1,sequence);
}

Alignment BWA::generate_final_alignment_from_sequence(bwa_seq_t* sequence) {
  // Calculate the local coordinate and local alignment.
  bwa_cal_pac_pos_core(bwts[0],bwts[1],sequence,options.max_diff,options.fnr);
  bwa_refine_gapped(bns, 1, sequence, reference, NULL);

  // Copy the local alignment data into the alignment object.
  Alignment alignment;

  // Populate basic path info
  alignment.edit_distance = sequence->nm;
  alignment.num_mismatches = sequence->n_mm;
  alignment.num_gap_opens = sequence->n_gapo;
  alignment.num_gap_extensions = sequence->n_gape;
  alignment.num_best = sequence->c1;
  alignment.num_second_best = sequence->c2;
  
  // Final alignment position.
  alignment.type = sequence->type;
  bns_coor_pac2real(bns, sequence->pos, pos_end(sequence) - sequence->pos, &alignment.contig);
  alignment.pos = sequence->pos - bns->anns[alignment.contig].offset + 1;
  alignment.negative_strand = sequence->strand;
  alignment.mapping_quality = sequence->mapQ;
  
  // Cigar step.
  alignment.cigar = NULL;
  if(sequence->cigar) {
    alignment.cigar = new uint16_t[sequence->n_cigar];
    memcpy(alignment.cigar,sequence->cigar,sequence->n_cigar*sizeof(uint16_t));
  }
  alignment.n_cigar = sequence->n_cigar;

  // MD tag with a better breakdown of differences in the cigar
  alignment.md = strdup(sequence->md);
  delete[] sequence->md;
  sequence->md = NULL;

  return alignment;
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

void BWA::initialize_random_seed()
{
  srand48(bns->seed);  
}

void BWA::set_max_edit_distance(float edit_distance) { 
  if(edit_distance > 0 && edit_distance < 1) {
    options.fnr = edit_distance;
    options.max_diff = -1;
  }
  else {
    options.fnr = -1.0;
    options.max_diff = (int)edit_distance;
  }
}

void BWA::set_max_gap_opens(int max_gap_opens) { options.max_gapo = max_gap_opens; }
void BWA::set_max_gap_extensions(int max_gap_extensions) { options.max_gape = max_gap_extensions; }
void BWA::set_disallow_indel_within_range(int indel_range) { options.indel_end_skip = indel_range; }
void BWA::set_mismatch_penalty(int penalty) { options.s_mm = penalty; }
void BWA::set_gap_open_penalty(int penalty) { options.s_gapo = penalty; }
void BWA::set_gap_extension_penalty(int penalty) { options.s_gape = penalty; }

/**
 * Create a sequence with a set of reasonable initial defaults.  
 * Will leave seq and rseq empty.
 */
bwa_seq_t* BWA::create_sequence(const char* bases, const unsigned read_length) 
{
  bwa_seq_t* sequence = new bwa_seq_t;

  sequence->tid = -1;

  sequence->name = 0;

  copy_bases_into_sequence(sequence, bases, read_length);

  sequence->qual = 0;
  sequence->aln = 0;
  sequence->md = 0;

  sequence->cigar = NULL;
  sequence->n_cigar = 0;

  sequence->multi = NULL;
  sequence->n_multi = 0;

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
