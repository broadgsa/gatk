#ifndef BWA_GATEWAY
#define BWA_GATEWAY

#include <cstdio>

#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"

class Alignment;

class BWA {
 private:
  bntseq_t *bns;
  ubyte_t* reference;
  bwt_t* bwts[2];
  gap_opt_t options;

  void load_default_options();
  bwa_seq_t* create_sequence();
  void copy_bases_into_sequence(bwa_seq_t* sequence, const char* bases, const unsigned read_length);
  void create_alignments(bwa_seq_t* sequence, Alignment*& alignments, unsigned& num_alignments);

 public:
  BWA(const char* ann_filename,
      const char* amb_filename,
      const char* pac_filename,
      const char* forward_bwt_filename, 
      const char* forward_sa_filename, 
      const char* reverse_bwt_filename, 
      const char* reverse_sa_filename);
  ~BWA();

  // Parameterize the aligner.
  void set_max_edit_distance(float edit_distance);
  void set_max_gap_opens(int max_gap_opens);
  void set_max_gap_extensions(int max_gap_extensions);
  void set_disallow_indel_within_range(int indel_range);
  void set_mismatch_penalty(int penalty);
  void set_gap_open_penalty(int penalty);
  void set_gap_extension_penalty(int penalty);

  // Perform the alignment
  void align(const char* bases, const unsigned read_length, Alignment*& alignments, unsigned& num_alignments);
};

class Alignment {
 public:
  uint32_t type;
  int contig;
  bwtint_t pos;
  bool negative_strand;

  uint16_t *cigar;
  int n_cigar;

  uint32_t mapQ;
};

#endif // BWA_GATEWAY
