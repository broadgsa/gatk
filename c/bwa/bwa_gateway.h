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
  bwt_t* bwts[2];
  gap_opt_t options;

  void load_default_options();
  void initialize_sequence(bwa_seq_t& sequence);
  void copy_bases_into_sequence(bwa_seq_t& sequence, const char* bases, const unsigned read_length, const bool reverse);
  void create_alignments(bwa_seq_t& sequence, Alignment*& alignments, unsigned& num_alignments);

 public:
  BWA(const char* ann_filename,
      const char* amb_filename,
      const char* pac_filename,
      const char* forward_bwt_filename, 
      const char* forward_sa_filename, 
      const char* reverse_bwt_filename, 
      const char* reverse_sa_filename);
  ~BWA();
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
