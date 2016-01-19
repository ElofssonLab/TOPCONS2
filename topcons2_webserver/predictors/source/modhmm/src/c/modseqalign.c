/* align two sequences via the most likely state path of a given hmm */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "structs.h"
#include "funcs.h"
#include "cmdline_modseqalign.h"


#define NORM_LOG_LIKE 0
#define LOGODDS 1

#define HMM 20
#define SEQ 21

#define LONGEST_SEQ -1 /* Note: this constant is also defined in read_seqs.c */
#define FIRST_SEQ 1

#define LEAD_SEQ 10
#define ALL_SEQS 11

#define MAX_LINE 500

//#define DEBUG_ALIGN

extern int verbose;


static struct hmm_multi_s hmm;
static struct msa_sequences_multi_s *msa_seq_infop_template, *msa_seq_infop_target;
static struct sequences_multi_s seq_info_template, seq_info_target;
static struct replacement_letter_multi_s replacement_letters;
long double *subst_mtxp;
long double *subst_mtxp_2;
long double *subst_mtxp_3;
long double *subst_mtxp_4;
long double *aa_freqs;
long double *aa_freqs_2;
long double *aa_freqs_3;
long double *aa_freqs_4;


void align_paths(int **path_templatep, int **path_targetp, int seq_len_target, int seq_len_template, FILE *outfile,
		 int seq_format);


int main(int argc, char* argv[])
{
  /* command line options */
  FILE *hmmfile; /* file to read hmms from */
  FILE *outfile; /* output file */
  FILE *target_seqfile; /* file to read sequences from */
  FILE *template_seqfile; /* file for reading sequence names */
  FILE *replfile; /* file for reading special replacement letters */
  FILE *priorfile; /* file to read priordistribution from */
  FILE *substmtxfile; /* file to read substitution matrix from */
  FILE *freqfile;
  int seq_format; /* input sequence format */
  int align_alg, label_alg; /* algorithm for calculating match score */
  int query; /* are hmms or seqs the queryside */
  char seq_name[200]; 
  char hmm_name[200];
  char cur_savefile_name[200];
  char *pure_seq_name_p, pure_seq_name_a[200]; /* sequences name without path name */
  char *pure_hmm_name_p, pure_hmm_name_a[200]; /* sequences name without path name */
  char *last_slash, *last_dot;
  int prf_mode; /* method for doing the msa-hmm search */
  int lead_seq; /* index of the sequence on which the search is performed */
  int use_gap_shares; /* account for/don't account for gaps use when calculating the share of a symbol */
  int use_prior; /* are we working with priors or not */
  int use_labels;
  int normalize;
  int scoring_method;
  int read_subst_mtx;
  int multi_scoring_method;
  int run_viterbi;
  int hmmfiletype;
  int use_lead_columns;

  struct viterbi_s *viterbi_mtx;
  long double viterbi_score;
  int *path_template, *path_incrementable_template;
  int *path_target, *path_incrementable_target;
  int seq_len_template, seq_len_target;
  int b;
  int k;
  
  struct gengetopt_args_info args_info;

  /* temporary variables */
  int i,j; /*standard loop indices */
 
  /*init some variables */
  outfile = NULL;
  hmmfile = NULL;
  target_seqfile = NULL;
  template_seqfile = NULL;
  replfile = NULL;
  priorfile = NULL;
  substmtxfile = NULL;
  freqfile = NULL;
  seq_format = STANDARD;
  align_alg = FORW_ALG;
  label_alg = ONE_BEST_ALG;
  query = SEQ;
  lead_seq = -1;
  prf_mode = ALL_SEQS;
  use_gap_shares = NO;
  use_prior = NO;
  use_labels = NO;
  normalize = NO;
  scoring_method = SJOLANDER;
  read_subst_mtx = NO;
  multi_scoring_method = JOINT_PROB;
  subst_mtxp = NULL;
  subst_mtxp_2 = NULL;
  subst_mtxp_3 = NULL;
  subst_mtxp_4 = NULL;
  aa_freqs = NULL;
  aa_freqs_2 = NULL;
  aa_freqs_3 = NULL;
  aa_freqs_4 = NULL;
  run_viterbi = NO;
  use_lead_columns = NO;


  /* parse command line */

  if(cmdline_parser(argc, argv, &args_info) != 0) {
    exit(1);
  }

  /* compulsory options */
  if(args_info.hmmfile_given) {
    if((hmmfile = fopen(args_info.hmmfile_arg, "r")) == NULL) {
      perror(args_info.hmmfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading hmm\n",args_info.hmmfile_arg);
    }
  }
  if(args_info.target_given) {
    if((target_seqfile = fopen(args_info.target_arg, "r")) == NULL) {
      perror(args_info.target_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading target sequence\n",args_info.target_arg);
    }
  }
  if(args_info.template_given) {
    if((template_seqfile = fopen(args_info.template_arg, "r")) == NULL) {
      perror(args_info.template_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading template sequence\n",args_info.template_arg);
    }
  }
  if(args_info.outfile_given) {
    if((outfile = fopen(args_info.outfile_arg, "w")) == NULL) {
      perror(args_info.outfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for writing\n",args_info.outfile_arg);
    }
  }
  
  if(args_info.seqformat_given) {
    if(strcmp(args_info.seqformat_arg, "fa") == 0) {
      seq_format = FASTA; 
    }
    else if(strcmp(args_info.seqformat_arg, "s") == 0) {
      seq_format = STANDARD;
    }
    else if(strcmp(args_info.seqformat_arg, "msa") == 0) {
      seq_format = MSA_STANDARD;
    }
    else if(strcmp(args_info.seqformat_arg, "prf") == 0) {
      seq_format = PROFILE;
    }
    else {
      printf("Incorrect sequence format: %s\n", args_info.seqformat_arg);
      exit(0);
    }
  }

  /* non compulsory options */
  if(args_info.smxfile_given) {
    if((substmtxfile = fopen(args_info.smxfile_arg, "r")) == NULL) {
      perror(args_info.smxfile_arg);
      exit(0);
    }
    else {
      read_subst_mtx = YES;
      printf("Opened file %s for reading substitution matrix\n",args_info.smxfile_arg);
    }
  }
  if(args_info.freqfile_given) {
    if((freqfile = fopen(args_info.freqfile_arg, "r")) == NULL) {
      perror(args_info.freqfile_arg);
      exit(0);
    }
    else {
      
      printf("Opened file %s for reading background frequencies\n",args_info.freqfile_arg);
    }
  }
  if(args_info.replfile_given) {
    if((replfile = fopen(args_info.replfile_arg, "r")) == NULL) {
      perror(args_info.replfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading replacement letters\n",args_info.replfile_arg);
    }
  }
  if(args_info.priorfile_given) {
    if((priorfile = fopen(args_info.priorfile_arg, "r")) == NULL) {
      perror(args_info.priorfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading priors\n",args_info.priorfile_arg);
    }
  }
 
 
  
  /* msa scoring options */
  if(args_info.msascoring_given) {
    if(strcmp(args_info.msascoring_arg, "DP") == 0) {
      scoring_method = DOT_PRODUCT;
    }
    else if(strcmp(args_info.msascoring_arg, "DPPI") == 0) {
      scoring_method = DOT_PRODUCT_PICASSO;
    }
    else if(strcmp(args_info.msascoring_arg, "PI") == 0) {
      scoring_method = PICASSO;
    }
    else if(strcmp(args_info.msascoring_arg, "PIS") == 0) {
      scoring_method = PICASSO_SYM;
    }
    else if(strcmp(args_info.msascoring_arg, "GM") == 0) {
      scoring_method = SJOLANDER;
    }
    else if(strcmp(args_info.msascoring_arg, "GMR") == 0) {
      scoring_method = SJOLANDER_REVERSED;
    }
    //else if(strcmp(args_info.msascoring_arg, "SMP") == 0) {
    //  scoring_method = SUBST_MTX_PRODUCT;
    //}
    //else if(strcmp(args_info.msascoring_arg, "SMDP") == 0) {
    //  scoring_method = SUBST_MTX_DOT_PRODUCT;
    //}
    //else if(strcmp(args_info.msascoring_arg, "SMDPP") == 0) {
    //  scoring_method = SUBST_MTX_DOT_PRODUCT_PRIOR;
    //}
    else {
      printf("Incorrect scoring method option: %s\n", args_info.msascoring_arg);
      exit(0);
    }
  }
  if(args_info.usecolumns_given) {
    if(strcmp(args_info.usecolumns_arg, "all") == 0) {
      prf_mode = ALL_SEQS;
      use_lead_columns = NO;
    }
    else {
      lead_seq = atoi(args_info.usecolumns_arg);
      prf_mode = LEAD_SEQ;
      use_lead_columns = YES;
      if(lead_seq <= 0) {
	printf("Incorrect use-column option: %s\n", args_info.usecolumns_arg);
	exit(0);
      }
    }
  }

  /* flags */
  if(args_info.nolabels_given) {
    /* checked after seqread */
  }
  if(args_info.verbose_given) {
    verbose = YES;
  }
 

  /* read subst mtx */
  if(substmtxfile != NULL) {  
    read_subst_matrix_multi(&subst_mtxp, &subst_mtxp_2, &subst_mtxp_3, &subst_mtxp_4, substmtxfile);
  }
  
  /* get frequency file */
  if(freqfile != NULL) {
    read_frequencies_multi(freqfile, &aa_freqs, &aa_freqs_2, &aa_freqs_3, &aa_freqs_4);
  }
  
  /* get prior file */

  /* get replacement letters */
  if(replfile != NULL) {
    get_replacement_letters_multi(replfile, &replacement_letters);
  }
  else {
    replacement_letters.nr_alphabets = 0;
    replacement_letters.nr_rl_1 = 0;
    replacement_letters.nr_rl_2 = 0;
    replacement_letters.nr_rl_3 = 0;
    replacement_letters.nr_rl_4 = 0;
    replacement_letters.letters_1 = NULL;
    replacement_letters.probs_1 = NULL;
    replacement_letters.letters_2 = NULL;
    replacement_letters.probs_2 = NULL;
    replacement_letters.letters_3 = NULL;
    replacement_letters.probs_3 = NULL;
    replacement_letters.letters_4 = NULL;
    replacement_letters.probs_4 = NULL;
  }

  if((scoring_method == SUBST_MTX_PRODUCT || scoring_method == SUBST_MTX_DOT_PRODUCT ||
      scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR)
     && read_subst_mtx == NO) {
    printf("Error: No substitution matrix supplied\n");
    exit(0);
  }


  /* read hmm */
  /* get hmm from file */
  if(hmmfile != NULL) {
    hmmfiletype = readhmm_check(hmmfile);
    if(hmmfiletype == SINGLE_HMM) {
      readhmm(hmmfile, &hmm);
    }
    else if(hmmfiletype == MULTI_HMM) {
      readhmm_multialpha(hmmfile, &hmm);
    }
    hmm.subst_mtx = subst_mtxp;
    hmm.subst_mtx_2 = subst_mtxp_2;
    hmm.subst_mtx_3 = subst_mtxp_3;
    hmm.subst_mtx_4 = subst_mtxp_4;
  }
  else {
    /* cannot happen */
    printf("Internal error: hmmfile not found, should have been detected earlier\n");
  }

  hmm.replacement_letters = &replacement_letters;

  
  /* read template seq */
  /* allocate space for msa_seq_info structs */
  if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    msa_seq_infop_template = (struct msa_sequences_multi_s*)(malloc_or_die(1 * sizeof(struct msa_sequences_multi_s)));
  }
  else if(seq_format == FASTA || seq_format == STANDARD) {
    seq_info_template.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * 1);
    seq_info_template.nr_alphabets = hmm.nr_alphabets;
    seq_info_template.nr_seqs = 1;
    seq_info_template.longest_seq = 0;
    seq_info_template.shortest_seq = INT_MAX;
    seq_info_template.avg_seq_len = 0;
  }
  /* check sequence file for labels + check nolabel flag */
  if(seqfile_has_labels(template_seqfile) == YES) {
    use_labels = YES;
  }
  else {
    use_labels = NO;
  }
  if(args_info.nolabels_given) {
    use_labels = NO;
  }
  
  if(seq_format == FASTA) {
    if(hmm.nr_alphabets > 1) {
      printf("fasta is a one alphabet format only\n");
      exit(0);
    }
    else {
      get_sequence_fasta_multi(template_seqfile, &seq_info_template, 0);
      hmm.alphabet_type = DISCRETE;
      if(use_labels == YES) {
	get_labels_multi(template_seqfile, &seq_info_template, &hmm, 0);
      }
    }
  }
  else if(seq_format == STANDARD) {
    get_sequence_std_multi(template_seqfile, &seq_info_template, &hmm, 0);
    if(use_labels == YES) {
      get_labels_multi(template_seqfile, &seq_info_template, &hmm, 0);
      }
  }
  else if(seq_format == MSA_STANDARD) {
    if(use_lead_columns == YES) {
      get_sequences_msa_std_multi(template_seqfile, NULL, msa_seq_infop_template, &hmm,
				  lead_seq, &replacement_letters);
    }
    else {
      get_sequences_msa_std_multi(template_seqfile, NULL, msa_seq_infop_template, &hmm, -1, &replacement_letters);
    }
    
    /* get labels */
    if(use_labels == YES) {
      if(use_lead_columns == YES) {
	get_msa_labels_multi(template_seqfile, msa_seq_infop_template, &hmm);
      }
      else {
	get_msa_labels_all_columns_multi(template_seqfile, msa_seq_infop_template, &hmm);
      }
    }
  }
  else if(seq_format == PROFILE) {
    get_sequences_msa_prf_multi(template_seqfile, NULL, msa_seq_infop_template, &hmm);
  }
  fclose(template_seqfile);
  if(seq_format == STANDARD || seq_format == FASTA) {
    seq_info_template.avg_seq_len = ((int)(seq_info_template.avg_seq_len / seq_info_template.nr_seqs));
    //dump_labeled_seqs_multi(&seq_info_template);
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    
  }
  
  /* score template seq, get path */
  /* get seq_length */
  if(seq_format == STANDARD || seq_format == FASTA) {
    seq_len_template = get_seq_length(seq_info_template.seqs->seq_1);
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    if(use_lead_columns == YES) {
      seq_len_template = msa_seq_infop_template->nr_lead_columns;
    }
    else {
      seq_len_template = msa_seq_infop_template->msa_seq_length;
    }
  }
  if(seq_format == STANDARD || seq_format == FASTA) {
    viterbi_multi(&hmm, seq_info_template.seqs->seq_1, seq_info_template.seqs->seq_2, seq_info_template.seqs->seq_3,
		  seq_info_template.seqs->seq_4,
		  &viterbi_mtx, use_labels, multi_scoring_method);
    viterbi_score = (viterbi_mtx + get_mtx_index(seq_len_template+1, hmm.nr_v-1, hmm.nr_v))->prob;
    if(viterbi_score == DEFAULT) {
      viterbi_score = log10(0.0);
    }
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    msa_viterbi_multi(&hmm, msa_seq_infop_template, use_lead_columns, use_gap_shares, use_prior, &viterbi_mtx, use_labels, normalize,
		      scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
    viterbi_score = (viterbi_mtx + get_mtx_index(seq_len_template+1, hmm.nr_v-1, hmm.nr_v))->prob;
    if(viterbi_score == DEFAULT) {
      viterbi_score = log10(0.0);
    }
  }

  if(seq_format == STANDARD || seq_format == FASTA) {
    path_template = (int*)malloc_or_die((seq_len_template+1) * 2 * sizeof(int));
    path_incrementable_template = path_template;
    b = 0;
    get_viterbi_path_multi(viterbi_mtx + get_mtx_index(seq_len_template+1,hmm.nr_v-1,hmm.nr_v), &hmm, viterbi_mtx,
			   seq_len_template + 1, hmm.nr_v, path_incrementable_template, &b);
  }
  
  if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    path_template = (int*)malloc_or_die((seq_len_template+1) * 2 * sizeof(int));
    path_incrementable_template = path_template;
      b = 0;
      get_viterbi_path_multi(viterbi_mtx + get_mtx_index(seq_len_template+1,hmm.nr_v-1,hmm.nr_v), &hmm, viterbi_mtx,
			     seq_len_template + 1, hmm.nr_v, path_incrementable_template, &b);
  }
  
  /* print*/
#ifdef DEBUG_ALIGN
  printf("Path:\n");
  for(k = 0; k < seq_len_template; k++) {
    if(k % 60 == 0 && k != 0) {
      printf("\n");
    }
    printf("%d ", path_template[k]);
  }
  printf("\n");
#endif  

  /* deallocate stuff */
  free(viterbi_mtx);
  

  
  /* read seq 2 */
  /* allocate space for msa_seq_info structs */
  if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    msa_seq_infop_target = (struct msa_sequences_multi_s*)(malloc_or_die(1 * sizeof(struct msa_sequences_multi_s)));
  }
  else if(seq_format == FASTA || seq_format == STANDARD) {
    seq_info_target.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * 1);
    seq_info_target.nr_alphabets = hmm.nr_alphabets;
    seq_info_target.nr_seqs = 1;
    seq_info_target.longest_seq = 0;
    seq_info_target.shortest_seq = INT_MAX;
    seq_info_target.avg_seq_len = 0;
  }
  /* check sequence file for labels + check nolabel flag */
  if(seqfile_has_labels(target_seqfile) == YES) {
    use_labels = YES;
  }
  else {
    use_labels = NO;
  }
  if(args_info.nolabels_given) {
    use_labels = NO;
  }
  
  if(seq_format == FASTA) {
    if(hmm.nr_alphabets > 1) {
      printf("fasta is a one alphabet format only\n");
      exit(0);
    }
    else {
      get_sequence_fasta_multi(target_seqfile, &seq_info_target, 0);
      hmm.alphabet_type = DISCRETE;
      if(use_labels == YES) {
	get_labels_multi(target_seqfile, &seq_info_target, &hmm, 0);
      }
    }
  }
  else if(seq_format == STANDARD) {
    get_sequence_std_multi(target_seqfile, &seq_info_target, &hmm, 0);
    if(use_labels == YES) {
      get_labels_multi(target_seqfile, &seq_info_target, &hmm, 0);
      }
  }
  else if(seq_format == MSA_STANDARD) {
    if(use_lead_columns == YES) {
      get_sequences_msa_std_multi(target_seqfile, NULL, msa_seq_infop_target, &hmm,
				  lead_seq, &replacement_letters);
    }
    else {
      get_sequences_msa_std_multi(target_seqfile, NULL, msa_seq_infop_target, &hmm, -1, &replacement_letters);
    }
    
    /* get labels */
    if(use_labels == YES) {
      if(use_lead_columns == YES) {
	get_msa_labels_multi(target_seqfile, msa_seq_infop_target, &hmm);
      }
      else {
	get_msa_labels_all_columns_multi(target_seqfile, msa_seq_infop_target, &hmm);
      }
    }
  }
  else if(seq_format == PROFILE) {
    get_sequences_msa_prf_multi(target_seqfile, NULL, msa_seq_infop_target, &hmm);
  }
  fclose(target_seqfile);
  if(seq_format == STANDARD || seq_format == FASTA) {
    seq_info_target.avg_seq_len = ((int)(seq_info_target.avg_seq_len / seq_info_target.nr_seqs));
    //dump_labeled_seqs_multi(&seq_info_target);
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    
  }

  
  /* score seq 2, get path */
  /* get seq_length */
  if(seq_format == STANDARD || seq_format == FASTA) {
    seq_len_target = get_seq_length(seq_info_target.seqs->seq_1);
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    if(use_lead_columns == YES) {
      seq_len_target = msa_seq_infop_target->nr_lead_columns;
    }
    else {
      seq_len_target = msa_seq_infop_target->msa_seq_length;
    }
  }
  if(seq_format == STANDARD || seq_format == FASTA) {
    viterbi_multi(&hmm, seq_info_target.seqs->seq_1, seq_info_target.seqs->seq_2, seq_info_target.seqs->seq_3,
		  seq_info_target.seqs->seq_4,
		  &viterbi_mtx, use_labels, multi_scoring_method);
    viterbi_score = (viterbi_mtx + get_mtx_index(seq_len_target+1, hmm.nr_v-1, hmm.nr_v))->prob;
    if(viterbi_score == DEFAULT) {
      viterbi_score = log10(0.0);
    }
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    msa_viterbi_multi(&hmm, msa_seq_infop_target, use_lead_columns, use_gap_shares, use_prior, &viterbi_mtx, use_labels, normalize,
		      scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
    viterbi_score = (viterbi_mtx + get_mtx_index(seq_len_target+1, hmm.nr_v-1, hmm.nr_v))->prob;
    if(viterbi_score == DEFAULT) {
      viterbi_score = log10(0.0);
    }
  }

  if(seq_format == STANDARD || seq_format == FASTA) {
    path_target = (int*)malloc_or_die((seq_len_target+1) * 2 * sizeof(int));
    path_incrementable_target = path_target;
    b = 0;
    get_viterbi_path_multi(viterbi_mtx + get_mtx_index(seq_len_target+1,hmm.nr_v-1,hmm.nr_v), &hmm, viterbi_mtx,
			   seq_len_target + 1, hmm.nr_v, path_incrementable_target, &b);
  }
  
  if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    path_target = (int*)malloc_or_die((seq_len_target+1) * 2 * sizeof(int));
    path_incrementable_target = path_target;
      b = 0;
      get_viterbi_path_multi(viterbi_mtx + get_mtx_index(seq_len_target+1,hmm.nr_v-1,hmm.nr_v), &hmm, viterbi_mtx,
			     seq_len_target + 1, hmm.nr_v, path_incrementable_target, &b);
  }
  
  /* print*/
#ifdef DEBUG_ALIGN
  printf("Path:\n");
  for(k = 0; k < seq_len_target; k++) {
    if(k % 60 == 0 && k != 0) {
      printf("\n");
    }
    printf("%d ", path_target[k]);
  }
  printf("\n");
#endif

  /* deallocate stuff */
  free(viterbi_mtx);


  /* align paths so that maximum number of letters are emitted by the same state */
  
  //get seq names and send along
  fprintf(outfile, "Aligment based on same-state-emission\n");
  fprintf(outfile, "Target file: %s\n", args_info.target_arg);
  fprintf(outfile, "Template file: %s\n", args_info.template_arg);
  align_paths(&path_template, &path_target, seq_len_target, seq_len_template, outfile, seq_format);
  

  /* garbage collection */
  if(replfile != NULL) {
    if(replacement_letters.nr_rl_1 > 0) {
      free(replacement_letters.letters_1);
      free(replacement_letters.probs_1);
    }
    if(replacement_letters.nr_rl_2 > 0) {
      free(replacement_letters.letters_2);
      free(replacement_letters.probs_2);
    }
    if(replacement_letters.nr_rl_3 > 0) {
      free(replacement_letters.letters_3);
      free(replacement_letters.probs_3);
    }
    if(replacement_letters.nr_rl_4 > 0) {
      free(replacement_letters.letters_4);
      free(replacement_letters.probs_4);
    }
    fclose(replfile);
  }
  

  /* free allocated sequence memory for msa/profile and single seqs respectively */ 
  if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    free((*(msa_seq_infop_template)).msa_seq_1);
    if(hmm.nr_alphabets > 1) {
      free((*(msa_seq_infop_template)).msa_seq_2);
    }
    if(hmm.nr_alphabets > 2) {
      free((*(msa_seq_infop_template)).msa_seq_3);
    }
    if(hmm.nr_alphabets > 3) {
      free((*(msa_seq_infop_template)).msa_seq_4);
    }
    free((*(msa_seq_infop_template)).gaps);
    free((*(msa_seq_infop_template)).gap_shares);
    free((*(msa_seq_infop_template)).lead_columns_start);
    free(msa_seq_infop_template);
    
    free((*(msa_seq_infop_target)).msa_seq_1);
    if(hmm.nr_alphabets > 1) {
      free((*(msa_seq_infop_target)).msa_seq_2);
    }
    if(hmm.nr_alphabets > 2) {
      free((*(msa_seq_infop_target)).msa_seq_3);
    }
    if(hmm.nr_alphabets > 3) {
      free((*(msa_seq_infop_target)).msa_seq_4);
    }
    free((*(msa_seq_infop_target)).gaps);
    free((*(msa_seq_infop_target)).gap_shares);
    free((*(msa_seq_infop_target)).lead_columns_start);
    free(msa_seq_infop_target);
  }
  else if(seq_format == FASTA || seq_format == STANDARD) {
    free(((seq_info_template.seqs))->seq_1);
    if(hmm.nr_alphabets > 1) {
      free((seq_info_template.seqs)->seq_2);
    }
    if(hmm.nr_alphabets > 2) {
      free((seq_info_template.seqs)->seq_3);
      }
    if(hmm.nr_alphabets > 3) {
      free((seq_info_template.seqs)->seq_4);
    }
    free(seq_info_template.seqs);
    
    free(((seq_info_target.seqs))->seq_1);
    if(hmm.nr_alphabets > 1) {
      free((seq_info_target.seqs)->seq_2);
    }
    if(hmm.nr_alphabets > 2) {
      free((seq_info_target.seqs)->seq_3);
      }
    if(hmm.nr_alphabets > 3) {
      free((seq_info_target.seqs)->seq_4);
    }
    free(seq_info_target.seqs);
  }


  if(substmtxfile != NULL) {
    free(subst_mtxp);
    if(hmm.nr_alphabets > 1) {
      free(subst_mtxp_2);
    }
    if(hmm.nr_alphabets > 2) {
      free(subst_mtxp_3);
    }
    if(hmm.nr_alphabets > 3) {
      free(subst_mtxp_4);
    }
    fclose(substmtxfile);
  }
  
  if(freqfile != NULL) {
    free(aa_freqs);
    if(hmm.nr_alphabets > 1) {
      free(aa_freqs_2);
    }
    if(hmm.nr_alphabets > 2) {
      free(aa_freqs_3);
    }
    if(hmm.nr_alphabets > 3) {
      free(aa_freqs_4);
    }
    fclose(freqfile);
  }
  hmm_garbage_collection_multi(hmmfile, &hmm);


  free(path_template);
  free(path_target);
 
}

void align_paths(int **path_templatep, int **path_targetp, int seq_len_target, int seq_len_template, FILE *outfile,
		 int seq_format)
{
  int *path_template, *path_target;
  struct align_mtx_element_s *align_mtx;
  int i,j;
  int score_up, score_left, score_left_up, score;
  char last_up, last_left, last_left_up, last;
  int alignment_length, cur_alignment_pos;
  int target_pos, template_pos;
  struct alignment_s *alignment;
  

  alignment_length = seq_len_target + seq_len_template;

  path_template = *path_templatep;
  path_target = *path_targetp;
  
  align_mtx = (struct align_mtx_element_s*)(malloc_or_die((seq_len_target+1) * (seq_len_template+1) *
							  sizeof(struct align_mtx_element_s)));
  
  (align_mtx + get_mtx_index(0,0,seq_len_target+1))->score = 0;
  (align_mtx + get_mtx_index(0,0,seq_len_target+1))->last = 'e';
  for(i = 1; i <= seq_len_target; i++) {
    (align_mtx + get_mtx_index(0,i,seq_len_target+1))->score = 0;
    (align_mtx + get_mtx_index(0,i,seq_len_target+1))->last = 'l';
  }

  for(j = 1; j <= seq_len_template; j++) {
    (align_mtx + get_mtx_index(j,0,seq_len_target+1))->score = 0;
    (align_mtx + get_mtx_index(j,0,seq_len_target+1))->last = 'u';
  }
  

  for(j = 1; j <= seq_len_template; j++) {
    for(i = 1; i <= seq_len_target; i++) {
      score_up = (align_mtx + get_mtx_index(j-1,i,seq_len_target+1))->score;
      score_left = (align_mtx + get_mtx_index(j,i-1,seq_len_target+1))->score;
      score_left_up = (align_mtx + get_mtx_index(j-1,i-1,seq_len_target+1))->score;
      last_up = (align_mtx + get_mtx_index(j-1,i,seq_len_target+1))->last;
      last_left = (align_mtx + get_mtx_index(j,i-1,seq_len_target+1))->last;
      last_left_up = (align_mtx + get_mtx_index(j-1,i-1,seq_len_target+1))->last;
      if(*(path_target + i) == *(path_template + j)) {
	score_left_up += 1;
      }

      if(score_left_up > score_left && score_left_up > score_up) {
	/* left_up */
	(align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left_up;
	(align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a';
      }
      else if(score_up > score_left && score_up > score_left_up) {
	/* up */
	(align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	(align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
      }
      else if(score_left > score_left_up && score_left > score_up) {
	/* left */
	(align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left;
	(align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'l';
      }
      else if(score_left_up > score_left) {
	/* left_up + up */
	if(last_left_up == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a';
	}
	else if(last_left == -1) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'l';
	}
	else if(last_left == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'l';
	}
	else {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a';
	}
      }
      else if(score_left_up > score_up) {
	/* left_up + left */
	if(last_left_up == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a' ;
	}
	else if(last_up == 1) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
	}
	else if(last_up == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
	}
	else {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a';
	}	
      }
      else if(score_up > score_left_up) {
	/* up + left */
	if(last_up == 1) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
	}
	else if(last_left == 1) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'l';
	}
	else if(last_up == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
	}
	else if(last_left == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a';
	}
	else {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
	}
      }
      else {
	/* all scores equal */
	if(last_left_up == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a';
	}
	else if(last_up == 1) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
	}
	else if(last_left == 1) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'l';
	}
	else if(last_up == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'u';
	}
	else if(last_left == 0) {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'l';
	}
	else {
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->score = score_left_up;
	  (align_mtx + get_mtx_index(j,i,seq_len_target+1))->last = 'a';
	}
      }
    }
  }
#ifdef DEBUG_ALIGN
  dump_align_matrix(seq_len_template+1, seq_len_target+1, align_mtx);
#endif
  
  target_pos = seq_len_target;
  template_pos = seq_len_template;
  
  alignment = (struct alignment_s*)(malloc_or_die(alignment_length * sizeof(struct alignment_s)));
  cur_alignment_pos = 0;
  while((align_mtx + get_mtx_index(template_pos, target_pos, seq_len_target+1))->last != 'e') {
    /* get letter combination */
    score = (align_mtx + get_mtx_index(template_pos, target_pos, seq_len_target+1))->score;
    last = (align_mtx + get_mtx_index(template_pos, target_pos, seq_len_target+1))->last;
    if(last == 'a') {
      (alignment + cur_alignment_pos)->target_pos = target_pos;
      (alignment + cur_alignment_pos)->template_pos = template_pos;
      if(seq_format == FASTA || seq_format == STANDARD) {
	strncpy((alignment + cur_alignment_pos)->target_letter, seq_info_target.seqs->seq_1 + target_pos-1, 5);
	strncpy((alignment + cur_alignment_pos)->template_letter, seq_info_template.seqs->seq_1 + template_pos-1, 5);
      }
      else if(seq_format == PROFILE || seq_format == MSA_STANDARD) {
	strncpy((alignment + cur_alignment_pos)->target_letter,
		(msa_seq_infop_target->msa_seq_1 + ((target_pos-1) * (hmm.a_size+1)))->query_letter, 5);
	strncpy((alignment + cur_alignment_pos)->template_letter,
		(msa_seq_infop_template->msa_seq_1 + ((template_pos-1) * (hmm.a_size+1)))->query_letter, 5);
      }
      target_pos--;
      template_pos--;
    }
    else if(last == 'u') {
      (alignment + cur_alignment_pos)->target_pos = -1;
      (alignment + cur_alignment_pos)->template_pos = template_pos;
      if(seq_format == FASTA || seq_format == STANDARD) {
	strncpy((alignment + cur_alignment_pos)->template_letter,
		(msa_seq_infop_target->msa_seq_1 + ((target_pos-1) * (hmm.a_size+1)))->query_letter, 5);
	(alignment + cur_alignment_pos)->target_letter[0] = '-';
	(alignment + cur_alignment_pos)->target_letter[1] = '\0';
      }
      else if(seq_format == PROFILE || seq_format == MSA_STANDARD) {
	strncpy((alignment + cur_alignment_pos)->template_letter,
		(msa_seq_infop_template->msa_seq_1 + ((template_pos-1) * (hmm.a_size+1)))->query_letter, 5);
	(alignment + cur_alignment_pos)->target_letter[0] = '-';
	(alignment + cur_alignment_pos)->target_letter[1] = '\0';
      }
      template_pos--;
    }
    else if(last == 'l') {
      (alignment + cur_alignment_pos)->target_pos = target_pos;
      (alignment + cur_alignment_pos)->template_pos = -1;
      if(seq_format == FASTA || seq_format == STANDARD) {
	strncpy((alignment + cur_alignment_pos)->target_letter, seq_info_target.seqs->seq_1 + target_pos-1, 5);
	(alignment + cur_alignment_pos)->template_letter[0] = '-';
	(alignment + cur_alignment_pos)->template_letter[1] = '\0';
      }
      else if(seq_format == PROFILE || seq_format == MSA_STANDARD) {
	strncpy((alignment + cur_alignment_pos)->target_letter,
		(msa_seq_infop_target->msa_seq_1 + ((target_pos-1) * (hmm.a_size+1)))->query_letter, 5);
	(alignment + cur_alignment_pos)->template_letter[0] = '-';
	(alignment + cur_alignment_pos)->template_letter[1] = '\0';
      }
      target_pos--;
    }
    else {

    }
    cur_alignment_pos++;
    
  }
#ifdef DEBUG_ALIGN
  for(i = cur_alignment_pos - 1; i >= 0; i--) {
    printf("%d\t%d\n", (alignment + i)->target_pos, (alignment + i)->template_pos);
  }
  for(i = cur_alignment_pos - 1; i >= 0; i--) {
    printf("%s\t%s\n", (alignment + i)->target_letter, (alignment + i)->template_letter);
  }
#endif
  
  
  fprintf(outfile, "target\ttemplate\n");
  for(i = cur_alignment_pos - 1; i >= 0; i--) {
    fprintf(outfile, "%s\t%s\n", (alignment + i)->target_letter, (alignment + i)->template_letter);
  }
  fclose(outfile);
  

  free(alignment);
  free(align_mtx);
}
