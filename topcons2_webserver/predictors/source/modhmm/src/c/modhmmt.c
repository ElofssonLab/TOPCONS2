#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include "structs.h"
#include "funcs.h"
#include "cmdline_modhmmt.h"

#define MAX_LINE 4000
#define MAX_SEQS 4000

#define FIRST_SEQ 1

extern int verbose;
extern int user_defined_emission_score;

/* memory for transition and emission matrices, from_vertex array and 
   to_vertex_array will be allocated in readhmm, but must be freed here */
static struct hmm_multi_s hmm;
static struct msa_sequences_multi_s *msa_seq_infop;
static struct msa_sequences_multi_s *msa_negseq_infop;
static struct sequences_multi_s seq_info;
static struct sequences_multi_s negseq_info;
static struct replacement_letter_multi_s replacement_letters;
long double *subst_mtxp;
long double *subst_mtxp_2;
long double *subst_mtxp_3;
long double *subst_mtxp_4;
long double *aa_freqs;
long double *aa_freqs_2;
long double *aa_freqs_3;
long double *aa_freqs_4;


int main(int argc, char* argv[])
{
  int i;
  FILE *hmmfile, *outfile, *seqfile, *negseqfile, *replfile, *seqnamefile, *negseqnamefile, *substmtxfile, *freqfile;
  long double d;
  int seq_format;
  int opt_alpha;
  int use_gap_shares;
  int use_lead_columns;
  int lead_seq, lead_negseq;
  int use_labels;
  int annealing;
  int use_transition_pseudo_counts, use_emission_pseudo_counts;
  int nr_seqs, seq_counter, nr_read_seqs;
  int nr_negseqs, negseq_counter, nr_read_negseqs;
  char seq_name[200];
  char negseq_name[200];
  int normalize;
  int scoring_method;
  int read_subst_mtx;
  int use_nr_occurences;
  int multi_scoring_method;
  int training_method;
  int hmmfiletype;
  int use_prior;
  int transonly;
  int emissonly;
 
  struct gengetopt_args_info args_info;

  seq_format = STANDARD;
  opt_alpha = 0;
  hmmfile = NULL;
  outfile = NULL;
  seqfile = NULL;
  replfile = NULL;
  seqnamefile = NULL;
  negseqnamefile = NULL;
  substmtxfile = NULL;
  freqfile = NULL;
  use_gap_shares = YES;
  use_lead_columns = YES;
  lead_seq = FIRST_SEQ;
  use_labels = NO;
  annealing = NO;
  use_transition_pseudo_counts = NO;
  use_emission_pseudo_counts = NO;
  normalize = NO;
  scoring_method = SJOLANDER;
  read_subst_mtx = NO;
  use_nr_occurences = NO;
  multi_scoring_method = JOINT_PROB;
  subst_mtxp = NULL;
  subst_mtxp_2 = NULL;
  subst_mtxp_3 = NULL;
  subst_mtxp_4 = NULL;
  aa_freqs = NULL;
  aa_freqs_2 = NULL;
  aa_freqs_3 = NULL;
  aa_freqs_4 = NULL;
  training_method = BW_STD;
  use_prior = YES;
  emissonly = NO;
  transonly = NO;
  

  /* parse command line */

  if(cmdline_parser(argc, argv, &args_info) != 0) {
    exit(1);
  }

  /* compulsory options */
  if(args_info.hmminfile_given) {
    if((hmmfile = fopen(args_info.hmminfile_arg, "r")) == NULL) {
      perror(args_info.hmminfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading model file\n",args_info.hmminfile_arg);
    }
  }
  if(args_info.seqnamefile_given) {
    if((seqnamefile = fopen(args_info.seqnamefile_arg, "r")) == NULL) {
      perror(args_info.seqnamefile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading sequence names\n",args_info.seqnamefile_arg);
    }
  }
  if(args_info.negseqnamefile_given) {
    if((negseqnamefile = fopen(args_info.negseqnamefile_arg, "r")) == NULL) {
      perror(args_info.negseqnamefile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading sequence names\n",args_info.negseqnamefile_arg);
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
  if(args_info.alg_given) {
    if(strcmp(args_info.alg_arg, "bw") == 0) {
      training_method = BW_STD;
    }
    else if(strcmp(args_info.alg_arg, "cml") == 0) {
      training_method = CML_STD;
    }
    else if(strcmp(args_info.alg_arg, "disc") == 0) {
      training_method = DISC_STD;
    }
    else {
      printf("Incorrect training method option: %s\n", args_info.alg_arg);
      exit(0);
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
      use_lead_columns = NO;
    }
    else {
      lead_seq = atoi(args_info.usecolumns_arg);
      use_lead_columns = YES;
      if(lead_seq <= 0) {
	printf("Incorrect use-column option: %s\n", args_info.usecolumns_arg);
	exit(0);
      }
    }
  }

  if(args_info.optalpha_given) {
    opt_alpha = atoi(args_info.optalpha_arg);
    if(opt_alpha <= 0 || opt_alpha > 4) {
      printf("Incorrect optalpha option: %s\n", args_info.optalpha_arg);
      exit(0);
    }
  }

  /* flags */
  if(args_info.nolabels_given) {
    /* checked after seqread */
  }
  if(args_info.noprior_given) {
    /* checked after hmm read */
  }
  if(args_info.tpcounts_given) {
    use_transition_pseudo_counts = YES;
  }
  if(args_info.epcounts_given) {
    use_emission_pseudo_counts = YES;
  }
  if(args_info.verbose_given) {
    verbose = YES;
  }
  if(args_info.transonly_given) {
    transonly = YES;;
  }
  if(args_info.emissonly_given) {
    emissonly = YES;
  }
  
  /* read subst mtx */
  if(substmtxfile != NULL) {  
    read_subst_matrix_multi(&subst_mtxp, &subst_mtxp_2, &subst_mtxp_3, &subst_mtxp_4, substmtxfile);
  }
  
  /* get frequency file */
  if(freqfile != NULL) {
    read_frequencies_multi(freqfile, &aa_freqs, &aa_freqs_2, &aa_freqs_3, &aa_freqs_4);
  }
  
  if((scoring_method == SUBST_MTX_PRODUCT || scoring_method == SUBST_MTX_DOT_PRODUCT ||
      scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR)
     && read_subst_mtx == NO) {
    printf("Error: No substitution matrix supplied\n");
    exit(0);
  }
   
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

  /* check if priors are present and if noprior option is set */
  if(args_info.noprior_given) {
    use_prior = NO;
  }
  else if(hmm.nr_ed > 0) {
    use_prior = YES;
  }
  
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
  hmm.replacement_letters = &replacement_letters;

  /* get training seqences */
  /* Note: memory for the sequences will be allocated in get_sequences method, but
   * must be freed here */
  
  /* count nr of sequences */
  nr_seqs = 0;
  while(fgets(seq_name, MAX_LINE, seqnamefile) != NULL) {
    nr_seqs++;
  }
  rewind(seqnamefile);

  /* count nr of negsequences */
  nr_negseqs = 0;
  if(args_info.negseqnamefile_given) {
    while(fgets(negseq_name, MAX_LINE, negseqnamefile) != NULL) {
      nr_negseqs++;
    }
    rewind(negseqnamefile);
  }
  /* allocate space for msa_seq_info structs */
  if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    msa_seq_infop = (struct msa_sequences_multi_s*)(malloc_or_die(nr_seqs * sizeof(struct msa_sequences_multi_s)));
    
  }
  else if(seq_format == FASTA || seq_format == STANDARD) {
    seq_info.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * nr_seqs);
    seq_info.nr_alphabets = hmm.nr_alphabets;
    seq_info.nr_seqs = nr_seqs;
    seq_info.longest_seq = 0;
    seq_info.shortest_seq = INT_MAX;
    seq_info.avg_seq_len = 0;
  }
  /* read all the sequences */
  seq_counter = 0;
  nr_read_seqs = 0;
  while(fgets(seq_name, MAX_LINE, seqnamefile) != NULL) {
    //printf("reading seq %s", seq_name);
    seq_name[strlen(seq_name)-1] = '\0';
    if((seqfile = fopen(seq_name, "r")) == NULL) {
      perror(seq_name);
      continue;
    }
    
    /* check sequence file for labels + check nolabel flag */
    if(seqfile_has_labels(seqfile) == YES) {
      use_labels = YES;
    }
    else if(seq_format == PROFILE) {
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
	get_sequence_fasta_multi(seqfile, &seq_info, nr_read_seqs);
	hmm.alphabet_type = DISCRETE;
	if(use_labels == YES) {
	  get_labels_multi(seqfile, &seq_info, &hmm, nr_read_seqs);
	}
      }
    }
    else if(seq_format == STANDARD) {
      get_sequence_std_multi(seqfile, &seq_info, &hmm, nr_read_seqs);
      if(use_labels == YES) {
	get_labels_multi(seqfile, &seq_info, &hmm, nr_read_seqs);
      }
    }
    else if(seq_format == MSA_STANDARD) {
      if(use_lead_columns == YES) {
	get_sequences_msa_std_multi(seqfile, NULL, msa_seq_infop + seq_counter, &hmm,
				    lead_seq, &replacement_letters);
      }
      else {
	get_sequences_msa_std_multi(seqfile, NULL, msa_seq_infop + seq_counter, &hmm, -1, &replacement_letters);
      }
      
      /* get labels */
      if(use_labels == YES) {
	if(use_lead_columns == YES) {
	  get_msa_labels_multi(seqfile, msa_seq_infop + seq_counter, &hmm);
	}
	else {
	  get_msa_labels_all_columns_multi(seqfile, msa_seq_infop + seq_counter, &hmm);
	}
      }
      else {
      }
    }
    else if(seq_format == PROFILE) {
      get_sequences_msa_prf_multi(seqfile, NULL, msa_seq_infop + seq_counter, &hmm);
    }
    nr_read_seqs++;
    seq_counter++;
    fclose(seqfile);
  }
  if(seq_format == STANDARD || seq_format == FASTA) {
    seq_info.avg_seq_len = ((int)(seq_info.avg_seq_len / seq_info.nr_seqs));
    //dump_labeled_seqs_multi(&seq_info);
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    
  }
 

  /* read all the negative sequences */
  if(args_info.negseqnamefile_given) {
    /* allocate space for msa_seq_info structs */
    if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      msa_negseq_infop = (struct msa_sequences_multi_s*)(malloc_or_die(nr_negseqs * sizeof(struct msa_sequences_multi_s)));
      
    }
    else if(seq_format == FASTA || seq_format == STANDARD) {
      negseq_info.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * nr_negseqs);
      negseq_info.nr_alphabets = hmm.nr_alphabets;
      negseq_info.nr_seqs = nr_negseqs;
      negseq_info.longest_seq = 0;
      negseq_info.shortest_seq = INT_MAX;
      negseq_info.avg_seq_len = 0;
    }
    
    negseq_counter = 0;   
    nr_read_negseqs = 0;  
    while(fgets(negseq_name, MAX_LINE, negseqnamefile) != NULL) {
      //printf("reading negseq %s", negseq_name);
      negseq_name[strlen(negseq_name)-1] = '\0';
      if((negseqfile = fopen(negseq_name, "r")) == NULL) { 
	perror(negseq_name);
	continue;
      }
      
      /* check sequence file for labels + check nolabel flag */
      if(seqfile_has_labels(negseqfile) == YES) {
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
	get_sequence_fasta_multi(negseqfile, &negseq_info, nr_read_negseqs);
	hmm.alphabet_type = DISCRETE;
	if(use_labels == YES) {
	  get_labels_multi(negseqfile, &negseq_info, &hmm, nr_read_negseqs);
	}
	}
      }
      else if(seq_format == STANDARD) {
	get_sequence_std_multi(negseqfile, &negseq_info, &hmm, nr_read_negseqs);
	if(use_labels == YES) {
	  get_labels_multi(negseqfile, &negseq_info, &hmm, nr_read_negseqs);
	}
      }
      else if(seq_format == MSA_STANDARD) {
	if(use_lead_columns == YES) {
	  get_sequences_msa_std_multi(negseqfile, NULL, msa_negseq_infop + negseq_counter, &hmm,
				      lead_negseq, &replacement_letters);
	}
	else {
	  get_sequences_msa_std_multi(negseqfile, NULL, msa_negseq_infop + negseq_counter, &hmm, -1, &replacement_letters);
	}
	

	/* get labels */
	if(use_labels == YES) {
	  if(use_lead_columns == YES) {
	    get_msa_labels_multi(negseqfile, msa_negseq_infop + negseq_counter, &hmm);
	  }
	  else {
	    get_msa_labels_all_columns_multi(negseqfile, msa_negseq_infop + negseq_counter, &hmm);
	  }
	}
	else {
	}
      }
      else if(seq_format == PROFILE) {
	get_sequences_msa_prf_multi(negseqfile, NULL, msa_negseq_infop + negseq_counter, &hmm);
      }
      nr_read_negseqs++;
      negseq_counter++;
      fclose(negseqfile);
    }
    if(seq_format == STANDARD || seq_format == FASTA) {
      negseq_info.avg_seq_len = ((int)(negseq_info.avg_seq_len / negseq_info.nr_seqs));
      //dump_labeled_seqs_multi(&seq_info);
    }
    else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      
    }
  }

  /* do training */
  if(nr_read_seqs > 0) {
    if(training_method == BW_STD) {
      /* do baum welch training */
      if(seq_format == STANDARD || seq_format == FASTA) {
	baum_welch_dirichlet_multi(&hmm, seq_info.seqs, seq_info.nr_seqs, annealing, use_labels, use_transition_pseudo_counts,
				   use_emission_pseudo_counts, multi_scoring_method, use_prior, transonly, emissonly);
      }
      else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	msa_baum_welch_dirichlet_multi(&hmm, msa_seq_infop, nr_seqs, annealing, use_gap_shares, use_lead_columns, use_labels,
				       use_transition_pseudo_counts, use_emission_pseudo_counts, normalize, scoring_method,
				       use_nr_occurences, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4,
				       use_prior, transonly, emissonly);
      }
    }
    else if(training_method == DISC_STD && args_info.negseqnamefile_given) {
      if(seq_format == STANDARD || seq_format == FASTA) {
	//discriminative_baum_welch_dirichlet_multi(&hmm, seq_info.seqs, seq_info.nr_seqs, annealing, use_labels,
	//use_transition_pseudo_counts,
	//use_emission_pseudo_counts, multi_scoring_method, use_prior, negseq_info.seqs,
	//negseq_info.nr_seqs, opt_alpha);
      }
      else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	discriminative_msa_baum_welch_dirichlet_multi(&hmm, msa_seq_infop, nr_seqs, annealing, use_gap_shares, use_lead_columns, use_labels,
						      use_transition_pseudo_counts, use_emission_pseudo_counts, normalize, scoring_method,
						      use_nr_occurences, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3,
						      aa_freqs_4, use_prior, msa_negseq_infop, nr_negseqs, opt_alpha, transonly, emissonly);
	
      }
      
    }
    else if(training_method == DISC_STD && !(args_info.negseqnamefile_given)) {
      printf("Error: cannot do discriminative training without negative training sequences\n");
      exit(0);
    }
    else if(training_method == CML_STD && use_labels == YES) {
      /* do cml training */
      if(seq_format == STANDARD || seq_format == FASTA) {
	extended_baum_welch_dirichlet_multi(&hmm, seq_info.seqs, seq_info.nr_seqs, annealing, use_labels,
					   use_transition_pseudo_counts,
					    use_emission_pseudo_counts, multi_scoring_method, use_prior, transonly, emissonly);
      }
      else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	extended_msa_baum_welch_dirichlet_multi(&hmm, msa_seq_infop, nr_seqs, annealing, use_gap_shares, use_lead_columns, use_labels,
						use_transition_pseudo_counts, use_emission_pseudo_counts, normalize, scoring_method,
						use_nr_occurences, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3,
						aa_freqs_4, use_prior, transonly, emissonly);
      }
    }
    else if(training_method == CML_STD && use_labels == NO) {
      printf("Error: cannot do CML-training on unlabeled sequences\n");
      exit(0);
    }
  }
  else {
    printf("Error: no sequences read: cannot train\n");
    exit(0);
  }
   
  /* save trained hmm */
  if(hmmfiletype == MULTI_HMM) {
    if(outfile != NULL) {
      savehmm_multialpha(outfile, &hmm);
    }
    else {
      outfile = fopen("a.hmg", "w");
      savehmm_multialpha(outfile, &hmm);
    }
  }
  else if(hmmfiletype == SINGLE_HMM) {
    if(outfile != NULL) {
      savehmm(outfile, &hmm);
    }
    else {
      outfile = fopen("a.hmg", "w");
      savehmm(outfile, &hmm);
    }
  }

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
    for(i = 0; i < nr_read_seqs; i++) {
      free((*(msa_seq_infop + i)).msa_seq_1);
      if(hmm.nr_alphabets > 1) {
	free((*(msa_seq_infop + i)).msa_seq_2);
      }
      if(hmm.nr_alphabets > 2) {
	free((*(msa_seq_infop + i)).msa_seq_3);
      }
      if(hmm.nr_alphabets > 3) {
	free((*(msa_seq_infop + i)).msa_seq_4);
      }
      free((*(msa_seq_infop + i)).gaps);
      free((*(msa_seq_infop + i)).gap_shares);
      free((*(msa_seq_infop + i)).lead_columns_start);
    }
    free(msa_seq_infop);
    if(args_info.negseqnamefile_given) {
      for(i = 0; i < nr_read_negseqs; i++) {
	free((*(msa_negseq_infop + i)).msa_seq_1);
	if(hmm.nr_alphabets > 1) {
	  free((*(msa_negseq_infop + i)).msa_seq_2);
	}
	if(hmm.nr_alphabets > 2) {
	  free((*(msa_negseq_infop + i)).msa_seq_3);
	}
	if(hmm.nr_alphabets > 3) {
	  free((*(msa_negseq_infop + i)).msa_seq_4);
	}
	free((*(msa_negseq_infop + i)).gaps);
	free((*(msa_negseq_infop + i)).gap_shares);
	free((*(msa_negseq_infop + i)).lead_columns_start);
      }
      free(msa_negseq_infop);
    }
  }
  else if(seq_format == FASTA || seq_format == STANDARD) {
    for(i = 0; i < nr_read_seqs; i++) {
      free(((seq_info.seqs) + i)->seq_1);
      if(hmm.nr_alphabets > 1) {
	free((seq_info.seqs + i)->seq_2);
      }
      if(hmm.nr_alphabets > 2) {
	free((seq_info.seqs + i)->seq_3);
      }
      if(hmm.nr_alphabets > 3) {
	free((seq_info.seqs + i)->seq_4);
      }
    }
    free(seq_info.seqs);
    if(args_info.negseqnamefile_given) {
      for(i = 0; i < nr_read_negseqs; i++) {
      free(((negseq_info.seqs) + i)->seq_1);
      if(hmm.nr_alphabets > 1) {
	free((negseq_info.seqs + i)->seq_2);
      }
      if(hmm.nr_alphabets > 2) {
	free((negseq_info.seqs + i)->seq_3);
      }
      if(hmm.nr_alphabets > 3) {
	free((negseq_info.seqs + i)->seq_4);
      }
    }
    free(negseq_info.seqs);
    }

  }
  
  fclose(outfile);


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
  return 0;
}


