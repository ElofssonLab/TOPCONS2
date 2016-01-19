#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>




#include "structs.h"
#include "funcs.h"

//#define DEBUG_LABELING_UPDATE
//#define DEBUG_DEALLOCATE_LABELINGS


#define REST_LETTER_INDEX 0.5

#define V_LIST_END  -99
#define V_LIST_NEXT  -9

long double pos_spec_DG(char aa, int i, int helix_length);


long double get_single_gaussian_statescore(long double mu, long double sigma_square, long double letter)
{

  long double res;

  if(sigma_square <= 0.0) {
    return 0.0;
  }
  else {
    res = exp(0.0 - (pow((letter-mu),2) / (2.0 * sigma_square))) / sqrt(sigma_square * 2.0 * 3.141592655);
    //printf("mu = %Lf, letter = %Lf, sigma_square = %Lf,  res = %Lf\n", mu, letter, sigma_square, res);
    return res;
  }
}

long double get_dp_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			 int p, long double *emissions,  int vertex, int normalize, long double *gap_shares)
{
  
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  
  t_res_3 = 0.0;
  /* scoring using dot-product method */
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share;
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      if( *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
	  (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share) < 0.0) {
	printf("seqpos = %d\n", p);
	printf("vertex = %d\n", vertex);
	printf("a_index = %d\n", a_index);
	printf("a_size = %d\n", a_size);
	printf("emission_p = %Lf\n", *(emissions + get_mtx_index(vertex, a_index, a_size)));
	printf("seq_share = %Lf\n", (msa_seq + get_mtx_index(p,a_index, a_size+1))->share);
	printf("gap_share = %Lf\n", (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share));
      }


      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
	(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share;
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }

  if(t_res_3 < 0.0) {
    printf("t_res_3 = %Lf\n", t_res_3);
    printf("Error: got strange dot product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


long double get_dp_picasso_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
				 int p, long double *emissions,  int vertex, int normalize, long double *gap_shares, long double *aa_freqs)
{
  
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  
  t_res_3 = 0.0;
  /* scoring using dot-product method */
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share /  *(aa_freqs + a_index);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
	((1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share) *  *(aa_freqs + a_index));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share /  *(aa_freqs + a_index);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %Lf\n", state_normalizer);
    printf("seq_normalizer = %Lf\n", seq_normalizer);
#endif     
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }

  if(t_res_3 < 0.0) {
    printf("t_res_3 = %Lf\n", t_res_3);
    printf("Error: got strange dot product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


long double get_sjolander_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
				int p, long double *emissions,  int vertex, int normalize, long double *gap_shares)
{
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;
  long double emiss_prob;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;

  t_res_3 = 1.0;
  /* scoring using sjolander score method */

  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      emiss_prob = *(emissions + (vertex * a_size + a_index));
      if(emiss_prob == 0.0) {
	//emiss_prob = 1.0;
      }
      t_res_3 *= pow(emiss_prob, 
		     (msa_seq + (p * (a_size+1) + a_index))->prior_share);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
     
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      emiss_prob = *(emissions + get_mtx_index(vertex, a_index, a_size));
      if(emiss_prob == 0.0) {
	//emiss_prob = 1.0;
      }
      t_res_3 *= pow(emiss_prob,
		     (msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
		     (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      emiss_prob = *(emissions + get_mtx_index(vertex, a_index, a_size));
      if(emiss_prob == 0.0) {
	//emiss_prob = 1.0;
      }
      t_res_3 *= pow(emiss_prob,
		     (msa_seq + get_mtx_index(p,a_index, a_size+1))->share);
      
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }

  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %Lf\n", state_normalizer);
    printf("seq_normalizer = %Lf\n", seq_normalizer);
#endif     
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }
  if(t_res_3 < 0.0) {
    printf("Error: got strange geometric mean state result value\n");
    printf("t_res_3 = %Lf\n", t_res_3);
    exit(0);
  }
  return t_res_3;
}



long double get_sjolander_reversed_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					 int p, long double *emissions,  int vertex, int normalize, long double *gap_shares)
{
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;

  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  

  t_res_3 = 1.0;
  /* scoring using sjolander score method */
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share != 0.0) {
      	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share,
      		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
		       (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share),
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share,
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %Lf\n", state_normalizer);
    printf("seq_normalizer = %Lf\n", seq_normalizer);
#endif     
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }
  if(t_res_3 < 0.0) {
    printf("Error: got strange geometric mean state result value\n");
    printf("t_res_3 = %Lf\n", t_res_3);
    exit(0);
  }
  
  return t_res_3;
}


long double get_picasso_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			      int p, long double *emissions,  int vertex, int normalize, long double *gap_shares, long double *aa_freqs)
{
  
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  

  t_res_3 = 1.0;
  /* scoring using picasso-product method */

  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share / *(aa_freqs + a_index),
		       *(emissions + (vertex * a_size + a_index)));
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow(((msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
			(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share)) / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
    }
    else {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0 &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
    }
  }
  
  if(t_res_3 < 0.0 || t_res_3 > 1000000000000.0) {
    printf("Error: got strange picasso product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


long double get_picasso_sym_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
				  int p, long double *emissions,  int vertex, int normalize, long double *gap_shares, long double *aa_freqs)
{
  
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  

  t_res_3 = 1.0;
  /* scoring using picasso-product method */

  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share / *(aa_freqs + a_index),
		       *(emissions + (vertex * a_size + a_index))) *
	  pow(*(emissions + (vertex * a_size + a_index)) / *(aa_freqs + a_index),
	      (msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow(((msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
			(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share)) / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size))) *
	  pow(*(emissions + (vertex * a_size + a_index)) / *(aa_freqs + a_index),
	      (msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
	      (1.0 - (msa_seq + get_mtx_index(p,a_size, a_size+1))->share));
      }
    }
    else {
      //printf(" *(aa_freqs + a_index) = %Lf\n",  *(aa_freqs + a_index));
      //printf("(msa_seq + (p * (a_size+1) + a_index))->share = %Lf\n", (msa_seq + get_mtx_index(p,a_index, a_size+1))->share);
      //printf("*(emissions + get_mtx_index(vertex, a_index, a_size)) = %Lf\n",*(emissions + get_mtx_index(vertex, a_index, a_size))); 
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0 &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size))) *
	  pow(*(emissions + (vertex * a_size + a_index)) / *(aa_freqs + a_index),
	      (msa_seq + get_mtx_index(p,a_index, a_size+1))->share);
      }
    }
  }

  if(t_res_3 < 0.0 || t_res_3 > 1000000000000.0) {
    printf("Error: got strange picasso product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


long double get_subst_mtx_product_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					int p, long double *emissions, int vertex, long double *subst_mtx)
{
  int a_index, a_index2;
  long double t_res_3;
  
  t_res_3 = 0.0;
  for(a_index = 0; a_index < a_size; a_index++) {
    for(a_index2 = 0; a_index2 < a_size; a_index2++) {
      if(use_gap_shares == YES) {
	t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	  (msa_seq + get_mtx_index(p,a_index2, a_size+1))->share /
	  (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share) *
	  *(subst_mtx + (a_index * a_size + a_index2));
      }
      else {
	t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	  (msa_seq + get_mtx_index(p,a_index2, a_size+1))->share /
	  *(subst_mtx + (a_index * a_size + a_index2));
      }
    }
  }
  
  if(t_res_3 < 0.0) {
    printf("Error: got strange subst mtx product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


long double get_subst_mtx_dot_product_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					    int p, long double *emissions,  int vertex, int normalize, long double *gap_shares,
					    int query_index, long double *subst_mtx)
{
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  
  t_res_3 = 0.0;
  /* scoring using subst_mtx_dot-product method */
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	  exit(0);
      }
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size)) /
	(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
      }
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
      }
    }
  }
  
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
    subst_mtx_normalizer = sqrt(subst_mtx_normalizer);
#ifdef DEBUG
    printf("state_normalizer = %Lf\n", state_normalizer);
      printf("seq_normalizer = %Lf\n", seq_normalizer);
#endif      
      if(t_res_3 != 0.0) {
	t_res_3 = t_res_3 / (seq_normalizer * state_normalizer * subst_mtx_normalizer);
      }
  }

  if(t_res_3 < 0.0) {
    printf("Error: got strange subst mtx dot product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


long double get_subst_mtx_dot_product_prior_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
						  int p, long double *emissions,  int vertex, int normalize, long double *gap_shares,
						  int query_index, long double *subst_mtx)
{
  long double seq_normalizer;
  long double state_normalizer;
  long double subst_mtx_normalizer;
  int a_index;
  long double t_res_3;
  long double rest_share, default_share;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  default_share = 1.0 / (long double)(a_size);

  t_res_3 = 0.0;
  
  /* scoring using dot-product method */
  rest_share = 1.0;
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)) != 0.0) {
	rest_share = rest_share - (msa_seq + (p * (a_size+1) + a_index))->prior_share;
	if(normalize == YES) {
	  seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	  state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	  subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
	}
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size)) /
	(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share);
      if(*(subst_mtx + get_mtx_index(a_index, a_index, a_size)) != 0.0) {
	rest_share = rest_share - (msa_seq + (p * (a_size+1) + a_index))->share / 
	  (1.0 - *(gap_shares + p));
	if(normalize == YES) {
	  seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
				(1.0 - *(gap_shares + p)), 2);
	  state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	  subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
	}
      }
      
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(*(subst_mtx + get_mtx_index(a_index, a_index, a_size)) != 0.0) {
	rest_share = rest_share - (msa_seq + (p * (a_size+1) + a_index))->share;
	if(normalize == YES) {
	  seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	  state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	  subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
	}
      }
    }
  }
  if(rest_share < 0.0) {
    rest_share = 0.0;
  }
  t_res_3 += default_share * rest_share;
  seq_normalizer += pow(rest_share, 2);
  state_normalizer += pow(default_share, 2);
  
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
    subst_mtx_normalizer = sqrt(subst_mtx_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %Lf\n", state_normalizer);
    printf("seq_normalizer = %Lf\n", seq_normalizer);
#endif
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer * subst_mtx_normalizer);
    }
  }

  if(t_res_3 < 0.0) {
    printf("Error: got strange subst mtx dot product prior state result value\n");
    exit(0);
  }
  
  return t_res_3;
}



/************************************* add to E methods *********************************************/
void add_to_E_continuous(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
			 int k, int a_size, long double *emissions)
{
  long double mean_value, varians;
  int i,j;
  long double continuous_score_all, continuous_score_j, gamma_p_j;

  mean_value = (msa_seq + get_mtx_index(p,0,a_size+1))->share;
  
  
  continuous_score_all = 0.0;
  for(j = 0; j < a_size / 3; j++) {
    continuous_score_all += 
      get_single_gaussian_statescore(*(emissions + get_mtx_index(k, (j * 3), a_size)),
				     *(emissions + get_mtx_index(k, (j * 3 + 1), a_size)),
				     mean_value) *
      *((emissions) + (k * (a_size)) + (j * 3 + 2));
  }
  
  for(j = 0; j < a_size / 3; j++) {
    continuous_score_j = 
      get_single_gaussian_statescore(*(emissions + get_mtx_index(k, (j * 3), a_size)),
				     *(emissions + get_mtx_index(k, (j * 3 + 1), a_size)),
				     mean_value) *
      *((emissions) + (k * (a_size)) + (j * 3 + 2));
    varians = pow((msa_seq + get_mtx_index(p,0,a_size+1))->share - *(emissions + get_mtx_index(k, j * 3, a_size)), 2);
    if(continuous_score_all > 0.0) {
      gamma_p_j = Eka_base * continuous_score_j / continuous_score_all;
    }
    else {
      gamma_p_j = 0.0;
    }
    *(E + get_mtx_index(k, j * 3, a_size + 1)) += mean_value * gamma_p_j;
    *(E + get_mtx_index(k, j * 3 + 1, a_size + 1)) += varians * gamma_p_j;
    *(E + get_mtx_index(k, j * 3 + 2, a_size + 1)) += gamma_p_j;
  }
  
  *(E + get_mtx_index(k, j * 3, a_size + 1)) += Eka_base;
}



void add_to_E_dot_product(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
			  int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_dot_product_picasso(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
				  int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_picasso(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
		      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_picasso_sym(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
		      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_score(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
			      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;



  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_reversed_score(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
				       int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_dot_product_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
				 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_dot_product_picasso_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
					 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_picasso_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
				 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_picasso_sym_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
				 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_score_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
				     int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_reversed_score_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
					      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_subst_mtx_product(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
			  int k, int a_size, int normalize, long double *subst_mtx)
{
  int a_index, a_index2;
  long double prf_column_length;
  long double subst_mtx_row_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      for(a_index2 = 0; a_index2 < a_size; a_index2++) {
	*(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	  (long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	  *(subst_mtx + get_mtx_index(a_index, a_index2, a_size));
      }
    }
  }
  else {
    printf("Error: no normalizing in subst_mtx_product, yet...\n");
    exit(0);
  }

}
void add_to_E_subst_mtx_product_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
				       int k, int a_size, int normalize, long double *subst_mtx)
{
  int a_index, a_index2;
  long double prf_column_length;
  long double subst_mtx_row_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      for(a_index2 = 0; a_index2 < a_size; a_index2++) {
	*(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	  (long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	  *(subst_mtx + get_mtx_index(a_index, a_index2, a_size));
      }
    }
  }
  else {
    printf("Error: no normalizing in subst_mtx_product, yet...\n");
    exit(0);
  }

}

void add_to_E_subst_mtx_dot_product(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
			  int k, int a_size, int normalize, long double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length, subst_mtx_row_length;
  int query_index;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}

void add_to_E_subst_mtx_dot_product_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
					   int k, int a_size, int normalize, long double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length, subst_mtx_row_length;
  int query_index;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}

void add_to_E_subst_mtx_dot_product_prior(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
					  int k, int a_size, int normalize, long double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length, subst_mtx_row_length;
  int query_index;
  long double rli;
  
  rli = REST_LETTER_INDEX;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) * rli;
    }
    
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) * rli / 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}


void add_to_E_subst_mtx_dot_product_prior_nr_occ(long double *E, long double Eka_base, struct msa_letter_s *msa_seq, int p,
						 int k, int a_size, int normalize, long double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  long double prf_column_length, subst_mtx_row_length;
  int query_index;
  long double rli;
  
  rli = REST_LETTER_INDEX;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) * rli;
    }
    
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(long double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) * rli / 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}



/* General versions of the functions needed for keeping track of the labeleings in the one-best algorithm */

void update_labelings(struct one_best_s *cur_rowp, char *vertex_labels, 
		      int *sorted_v_list, int seq_len, int c, char *labels, int nr_of_labels, int nr_v)
{
  int v,w;
  int v_list_index;
  int cur_address;
  int first;
  char *tmp_labeling;
  char cur_label;
  int **same_labeling_lists;
  int *same_labeling_list_indices;
  int i;
  

#ifdef DEBUG_LABELING_UPDATE
  dump_v_list(sorted_v_list);
  printf("nr of labels = %d\n", nr_of_labels);
  printf("dump of nr three\n");
#endif

  v_list_index = 0;

  same_labeling_list_indices = (int*)(malloc_or_die(nr_of_labels * sizeof(int)));
  
  same_labeling_lists = (int**)(malloc_or_die(nr_of_labels * sizeof(int*)));


  for(i = 0; i < nr_of_labels; i++) {
    same_labeling_lists[i] = (int*)(malloc_or_die((nr_v * 2 + 1) * sizeof(int)));
  }
 

  for(v = 0; v < nr_v; v++) {
    (cur_rowp + v)->is_updated = NO;
  }
  
  for(i = 0; i < nr_of_labels; i++) {
    same_labeling_list_indices[i] = 0;
  }


  for(v = 1; v < nr_v-1; v++) {
    /* find all states with same labeling as this state up to previous row */
    
    if((cur_rowp + v)->is_updated == NO && (cur_rowp + v)->labeling != NULL) {
      cur_address = (int)((cur_rowp+v)->labeling);
#ifdef DEBUG_LABELING_UPDATE
      printf("searching vertex %d\n", v);
#endif
      for(i = 0; i < nr_of_labels; i++) {
	if(*(vertex_labels + v) == *(labels + i)) {
	  *(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = v;
	  same_labeling_list_indices[i] += 1;
	  break;
	}
      }
      
      (cur_rowp+v)->is_updated = YES;
      for(w = v+1; w < nr_v-1; w++) {
	if((int)((cur_rowp+w)->labeling) == cur_address && (cur_rowp + w)->is_updated == NO) {
#ifdef DEBUG_LABELING_UPDATE	  
	  printf("found same address, vertex nr = %d\n", w);
#endif
	  for(i = 0; i < nr_of_labels; i++) {
	    if(*(vertex_labels + w) == *(labels + i)) {
	      *(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = w;
	      same_labeling_list_indices[i] += 1;
	      break;
	    }
	  }

	  (cur_rowp+w)->is_updated = YES;
	}
      }

      for(i = 0; i < nr_of_labels; i++) {
	*(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = END;
	same_labeling_list_indices[i] += 1;
      }
    }
  }
  for(i = 0; i < nr_of_labels; i++) {
    *(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = TOT_END;
  }


#ifdef DEBUG_LABELING_UPDATE
  for(i = 0; i < nr_of_labels; i++) {
    printf("same_labeling_lists, label: %c\n", labels[i]);
    dump_label_tmp_list(*(same_labeling_lists + i));
  }
  //exit(0);
#endif
  
  for(i = 0; i < nr_of_labels; i++) {
    same_labeling_list_indices[i] = 0;
    while(*(*(same_labeling_lists + i) + same_labeling_list_indices[i]) != TOT_END) {
      first = YES;
      while(*(*(same_labeling_lists + i) + same_labeling_list_indices[i]) != END) {
	/* update sorted_v_list */
	*(sorted_v_list + v_list_index) = *(*(same_labeling_lists + i) + same_labeling_list_indices[i]);
	v_list_index++;

	/* update pointers and label paths */
	if(first == YES) {
	  tmp_labeling = (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling;
	  (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling =
	    (char*)malloc_or_die((c+1) * sizeof(char));
	  memcpy((cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling,
		 tmp_labeling, (c) * sizeof(char));
	  ((cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling)[c] = labels[i];
#ifdef DEBUG_LABELING_UPDATE	 
	  printf("added label;labels[%d] = %c\n", i, labels[i]);
#endif
	  first = NO;
	  tmp_labeling = (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling;
	}
	else {
	  (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling = tmp_labeling;
	}
#ifdef DEBUG_LABELING_UPDATE
      printf("label length c = %d\n", c);
      dump_labeling((cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling, c);
#endif


	same_labeling_list_indices[i] += 1;
      }
      if(first == NO) {
	*(sorted_v_list + v_list_index) = V_LIST_NEXT;
	v_list_index++;
      }
      same_labeling_list_indices[i] += 1;
    }
  }
  
  *(sorted_v_list + v_list_index) = V_LIST_END;
  
  for(v = 1; v < nr_v; v++) {
    (cur_rowp + v)->is_updated = NO;
  }
  

  /* garbage collection */
  for(i = 0; i < nr_of_labels; i++) {
    free(same_labeling_lists[i]);
  }
  free(same_labeling_list_indices);
  free(same_labeling_lists);
}

void deallocate_row_labelings(struct one_best_s *prev_rowp, int nr_v)
{
  int dealloc_index;
  int cur_address;
  int v,w;
  int *dealloc_list;
  
#ifdef DEBUG_DEALLOCATE_LABELINGS
  printf("starting dealloc\n");
  printf("nr_v = %d\n", nr_v);
#endif

  dealloc_list = (int*)(malloc_or_die((nr_v+1) * sizeof(int)));

  for(v = 0; v < nr_v; v++) {
    (prev_rowp + v)->is_updated = NO;
  }
  
  
  /* deallocate last row's labelings */
  dealloc_index = 0;
  for(v = 0; v < nr_v; v++) {
    /* find all states with same labeling as this state up to previous row */
    if((prev_rowp + v)->is_updated == NO && (prev_rowp + v)->labeling != NULL) {
      cur_address = (int)((prev_rowp+v)->labeling);
      dealloc_list[dealloc_index] = v;
      dealloc_index++;
      (prev_rowp+v)->is_updated = YES;
      for(w = v+1; w < nr_v; w++) {
	if((int)((prev_rowp+w)->labeling) == cur_address) {
#ifdef DEBUG_DEALLOCATE_LABELINGS
	  printf("found same address, vertices %d and %d: %x\n", v, w, (prev_rowp + w)->labeling);
#endif	
	  (prev_rowp+w)->is_updated = YES;
	  (prev_rowp+w)->labeling = NULL;
	}
      }
    }
  }
  dealloc_list[dealloc_index] = END;
  
  for(dealloc_index = 0; dealloc_list[dealloc_index] != END; dealloc_index++) {
#ifdef DEBUG_DEALLOCATE_LABELINGS    
    printf("dealloc_index = %d\n", dealloc_index);
    printf("freeing labeling of vertex %d\n", dealloc_list[dealloc_index]);
#endif
    free((prev_rowp + dealloc_list[dealloc_index])->labeling);
#ifdef DEBUG_DEALLOCATE_LABELINGS    
    printf("done\n");
#endif
  } 
  
  free(dealloc_list);
#ifdef DEBUG_DEALLOCATE_LABELINGS  
  printf("exiting dealloc\n");
#endif
}



/*** help function for user_defined_emission_score ***/
long double pos_spec_DG (char aa, int i, int helix_length) {
  const long double profiles[26][5] = { {  0.1267255,  0.0215152,  0        ,  0        ,  0          },   //A
				   {  1.5556351,  0.0133523,  0        ,  0        ,  0          },   //B
				   { -0.0765051,  0.0994228,  0        ,  0        ,  0          },   //C
				   {  1.7939795,  0.0172643,  0        ,  0        ,  0          },   //D
				   {  1.4193720,  0.0089351,  0        ,  0        ,  0          },   //E
				   { -0.2766953,  0.0010297,  0        ,  0        ,  0          },   //F
				   {  0.4813492,  0.0047210,  0        ,  0        ,  0          },   //G
				   {  1.1998590,  0.0080127,  0        ,  0        ,  0          },   //H
				   { -0.4597384,  0.0181495,  0        ,  0        ,  0          },   //I
				   {  0        ,  0        ,  0        ,  0        ,  0          },
				   {  1.8485768,  0.0218446,  0        ,  0        ,  0          },   //K
				   { -0.4282992,  0.0023804,  0        ,  0        ,  0          },   //L
				   { -0.0774786,  0.0984413,  0        ,  0        ,  0          },   //M
				   {  1.3266132,  0.0092375,  0        ,  0        ,  0          },   //N
				   {  0        ,  0        ,  0        ,  0        ,  0          },
				   {  1.0860888,  0.0100568,  0        ,  0        ,  0          },   //P
				   {  1.3336109,  0.0111996,  0        ,  0        ,  0          },   //Q
				   {  1.6492534,  0.0512044,  0        ,  0        ,  0          },   //R
				   {  0.7023921,  0.0077661,  0        ,  0        ,  0          },   //S
				   {  0.5266550,  0.0311973,  0        ,  0        ,  0          },   //T
				   {  0        ,  0        ,  0        ,  0        ,  0          },
				   { -0.2447218,  0.0979201,  0        ,  0        ,  0          },   //V
				   {  0.2909390,  0.0189282, -0.5479140,  0.0930222,  6.4736619  },   //W
				   {  0.6348884,  0.0180273,  0        ,  0        ,  0          },   //X
				   {  0.6275249,  0.0103896, -0.5744404,  0.0947821,  6.9164963  },   //Y
				   {  1.3761092,  0.0099898,  0        ,  0        ,  0          } }; //Z
  

  long double pos = 9 * (2 * (i-1) / (helix_length - 1.0) - 1);
  long double DG;
  int aa_nr = (int) (aa-65);
  if ((aa == 'W') || (aa == 'Y')) {
    DG = profiles[aa_nr][0] * exp(-1*profiles[aa_nr][1]*pow(pos,2)) 
         + profiles[aa_nr][2] * (    exp(-1*profiles[aa_nr][3]*pow((pos-profiles[aa_nr][4]),2)) 
				   + exp(-1*profiles[aa_nr][3]*pow((pos+profiles[aa_nr][4]),2)) );
  } else {
    DG = profiles[aa_nr][0] * exp(-1*profiles[aa_nr][1]*pow(pos,2));
    //printf("aa = %c, aa_nr = %d, profiles[aa_nr][0] = %Lf, profiles[aa_nr][1] = %Lf\n", aa, aa_nr, profiles[aa_nr][0],profiles[aa_nr][1]);
  }
  //printf("DG = %Lf\n", aa, DG, );
  return DG;
}
  



#if TARGEToctopus || TARGETspoctopus

#if TARGEToctopus
#define VERTEX_NR_39 39
#define VERTEX_NR_40 40
#define VERTEX_NR_74 7
#define VERTEX_NR_75 75
#define VERTEX_NR_90 90
#define VERTEX_NR_109 109
#define VERTEX_NR_109 109
#define VERTEX_NR_110 110
#define VERTEX_NR_125 125
#define VERTEX_NR_144 144
#elif TARGETspoctopus
#define VERTEX_NR_39 35
#define VERTEX_NR_40 36
#define VERTEX_NR_74 66
#define VERTEX_NR_75 67
#define VERTEX_NR_90 82
#define VERTEX_NR_109 97
#define VERTEX_NR_110 98
#define VERTEX_NR_125 113
#define VERTEX_NR_144 128
#endif 


long double get_user_defined_emission_score(char label, struct letter_s *seq, int seq_pos)
{
  const long double bgr_freq[] = { 0.0785,    //A
			      0.0474,    //B
			      0.0150,    //C
			      0.0534,    //D
			      0.0666,    //E
			      0.0395,    //F
			      0.0694,    //G
			      0.0229,    //H
			      0.0589,    //I
			      0     , 
			      0.0593,    //K
			      0.0966,    //L
			      0.0239,    //M
			      0.0413,    //N
			      0     , 
			      0.0484,    //P
			      0.0396,    //Q
			      0.0542,    //R
			      0.0687,    //S
			      0.0540,    //T
			      0.0239, 	 //U  
			      0.0672,    //V
			      0.0113,    //W
			      0.0500, 	 //X  
			      0.0301,    //Y
			      0.0531};   //Z
  
 
  const long double c0 =             0.2704503;
  const long double c1 =             9.81504713559865;
  const long double c2 =            -0.69588262140385;
  const long double c3 =             0.00944060947990;
  const long double pi_over_180 =    0.01745329251994;

  const long double RETURN_IO = 0.048;
  const long double RETURN_J = 0.05;
  const long double DG2PROB = 0.582;
  const long double DG_SHIFT = 0;

  long double DG, DG_sum=0, DG_sin_sum=0, DG_cos_sum=0, segment_DG, prob;

  int helix_length=22;
  //printf("%c\n", ((seq + seq_pos)->letter)[0]);
  //printf("%c\n", label);
  //printf("%d\n", (int) (((seq + seq_pos)->letter)[0]-65));
  //printf("%Lf\n", bgr_freq[(int) (((seq + seq_pos)->letter)[0]-65)]);
  
  if (label == 'O' || label == 'I') {
    return RETURN_IO;
    //return bgr_freq[(int) (((seq + seq_pos)->letter)[0]-65)];
  }
  else if (label == 'J') {
    if(((seq + seq_pos)->letter)[0] == 'K' || ((seq + seq_pos)->letter)[0] == 'R') {
      return RETURN_J;
    }
    else {
      return RETURN_IO;
    }
  }
  else if (label == 'M') {
    return 0.05;
  }
  else if (label == 'T') {
    return RETURN_IO;
  }
  else if (label >= 'a' && label <= 'n') {
    int i, j=1;
    if(label <= 'g') {
      helix_length = label - 80;
      //printf("label: %c,  helix_length = %d\n", label, helix_length);
    }
    else {
      helix_length = label - 87;
      //printf("label: %c,  helix_length = %d\n", label, helix_length);
    }

    DG_sum=0.0;
    DG_sin_sum=0.0;
    DG_cos_sum=0.0;
    if(seq_pos <= helix_length) {
      return 0.0;
    }
    for (i=seq_pos-helix_length+1;i<=seq_pos;i++) {
      DG = pos_spec_DG(((seq+i)->letter)[0],j,helix_length);
      DG_sum += DG;
      DG_sin_sum += (DG * sin(100*(j-1)*pi_over_180));
      DG_cos_sum += (DG * cos(100*(j-1)*pi_over_180));
      j++;
      //printf("DG = %Lf, DG_sin_sum = %Lf, DG_cos_sum = %Lf\n", DG, DG_sin_sum, DG_cos_sum);
    }
    segment_DG = DG_sum + c0 * sqrt ( pow(DG_sin_sum,2) + pow(DG_cos_sum,2) ) + c1 + c2 * helix_length + c3 * pow(helix_length,2);
    prob = 1 - 1 / ( 1 + exp( -1 * (segment_DG - DG_SHIFT) / DG2PROB ) );
    //printf("pos = %d: prob = %Lf, segment_DG = %Lf\n", seq_pos, prob, segment_DG);
    return prob;
    } 
  else {
    //printf("%c\n", ((seq + seq_pos)->letter)[0]);
    //printf("false: %c\n", label);
    //printf("%d\n", (int) (((seq + seq_pos)->letter)[0]-65));
    //printf("%Lf\n", bgr_freq[(int) (((seq + seq_pos)->letter)[0]-65)]);
    //printf("Error\n");
  }
  return 0.0;
}


long double get_user_defined_emission_score_NNZHMM_msa(int a_size, int vertex_nr, struct msa_letter_s *msa_seq, int seq_pos, int seq_length)
{
  long double I_frac = 0;
  long double O_frac = 0;
  long double IO_frac = 0.5;
  int i,j;
  const int M_DEPTH = 4;
  const int IO_CLOSENESS = 4;
  long double frac = (long double)(M_DEPTH + IO_CLOSENESS);
  int nr_counted = 0;

  
  if(vertex_nr == 5) {
    IO_frac = 0;
    for(i = seq_pos + M_DEPTH - 1; i >= seq_pos - IO_CLOSENESS; i--) {
      if(i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,0, a_size+1))->share;
      }
    }
    
  }
  else if(vertex_nr == VERTEX_NR_39) {
    IO_frac = 0;
    for(i = seq_pos - M_DEPTH + 1; i <= seq_pos + IO_CLOSENESS; i++) {
      if(i <= seq_length-1 && i > 0) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,1, a_size+1))->share;
      }
    }
   
  }
  else if(vertex_nr == VERTEX_NR_40) {
    IO_frac = 0;
    for(i = seq_pos + M_DEPTH - 1; i >= seq_pos - IO_CLOSENESS; i--) {
      if(i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,1, a_size+1))->share;
      }
    }
    
  }
  else if(vertex_nr == VERTEX_NR_74) {
    IO_frac = 0;
    for(i = seq_pos - M_DEPTH + 1; i <= seq_pos + IO_CLOSENESS; i++) {
      if(i <= seq_length-1 && i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,0, a_size+1))->share;
      }
    }
  }
  else if(vertex_nr == VERTEX_NR_75) {
    IO_frac = 0;
    for(i = seq_pos - 1; i >= seq_pos - 10; i--) {
      if(i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,0, a_size+1))->share;
      }
    }
    
  }
  else if(vertex_nr == VERTEX_NR_90) {
    IO_frac = 0;
    for(i = seq_pos - 2; i <= seq_pos + 2; i++) {
      if(i <= seq_length-1 && i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,1, a_size+1))->share;
      }
    }
   
  }
  else if(vertex_nr == VERTEX_NR_109) {
    IO_frac = 0;
    for(i = seq_pos + 1; i <= seq_pos + 10; i++) {
      if(i <= seq_length-1 && i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,0, a_size+1))->share;
      }
    }
   
  }
  else if(vertex_nr == VERTEX_NR_110) {
    IO_frac = 0;
    for(i = seq_pos - 1; i >= seq_pos - 10; i--) {
      if(i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,1, a_size+1))->share;
      }
    }
    
  }
  else if(vertex_nr == VERTEX_NR_125) {
    IO_frac = 0;
    for(i = seq_pos - 2; i <= seq_pos + 2; i++) {
      if(i <= seq_length-1 && i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,0, a_size+1))->share;
      }
    }
   
  }
  else if(vertex_nr == VERTEX_NR_144) {
    IO_frac = 0;
    for(i = seq_pos + 1; i <= seq_pos + 10; i++) {
      if(i <= seq_length-1 && i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,1, a_size+1))->share;
      }
    }
  }


  if(nr_counted > 0) {
    IO_frac = IO_frac / nr_counted;
  }
 
  return IO_frac;
}

long double get_user_defined_emission_score_msa(int a_size, char label, struct msa_letter_s *msa_seq, int seq_pos)
{
  const long double bgr_freq[] = { 0.0785,    //A
				   0.0474,    //B
				   0.0150,    //C
				   0.0534,    //D
				   0.0666,    //E
				   0.0395,    //F
				   0.0694,    //G
				   0.0229,    //H
				   0.0589,    //I
				   0     , 
				   0.0593,    //K
				   0.0966,    //L
				   0.0239,    //M
				   0.0413,    //N
				   0     , 
				   0.0484,    //P
				   0.0396,    //Q
				   0.0542,    //R
				   0.0687,    //S
				   0.0540,    //T
				   0.0239,    //U  
				   0.0672,    //V
				   0.0113,    //W
				   0.0500,    //X  
				   0.0301,    //Y
				   0.0531};   //Z
  
 
  const long double c0 =             0.2704503;
  const long double c1 =             9.81504713559865;
  const long double c2 =            -0.69588262140385;
  const long double c3 =             0.00944060947990;
  const long double pi_over_180 =    0.01745329251994;

  const long double RETURN_IO = 0.047;
  const long double RETURN_J = 0.05;
  const long double DG2PROB = 0.582;
  const long double DG_SHIFT = 0;
  const int MAX_MEM_KR = 1;

  const char aa_alphabet[] = "ACDEFGHIKLMNPQRSTVWY";

  long double aa_frac = 0.0;
  long double DG, DG_sum=0, DG_sin_sum=0, DG_cos_sum=0, segment_DG, prob;
  int a_index;
  int helix_length = 0;

  //printf("label %c\n", label);

  if (label == 'O' || label == 'I') {
    return RETURN_IO;
  }
  else if (label == 'J') {
    /* loop getting the amino acid fraction for each aa in this position */
    long double KR_frac = (msa_seq + get_mtx_index(seq_pos,8, a_size+1))->share + (msa_seq + get_mtx_index(seq_pos,14, a_size+1))->share;
    //  for(a_index = 0; a_index < a_size; a_index++) {
    //    aa_frac = (msa_seq + get_mtx_index(seq_pos,a_index, a_size+1))->share;
    //    printf("label: %c, a_index: %d, aa_frac: %.3Lf\n", label, a_index, aa_frac);
    //  }
    return (RETURN_IO + RETURN_J * KR_frac);
  }
  else if (label == 'M') {
    return 0.05;
  }
  else if (label >= 'a' && label <= 'n') {
    int i, j=1;
    if(label <= 'g') {
      helix_length = label - 80;
      
    }
    else {
      helix_length = label - 87;
      
    }

    DG_sum=0.0;
    DG_sin_sum=0.0;
    DG_cos_sum=0.0;
    if(seq_pos <= helix_length) {
      //printf("seqpos: %d\n", seq_pos);
      return 0.0;
    }
    for (i=seq_pos-helix_length+1;i<=seq_pos;i++) {
      /* loop getting the amino acid fraction for each aa in this position */
      DG=0.0;
      for(a_index = 0; a_index < a_size; a_index++) {
	aa_frac = (msa_seq + get_mtx_index(i,a_index, a_size+1))->share;
	DG += aa_frac * pos_spec_DG(aa_alphabet[a_index],j,helix_length);
	//printf("label: %c, a_index: %d, aa_frac: %.3Lf\n", label, a_index, aa_frac);
      }

      //DG = pos_spec_DG(((seq+i)->letter)[0],j,helix_length);
      DG_sum += DG;
      DG_sin_sum += (DG * sin(100*(j-1)*pi_over_180));
      DG_cos_sum += (DG * cos(100*(j-1)*pi_over_180));
      j++;
     
    }
    segment_DG = DG_sum + c0 * sqrt ( pow(DG_sin_sum,2) + pow(DG_cos_sum,2) ) + c1 + c2 * helix_length + c3 * pow(helix_length,2);
    prob = 1 - 1 / ( 1 + exp( -1 * (segment_DG - DG_SHIFT) / DG2PROB ) );
    return prob;
    } 
  else {
    printf("Error\n");
    exit(0);
  }
  return 0.0;
}

#endif /* TARGEToctopus || TARGETspoctopus */


#if TARGETscampi

long double get_user_defined_emission_score(char label, struct letter_s *seq, int seq_pos)
{
const long double bgr_freq[] = { 0.0785,    //A
                                 0.0474,    //B
			      0.0150,    //C
			      0.0534,    //D
			      0.0666,    //E
			      0.0395,    //F
			      0.0694,    //G
			      0.0229,    //H
			      0.0589,    //I
			      0     , 
			      0.0593,    //K
			      0.0966,    //L
			      0.0239,    //M
			      0.0413,    //N
			      0     , 
			      0.0484,    //P
			      0.0396,    //Q
			      0.0542,    //R
			      0.0687,    //S
			      0.0540,    //T
			      0.0239, 	 //U  
			      0.0672,    //V
			      0.0113,    //W
			      0.0500, 	 //X  
			      0.0301,    //Y
			      0.0531};   //Z
  
 
  const long double c0 =             0.2704503;
  const long double c1 =             9.81504713559865;
  const long double c2 =            -0.69588262140385;
  const long double c3 =             0.00944060947990;
  const long double pi_over_180 =    0.01745329251994;

  const long double RETURN_IO = 0.05;
  const long double RETURN_J = 0.07;
  const long double DG2PROB = 0.582;
  const long double DG_SHIFT = 0.8;

  long double DG, DG_sum=0, DG_sin_sum=0, DG_cos_sum=0, segment_DG, prob;

  int helix_length=21;
  //printf("%c\n", ((seq + seq_pos)->letter)[0]);
  //printf("%c\n", label);
  //printf("%d\n", (int) (((seq + seq_pos)->letter)[0]-65));
  //printf("%Lf\n", bgr_freq[(int) (((seq + seq_pos)->letter)[0]-65)]);
  
  //  printf("pos = %d, label = %c\n",seq_pos, label);
  if (label == 'O' || label == 'I') {
    return RETURN_IO;
    //return bgr_freq[(int) (((seq + seq_pos)->letter)[0]-65)];
  }
  else if (label == 'J') {
    if(((seq + seq_pos)->letter)[0] == 'K' || ((seq + seq_pos)->letter)[0] == 'R') {
      return RETURN_J;
    }
    else {
      return RETURN_IO;
    }
  }
  else if (label == 'M') {
    return 0.05;
  }
  else if (label == 'T') {
    return RETURN_IO;
  }
  else if (label >= 'a' && label <= 'n') {
    int i, j=1;
    if(label <= 'g') {
      helix_length = label - 80;
      //printf("label: %c,  helix_length = %d\n", label, helix_length);
    }
    else {
      helix_length = label - 87;
      //printf("label: %c,  helix_length = %d\n", label, helix_length);
    }

    DG_sum=0.0;
    DG_sin_sum=0.0;
    DG_cos_sum=0.0;
    if(seq_pos <= helix_length) {
      return 0.0;
    }
    for (i=seq_pos-helix_length+1;i<=seq_pos;i++) {
      DG = pos_spec_DG(((seq+i)->letter)[0],j,helix_length);
      DG_sum += DG;
      DG_sin_sum += (DG * sin(100*(j-1)*pi_over_180));
      DG_cos_sum += (DG * cos(100*(j-1)*pi_over_180));
      j++;
      //printf("DG = %Lf, DG_sin_sum = %Lf, DG_cos_sum = %Lf\n", DG, DG_sin_sum, DG_cos_sum);
    }
    segment_DG = DG_sum + c0 * sqrt ( pow(DG_sin_sum,2) + pow(DG_cos_sum,2) ) + c1 + c2 * helix_length + c3 * pow(helix_length,2);
    prob = 1 - 1 / ( 1 + exp( -1 * (segment_DG - DG_SHIFT) / DG2PROB ) );
    //printf("pos = %d: prob = %Lf, segment_DG = %Lf\n", seq_pos, prob, segment_DG);
    return prob;
    } 
  else {
    //printf("%c\n", ((seq + seq_pos)->letter)[0]);
    //printf("false: %c\n", label);
    //printf("%d\n", (int) (((seq + seq_pos)->letter)[0]-65));
    //printf("%Lf\n", bgr_freq[(int) (((seq + seq_pos)->letter)[0]-65)]);
    //printf("Error\n");
  }
  return 0.0;
}


long double get_user_defined_emission_score_NNZHMM_msa(int a_size, int vertex_nr, struct msa_letter_s *msa_seq, int seq_pos, int seq_length)
{
  long double I_frac = 0;
  long double O_frac = 0;
  long double IO_frac = 0.5;
  int i,j;
  const int M_DEPTH = 4;
  const int IO_CLOSENESS = 4;
  long double frac = (long double)(M_DEPTH + IO_CLOSENESS);
  int nr_counted = 0;

  if(vertex_nr == 5) {
    IO_frac = 0;
    for(i = seq_pos + M_DEPTH - 1; i >= seq_pos - IO_CLOSENESS; i--) {
      if(i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,0, a_size+1))->share;
      }
    }
    
  }
  else if(vertex_nr == 27) {
    IO_frac = 0;
    for(i = seq_pos - M_DEPTH + 1; i <= seq_pos + IO_CLOSENESS; i++) {
      if(i <= seq_length-1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,1, a_size+1))->share;
      }
    }
   
  }
  else if(vertex_nr == 28) {
    IO_frac = 0;
    for(i = seq_pos + M_DEPTH - 1; i >= seq_pos - IO_CLOSENESS; i--) {
      if(i >= 1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,1, a_size+1))->share;
      }
    }
    
  }
  else if(vertex_nr == 50) {
    IO_frac = 0;
    for(i = seq_pos - M_DEPTH + 1; i <= seq_pos + IO_CLOSENESS; i++) {
      if(i <= seq_length-1) {
	nr_counted++;
	IO_frac += (msa_seq + get_mtx_index(i,0, a_size+1))->share;
      }
    }
    
  }
  if(nr_counted > 0) {
    IO_frac = IO_frac / nr_counted;
  }
  //IO_frac = IO_frac / frac;
  //IO_frac = IO_frac;

  if(vertex_nr == 27) {
    //printf("%Lf,   %d\n", IO_frac, seq_pos);
  }
  //return 1.0;
  //IO_frac = 0.5;
  return IO_frac;
}

long double get_user_defined_emission_score_msa(int a_size, char label, struct msa_letter_s *msa_seq, int seq_pos)
{
  const long double bgr_freq[] = { 0.0785,    //A
				   0.0474,    //B
				   0.0150,    //C
				   0.0534,    //D
				   0.0666,    //E
				   0.0395,    //F
				   0.0694,    //G
				   0.0229,    //H
				   0.0589,    //I
				   0     , 
				   0.0593,    //K
				   0.0966,    //L
				   0.0239,    //M
				   0.0413,    //N
				   0     , 
				   0.0484,    //P
				   0.0396,    //Q
				   0.0542,    //R
				   0.0687,    //S
				   0.0540,    //T
				   0.0239,    //U  
				   0.0672,    //V
				   0.0113,    //W
				   0.0500,    //X  
				   0.0301,    //Y
				   0.0531};   //Z
  
 
  const long double c0 =             0.2704503;
  const long double c1 =             9.81504713559865;
  const long double c2 =            -0.69588262140385;
  const long double c3 =             0.00944060947990;
  const long double pi_over_180 =    0.01745329251994;

  const long double RETURN_IO = 0.05;
  const long double RETURN_J = 0.04;
  const long double DG2PROB = 0.582;
  const long double DG_SHIFT = 1.1;
  const int MAX_MEM_KR = 1;

  const char aa_alphabet[] = "ACDEFGHIKLMNPQRSTVWY";

  long double aa_frac = 0.0;
  long double DG, DG_sum=0, DG_sin_sum=0, DG_cos_sum=0, segment_DG, prob;
  int a_index;
  int helix_length = 0;

  //printf("label %c\n", label);

  if (label == 'O' || label == 'I') {
    return RETURN_IO;
  }
  else if (label == 'J') {
    /* loop getting the amino acid fraction for each aa in this position */
    long double KR_frac = (msa_seq + get_mtx_index(seq_pos,8, a_size+1))->share + (msa_seq + get_mtx_index(seq_pos,14, a_size+1))->share;
    //  for(a_index = 0; a_index < a_size; a_index++) {
    //    aa_frac = (msa_seq + get_mtx_index(seq_pos,a_index, a_size+1))->share;
    //    printf("label: %c, a_index: %d, aa_frac: %.3Lf\n", label, a_index, aa_frac);
    //  }
    return (RETURN_IO + RETURN_J * KR_frac);
  }
  else if (label == 'M') {
    return 0.05;
  }
  else if (label >= 'a' && label <= 'n') {
    int i, j=1;
    if(label <= 'g') {
      helix_length = label - 80;
      
    }
    else {
      helix_length = label - 87;
      
    }

    DG_sum=0.0;
    DG_sin_sum=0.0;
    DG_cos_sum=0.0;

    if(seq_pos < helix_length) {
      //printf("seqpos: %d\n", seq_pos);
      return 0.0;
    }
    //    printf("seqpos: %d, helix length: %d\n", seq_pos,helix_length);
    for (i=seq_pos-helix_length+1;i<=seq_pos;i++) {
      /* loop getting the amino acid fraction for each aa in this position */
      DG=0.0;
      for(a_index = 0; a_index < a_size; a_index++) {
	aa_frac = (msa_seq + get_mtx_index(i,a_index, a_size+1))->share;
	DG += aa_frac * pos_spec_DG(aa_alphabet[a_index],j,helix_length);
	//printf("label: %c, a_index: %d, aa_frac: %.3Lf\n", label, a_index, aa_frac);
      }

      //DG = pos_spec_DG(((seq+i)->letter)[0],j,helix_length);
      DG_sum += DG;
      DG_sin_sum += (DG * sin(100*(j-1)*pi_over_180));
      DG_cos_sum += (DG * cos(100*(j-1)*pi_over_180));
      j++;
     
    }
    segment_DG = DG_sum + c0 * sqrt ( pow(DG_sin_sum,2) + pow(DG_cos_sum,2) ) + c1 + c2 * helix_length + c3 * pow(helix_length,2);
    prob = 1 - 1 / ( 1 + exp( -1 * (segment_DG - DG_SHIFT) / DG2PROB ) );
    //    printf("prob: %.3Lf\n",prob);
    return prob;
    } 
  else {
    printf("Error\n");
    exit(0);
  }
  return 0.0;
}

  
#endif /* TARGETscampi */
