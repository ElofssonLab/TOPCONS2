#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "funcs.h"
#include "cmdline_add_alphabet.h"

#define MAX_LINE 4000
#define MAX_SEQS 4000

extern int verbose;

void read_additional_alphabet(FILE *alphafile, struct hmm_multi_s *hmmp, int alphabet_nr);
void printhelp();

/* memory for transition and emission matrices, from_vertex array and 
   to_vertex_array will be allocated in readhmm, but must be freed here */
static struct hmm_multi_s hmm;


int main(int argc, char* argv[])
{
  int i;
  FILE *hmmfile, *outfile, *alphafile, *outfile_tmp;
  int nr_alphabets;
  int alphabet_nr;
  int hmmfiletype;

  struct gengetopt_args_info args_info;

  hmmfile = NULL;
  outfile = NULL;
  alphafile = NULL;
  outfile_tmp = NULL;
  verbose = NO;
  
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
  if(args_info.outfile_given) {
    if((outfile = fopen(args_info.outfile_arg, "w")) == NULL) {
      perror(args_info.outfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for writing\n",args_info.outfile_arg);
    }
  }
  

  
  if(args_info.alphafile_given) {
    if((alphafile = fopen(args_info.alphafile_arg, "r")) == NULL) {
      perror(args_info.alphafile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading alphabet file\n",args_info.alphafile_arg);
    }
  }
  
  if(args_info.verbose_given) {
    verbose = YES;
  }
  
  
  /* find out nr of alphabets */
  nr_alphabets = get_nr_alphabets(hmmfile);
  rewind(hmmfile);
  

  
  /* get hmm from file */
  printf("starting reading\n");  
  if(hmmfile != NULL) {
    hmmfiletype = readhmm_check(hmmfile);
    if(hmmfiletype == SINGLE_HMM) {
      readhmm(hmmfile, &hmm);
    }
    else if(hmmfiletype == MULTI_HMM) {
      readhmm_multialpha(hmmfile, &hmm);
    }
  }
  else {
    /* cannot happen */
    printf("Internal error: hmmfile not found, should have been detected earlier\n");
  }



  /* read additional alphabet info */
  if(hmm.nr_alphabets > 3) {
    printf("Error: maximum number of alphabets is 4\n");
    exit(0);
  }
  else {
    alphabet_nr = hmm.nr_alphabets + 1;
    read_additional_alphabet(alphafile, &hmm, alphabet_nr);
  }
  


  /* save adjusted hmm */
  printf("finished reading, starting save\n");
  savehmm_multialpha(outfile, &hmm);
  printf("finished save\n");
  

  /* garbage collection */
  hmm_garbage_collection_multi(hmmfile, &hmm);
  fclose(outfile);
  fclose(alphafile);
 
  /*****************************************************************************************/
  /* remember that garbage collection has not been done here */

}


void read_additional_alphabet(FILE *alphafile, struct hmm_multi_s *hmmp, int alphabet_nr)
{
  char s[4000], *prob_name;
  long double *emissions, *log_emissions, *emiss_prior_scalers;
  long double prob;
  struct module_multi_s *module;
  char null[5] = "null";

  int a_size;
  int i, j, k;
  

  while(fgets(s, MAX_LINE, alphafile) != NULL) {
    if(strncmp(s, "ALPHABET:", 9) == 0) {
      break;
    }
    
  }
  

  /* read alphabet */
  switch(alphabet_nr) {
  case 1:
    strcpy(hmmp->alphabet, &s[10]);
    break;
  case 2:
    strcpy(hmmp->alphabet_2, &s[10]);
    break;
  case 3:
    strcpy(hmmp->alphabet_3, &s[10]);
    break;
  case 4:
    strcpy(hmmp->alphabet_4, &s[10]);
    break;  
  }
  
  while(fgets(s, MAX_LINE, alphafile) != NULL) {
    if(strncmp(s, "ALPHABET LENGTH:", 16) == 0) {
      break;
    }
  }
  
  /* read a_size */
  a_size = atoi(&s[17]);
  switch(alphabet_nr) {
  case 1:
    hmmp->a_size = atoi(&s[17]);
    break;
  case 2:
    hmmp->a_size_2 = atoi(&s[17]);
    break;
  case 3:
    hmmp->a_size_3 = atoi(&s[17]);
    break;
  case 4:
    hmmp->a_size_4 = atoi(&s[17]);
    break;  
  }
  
  
  emissions = (long double*)(malloc_or_die(hmmp->nr_v * a_size * 
				      sizeof(long double)));
  log_emissions = (long double*)(malloc_or_die(hmmp->nr_v * a_size * 
				      sizeof(long double)));
  emiss_prior_scalers = (long double*)(malloc_or_die(hmmp->nr_v * 
				      sizeof(long double)));
  
  /*read emission probs */
  i = 0;
  while(fgets(s, MAX_LINE, alphafile) != NULL) {
    if(strncmp(s, "VERTEX ", 7) == 0) {
      if(i >= hmmp->nr_v) {
	printf("Warning: Alphabet file has more vertex rows than HMM, reading stopped at vertex %d\n", i);
	break;
      }

      strtok(s," ");
      strtok(NULL, " ");
      

      for(j = 0; j < a_size; j++) {
	/* open the priorfile */
	if((prob_name = strtok(NULL, " ")) == NULL) {
	  printf("Error: Too few probabilities on vertex prob row.\n");
	  exit(0);
	}
	else {
	  prob = atof(prob_name);
	  *(emissions + get_mtx_index(i,j,a_size)) = prob;
	} 
      }     
      
      
      i++;
    }
    else {
      
    }
  }
  
  for(i = 0; i < hmmp->nr_v; i++) {
    *(emiss_prior_scalers + i) = 1.0; 
    for(j = 0; j < a_size; j++) {
      if(*(emissions + get_mtx_index(i,j,a_size)) == 0.0) {
	*(log_emissions + get_mtx_index(i,j,a_size)) = DEFAULT;
      }
      else {
	*(log_emissions + get_mtx_index(i,j,a_size)) = log10(*(emissions + get_mtx_index(i,j,a_size)));
      }
    }
  }

  for(k = 0; k < hmmp->nr_m; k++) {
    module = *(hmmp->modules + k);
    switch(alphabet_nr) {
    case 1:
      memcpy(module->priorfile_name, null, sizeof(char) * 5);
      break;
    case 2:
      memcpy(module->priorfile_name_2, null, sizeof(char) * 5);
      break;
    case 3:
      memcpy(module->priorfile_name_3, null, sizeof(char) * 5);
      break;
    case 4:
      memcpy(module->priorfile_name_4, null, sizeof(char) * 5);
      break;  
    }
  }

  
  /* struct info (some arbitrary) */
  switch(alphabet_nr) {
  case 1:
    hmmp->emissions = emissions;
    hmmp->log_emissions = log_emissions;
    hmmp->vertex_emiss_prior_scalers = emiss_prior_scalers;
    hmmp->nr_ed = 0;
    hmmp->emission_dirichlets = NULL;
    hmmp->ed_ps = NULL;
    hmmp->subst_mtx = NULL;
    break;
  case 2:
    hmmp->emissions_2 = emissions;
    hmmp->log_emissions_2 = log_emissions;
    hmmp->vertex_emiss_prior_scalers_2 = emiss_prior_scalers;
    hmmp->nr_ed_2 = 0;
    hmmp->emission_dirichlets_2 = NULL;
    hmmp->ed_ps_2 = NULL;
    hmmp->subst_mtx_2 = NULL;
    break;
  case 3:
    hmmp->emissions_3 = emissions;
    hmmp->log_emissions_3 = log_emissions;
    hmmp->vertex_emiss_prior_scalers_3 = emiss_prior_scalers;
    hmmp->nr_ed_3 = 0;
    hmmp->emission_dirichlets_3 = NULL;
    hmmp->ed_ps_3 = NULL;
    hmmp->subst_mtx_3 = NULL;
    break;
  case 4:
    hmmp->emissions_4 = emissions;
    hmmp->log_emissions_4 = log_emissions;
    hmmp->vertex_emiss_prior_scalers_4 = emiss_prior_scalers;
    hmmp->nr_ed_4 = 0;
    hmmp->emission_dirichlets_4 = NULL;
    hmmp->ed_ps_4 = NULL;
    hmmp->subst_mtx_4 = NULL;
    break;
  }

  hmmp->nr_alphabets = alphabet_nr;

}
