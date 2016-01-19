#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "structs.h" /* data structures etc */
#include "funcs.h" /* function header */

#define POS 0

//#define SAVE_DEBUG

extern int verbose;

void savemodule_multi(FILE*, struct hmm_multi_s*, struct module_multi_s*);


savehmm_multialpha(FILE *outfile, struct hmm_multi_s *hmmp)
{
  char s[1000];
  char nr[32];
  time_t t;
  int t1, t2;
  int i,j;
  struct module_multi_s *modulep;
  struct transition_s trans, *trans_ties;
  
  if(outfile == NULL) {
    printf("Could not write to outfile when saving hmm\n");
    exit(0);
  }

  if(verbose == YES) {
    printf("saving hmm %s ... ", hmmp->name);
    fflush(stdout);
  }

  /* header */
  if(fputs("***********************Header***************************\n", outfile) ==
     EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("********Header********\n");
#endif
  /* name */
  strcpy(s, "NAME: ");
  if(fputs(strcat(s , hmmp->name), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* time */
  time(&t);
  strcpy(s, "TIME OF CREATION: ");
  if(fputs(strcat(s, ctime(&t)), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif

  /* nr of alphabets */
  strcpy(s, "NR OF ALPHABETS: ");
  itoa(nr, hmmp->nr_alphabets); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  
  /* alphabet */
  strcpy(s, "ALPHABET 1: ");
  if(fputs(strcat(s, hmmp->alphabet), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* alphabet length */
  strcpy(s, "ALPHABET LENGTH 1: ");
  itoa(nr, hmmp->a_size); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif

  if(hmmp->nr_alphabets > 1) {
    /* alphabet */
    strcpy(s, "ALPHABET 2: ");
    if(fputs(strcat(s, hmmp->alphabet_2), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    /* alphabet length */
    strcpy(s, "ALPHABET LENGTH 2: ");
    itoa(nr, hmmp->a_size_2); 
    if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }
  
  if(hmmp->nr_alphabets > 2) {
    /* alphabet */
    strcpy(s, "ALPHABET 3: ");
    if(fputs(strcat(s, hmmp->alphabet_3), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    /* alphabet length */
    strcpy(s, "ALPHABET LENGTH 3: ");
    itoa(nr, hmmp->a_size_3); 
    if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }

  if(hmmp->nr_alphabets > 3) {
    /* alphabet */
    strcpy(s, "ALPHABET 4: ");
    if(fputs(strcat(s, hmmp->alphabet_4), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    /* alphabet length */
    strcpy(s, "ALPHABET LENGTH 4: ");
    itoa(nr, hmmp->a_size_4); 
    if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }
  
  /* nr of modules */
  strcpy(s, "NR OF MODULES: ");
  itoa(nr, hmmp->nr_m); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* nr of vertices */
  strcpy(s, "NR OF VERTICES: ");
  itoa(nr, hmmp->nr_v); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* nr of TRANSITIONS */
  strcpy(s, "NR OF TRANSITIONS: ");
  itoa(nr, hmmp->nr_t); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* nr of DISTRIBUTION GROUPS */
  strcpy(s, "NR OF DISTRIBUTION GROUPS: ");
  itoa(nr, hmmp->nr_d); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* nr of TRANS TIE GROUPS */
  strcpy(s, "NR OF TRANSITION TIE GROUPS: ");
  itoa(nr, hmmp->nr_ttg); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* nr of PRIOR FILES */
  strcpy(s, "NR OF EMISSION PRIORFILES 1: ");
  itoa(nr, hmmp->nr_ed); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* the PRIOR FILES */
  strcpy(s, "EMISSION PRIORFILES_1:");
  for(i = 0; i < hmmp->nr_ed; i++) {
    strcat(s, " ");
    strcat(s, (hmmp->emission_dirichlets + i)->name);
  }
  if(fputs(strcat(s, " \n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif

  if(hmmp->nr_alphabets > 1) {
    /* nr of PRIOR FILES */
    strcpy(s, "NR OF EMISSION PRIORFILES 2: ");
    itoa(nr, hmmp->nr_ed_2); 
    if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    /* the PRIOR FILES */
    strcpy(s, "EMISSION PRIORFILES_2:");
    for(i = 0; i < hmmp->nr_ed_2; i++) {
      strcat(s, " ");
      strcat(s, (hmmp->emission_dirichlets_2 + i)->name);
    }
    if(fputs(strcat(s, " \n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }
  
  if(hmmp->nr_alphabets > 2) {
    /* nr of PRIOR FILES */
    strcpy(s, "NR OF EMISSION PRIORFILES 3: ");
    itoa(nr, hmmp->nr_ed_3); 
    if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    /* the PRIOR FILES */
    strcpy(s, "EMISSION PRIORFILES_3:");
    for(i = 0; i < hmmp->nr_ed_3; i++) {
      strcat(s, " ");
      strcat(s, (hmmp->emission_dirichlets_3 + i)->name);
    }
    if(fputs(strcat(s, " \n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }

  if(hmmp->nr_alphabets > 3) {
    /* nr of PRIOR FILES */
    strcpy(s, "NR OF EMISSION PRIORFILES 4: ");
    itoa(nr, hmmp->nr_ed_4); 
    if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    /* the PRIOR FILES */
    strcpy(s, "EMISSION PRIORFILES_4:");
    for(i = 0; i < hmmp->nr_ed_4; i++) {
      strcat(s, " ");
      strcat(s, (hmmp->emission_dirichlets_4 + i)->name);
    }
    if(fputs(strcat(s, " \n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }

  /* nr of PRIOR FILES */
  strcpy(s, "NR OF TRANSITION PRIORFILES: "); 
  if(fputs(strcat(s, "0\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* the PRIOR FILES */
  strcpy(s, "TRANSITION PRIORFILES:");
  for(i = 0; i < 0; i++) {
    strcat(s, " ");
    //strcat(s, (hmmp->emission_dirichlets)->name);
  }
  if(fputs(strcat(s, "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* two empty rows */
  if(fputs("\n\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif


  /* modules */
  if(fputs("***********************Modules*************************\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* save the modules */
  for(i = 0; i < hmmp->nr_m; i++) {
    //printf("inne\n");
    modulep = *(hmmp->modules + i);
    savemodule_multi(outfile, hmmp, modulep);
  }
  
  /* empty row */
  if(fputs("\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  
  /* save distribution groups */
  /****************Distribution groups***************************/
  if(fputs("**********************Emission distribution groups*******************\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("************************Emission distribution groups******************\n", s);
#endif

  /* the groups */
  j = 0;
  for(i = 0; i < hmmp->nr_d; i++) {
    strcpy(s, "Group ");
    itoa(nr, i+1);
    strcat(s, nr);
    strcat(s, ": ");
    while((t1 = *(hmmp->distrib_groups+j)) != END) {
      itoa(nr, t1);
      strcat(s, nr);
      strcat(s, " ");
      j++;
    }
    j++;
    strcat(s, "\n");
    if(fputs(s, outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }
  
  /* empty row */
  if(fputs("\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  
  /* save transition tie groups */
  /****************Transition tie groups***************************/
  if(fputs("**********************Transition tie groups*******************\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
  printf("hmmp->nr_ttg = %d\n", hmmp->nr_ttg);
#endif

  /* the groups */
  j = 0;
  for(i = 0; i < hmmp->nr_ttg; i++) {
    strcpy(s, "Tie ");
    itoa(nr, i+1);
    strcat(s, nr);
    strcat(s, ": ");
    while(((hmmp->trans_tie_groups)+j)->from_v != END) {
      t1 = ((hmmp->trans_tie_groups)+j)->from_v;
      itoa(nr, t1);
      strcat(s, nr);
      strcat(s, "->");
      t1 = ((hmmp->trans_tie_groups)+j)->to_v;
      itoa(nr, t1);
      strcat(s, nr);
      strcat(s, " ");
      j++;
    }
    j++;
    strcat(s, "\n");
    if(fputs(s, outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
  }

  /* flush buffer */
  fflush(outfile);
  
  if(verbose == YES) {
    printf("done\n");
  }
}

void savemodule_multi(FILE *outfile,struct hmm_multi_s *hmmp, struct module_multi_s *modulep)
{
  char s[500];
  char nr[30];
  int i, j, k ,l, m;
  int nr_t, nr_et, nr_e;
  int cur_v;
  char *alphp;
  
  //printf("name %s\n", modulep->name);

  /* name */
  strcpy(s, "Module: ");
  if(fputs(strcat(s , modulep->name), outfile) == EOF) {
    perror("");
  }
  
#ifdef SAVE_DEBUG
  printf("\n%s", s);
#endif
  /* type */
  strcpy(s, "Type: ");
  switch(modulep->type) {
  case SINGLENODE: strcat(s, "Singlenode\n"); break;
  case CLUSTER: strcat(s, "Cluster\n"); break;
  case HIGHWAY: strcat(s, "Highway\n"); break;
  case U_TURN: strcat(s, "U_Turn\n"); break;
  case SINGLELOOP: strcat(s, "Singleloop\n"); break;
  case FORWARD_STD: strcat(s, "Forward_std\n"); break;
  case FORWARD_ALT: strcat(s, "Forward_alt\n"); break; 
  case PROFILE7:  strcat(s, "Profile7\n"); break;
  case PROFILE9:  strcat(s, "Profile9\n"); break;
  default: printf("Error: unknown type %d, can't save hmm\n"); exit(-1);
  }
  if(fputs(s, outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* nr of vertices */
  strcpy(s, "NrVertices: ");
  itoa(nr, modulep->nr_v); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* emission prior file */
  strcpy(s, "Emission prior file 1: ");
  if(fputs(strcat(strcat(s, modulep->priorfile_name), "\0"), outfile) == EOF) {
    perror("");
  }
  if(fputs("\n", outfile) == EOF) {
    perror("");
  }

  if(hmmp->nr_alphabets > 1) {
    /* emission prior file */
    strcpy(s, "Emission prior file 2: ");
    if(fputs(strcat(strcat(s, modulep->priorfile_name_2), "\0"), outfile) == EOF) {
      perror("");
    }
    if(fputs("\n", outfile) == EOF) {
      perror("");
    }
  }

  if(hmmp->nr_alphabets > 2) {
    /* emission prior file */
    strcpy(s, "Emission prior file 3: ");
    if(fputs(strcat(strcat(s, modulep->priorfile_name_3), "\0"), outfile) == EOF) {
      perror("");
    }
    if(fputs("\n", outfile) == EOF) {
      perror("");
    }
  }

  if(hmmp->nr_alphabets > 3) {
    /* emission prior file */
    strcpy(s, "Emission prior file 4: ");
    if(fputs(strcat(strcat(s, modulep->priorfile_name_4), "\0"), outfile) == EOF) {
      perror("");
    }
    if(fputs("\n", outfile) == EOF) {
      perror("");
    }
  }
  
  /* transition prior file */
  strcpy(s, "Transition prior file: null\0");
  if(fputs(s, outfile) == EOF) {
    perror("");
  }
  /* empty row */
  if(fputs("\n\n", outfile) == EOF) {
    perror("");
  }
  for(i = 0; i < modulep->nr_v; i++) {
    
#ifdef SAVE_DEBUG
    printf("saving vertex %d...\n", *(modulep->vertices + i));
    printf("nr_v left: %d\n", modulep->nr_v - i);
#endif
    /* Vertex nr */
    cur_v = *(modulep->vertices + i);
    strcpy(s, "Vertex ");
    itoa(nr, cur_v); 
    if(fputs(strcat(strcat(s, nr), ":\n"), outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    /* vertex type */
    strcpy(s, "Vertex type: ");
    switch(modulep->v_type) {
    case STARTV: strcat(s, "start\n"); break;
    case ENDV: strcat(s, "end\n"); break;
    case STANDARDV: strcat(s, "standard\n"); break;
    case LOCKEDV: strcat(s, "locked\n"); break;
    case SILENTV: strcat(s, "silent\n"); break;
    case PROFILEV: strcat(s, get_profile_vertex_type(cur_v, hmmp->silent_vertices)); break;
    default: printf("Error: unknown vertex type, can't save hmm\n"); exit(-1);
    }
    if(fputs(s, outfile) == EOF) {
      perror("");
    }
#ifdef SAVE_DEBUG
    printf("%s", s);
#endif
    
    /* vertex label */
    strcpy(s, "Vertex label:  \n");
    s[14] = *(hmmp->vertex_labels + cur_v);
    if(fputs(s, outfile) == EOF) {
      perror("");
    }

    /* Transition prior scaler */
    strcpy(s, "Transition prior scaler: ");
    ftoa(nr, *(hmmp->vertex_trans_prior_scalers + cur_v), 6);
    strcat(s, nr);
    strcat(s, "\n");
    if(fputs(s, outfile) == EOF) {
      perror("");
    }
    
    /* Emission prior scaler */
    strcpy(s, "Emission prior scaler 1: ");
    ftoa(nr, *(hmmp->vertex_emiss_prior_scalers + cur_v), 6);
    strcat(s, nr);
    strcat(s, "\n");
    if(fputs(s, outfile) == EOF) {
      perror("");
    }
    
    if(hmmp->nr_alphabets > 1) {
      strcpy(s, "Emission prior scaler 2: ");
      ftoa(nr, *(hmmp->vertex_emiss_prior_scalers_2 + cur_v), 6);
      strcat(s, nr);
      strcat(s, "\n");
      if(fputs(s, outfile) == EOF) {
	perror("");
      }
    }

    if(hmmp->nr_alphabets > 2) {
      strcpy(s, "Emission prior scaler 3: ");
      ftoa(nr, *(hmmp->vertex_emiss_prior_scalers_3 + cur_v), 6);
      strcat(s, nr);
      strcat(s, "\n");
      if(fputs(s, outfile) == EOF) {
	perror("");
      }
    }

    if(hmmp->nr_alphabets > 3) {
      strcpy(s, "Emission prior scaler 4: ");
      ftoa(nr, *(hmmp->vertex_emiss_prior_scalers_4 + cur_v), 6);
      strcat(s, nr);
      strcat(s, "\n");
      if(fputs(s, outfile) == EOF) {
	perror("");
      }
    }


    /* nr transitions */
     strcpy(s, "Nr transitions = ");
     nr_t = 0;
     for(j = 0; j < hmmp->nr_v-1; j++) {
       if(*(hmmp->transitions + get_mtx_index(cur_v, j, hmmp->nr_v)) != 0) {
	 nr_t++;
       }
     }
     itoa(nr, nr_t);
     if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("%s", s);
#endif
     /* nr end transitions */
     strcpy(s, "Nr end transitions = ");
     nr_et = 0;
     if(*(hmmp->transitions + get_mtx_index(cur_v, hmmp->nr_v-1, hmmp->nr_v)) != 0.0) {
       nr_et++;
     }
     itoa(nr, nr_et);
     if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("%s", s);
#endif
     /* nr emissions */
     strcpy(s, "Nr emissions 1 = ");
     itoa(nr, hmmp->a_size);
     if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("%s", s);
#endif

     if(hmmp->nr_alphabets > 1) {
       /* nr emissions */
       strcpy(s, "Nr emissions 2 = ");
       itoa(nr, hmmp->a_size_2);
       if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
	 perror("");
       }
#ifdef SAVE_DEBUG
       printf("%s", s);
#endif
     }

     if(hmmp->nr_alphabets > 2) {
       /* nr emissions */
       strcpy(s, "Nr emissions 3 = ");
       itoa(nr, hmmp->a_size_3);
       if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
	 perror("");
       }
#ifdef SAVE_DEBUG
       printf("%s", s);
#endif
     }

     if(hmmp->nr_alphabets > 3) {
       /* nr emissions */
       strcpy(s, "Nr emissions 4 = ");
       itoa(nr, hmmp->a_size_4);
       if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
	 perror("");
       }
#ifdef SAVE_DEBUG
       printf("%s", s);
#endif
     }
     

     /* transition probabilities */
     if(fputs("Transition probabilities\n", outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("Transition probabilities\n");
     //dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
     printf("saving transp probs:\n");
#endif
    
     for(k = 0; k < hmmp->nr_v-1; k++) {
       if(*(hmmp->transitions + get_mtx_index(cur_v, k, hmmp->nr_v)) != 0.0) {
	 strcpy(s, "\tVertex ");
	 itoa(nr, k);
	 strcat(s, nr);
	 strcat(s, ": ");
	 ftoa(nr, *(hmmp->transitions + get_mtx_index(cur_v, k, hmmp->nr_v)), 6);
	 strcat(s, nr);
	 if(fputs(strcat(s, "\n"), outfile) == EOF) {
	   perror("");
	 }
#ifdef SAVE_DEBUG
	 printf("%s", s);
#endif
       }
     }
    
     
     /* end transition probabilities */
     if(fputs("End transition probabilities\n", outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("saving end trans probs\n");
#endif
     if(*(hmmp->transitions + get_mtx_index(cur_v, hmmp->nr_v-1, hmmp->nr_v)) != 0.0) {
       strcpy(s, "\tVertex ");
       itoa(nr, hmmp->nr_v - 1);
       strcat(s, nr);
       strcat(s, ": ");
       ftoa(nr, *(hmmp->transitions + get_mtx_index(cur_v, hmmp->nr_v-1, hmmp->nr_v)), 6);
       strcat(s, nr);
       if(fputs(strcat(s, "\n"), outfile) == EOF) {
	 perror("");
       }
#ifdef SAVE_DEBUG
       printf("%s", s);
#endif
     }
     
     /* emission probabilities */
     if(fputs("Emission probabilities 1\n", outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("saving emiss probs\n");
#endif
     alphp = hmmp->alphabet;
     strcpy(s, "\t");
     for(l = 0; l < hmmp->a_size; l++) {
       
       m = 1;
       while(*alphp != ';') {
	 s[m] = *alphp;
	 m++;
	 alphp++;
       }
       alphp++;
       s[m] =':';
       s[m+1] = ' ';
       s[m+2] ='\0';
       if(*(hmmp->emissions + get_mtx_index(cur_v, l, hmmp->a_size)) == SILENT) {
	 ftoa(nr, 0.0, 6);
       }
       else {
	 ftoa(nr, *(hmmp->emissions + get_mtx_index(cur_v, l, hmmp->a_size)), 6);
       }
       strcat(s, nr);
#ifdef SAVE_DEBUG
       printf("%s\n", s);
#endif
       strcat(s,"\n");
       if(fputs(s, outfile) == EOF) {
	 perror("");
       }
     }
     /* empty row */
     if(fputs("\n", outfile) == EOF) {
       perror("");
     }
     
#ifdef SAVE_DEBUG
     printf("reached end of module vertices loop\n");
     printf("loop index i = %d\n", i);
#endif
  
     if(hmmp->nr_alphabets > 1) {
       /* emission probabilities */
       if(fputs("Emission probabilities 2\n", outfile) == EOF) {
	 perror("");
       }
#ifdef SAVE_DEBUG
       printf("saving emiss probs\n");
#endif
       alphp = hmmp->alphabet_2;
       strcpy(s, "\t");
       for(l = 0; l < hmmp->a_size_2; l++) {
	 m = 1;
	
	 while(*alphp != ';') {
	   s[m] = *alphp;
	   m++;
	   alphp++;
	 }
	 
	 
	 alphp++;
	 s[m] =':';
	 s[m+1] = ' ';
	 s[m+2] ='\0';
	
	 if(*(hmmp->emissions_2 + get_mtx_index(cur_v, l, hmmp->a_size_2)) == SILENT) {
	   ftoa(nr, 0.0, 6);
	 }
	 else {
	   ftoa(nr, *(hmmp->emissions_2 + get_mtx_index(cur_v, l, hmmp->a_size_2)), 6);
	  
	 }
	 strcat(s, nr);
#ifdef SAVE_DEBUG
	 printf("%s\n", s);
#endif
	 strcat(s,"\n");
	 if(fputs(s, outfile) == EOF) {
	   perror("");
	 }
       }
       /* empty row */
       if(fputs("\n", outfile) == EOF) {
	 perror("");
       }
       
#ifdef SAVE_DEBUG
       printf("reached end of module vertices loop\n");
       printf("loop index i = %d\n", i);
#endif
     }

     if(hmmp->nr_alphabets > 2) {
       /* emission probabilities */
       if(fputs("Emission probabilities 3\n", outfile) == EOF) {
	 perror("");
       }
#ifdef SAVE_DEBUG
       printf("saving emiss probs\n");
#endif
       alphp = hmmp->alphabet_3;
       strcpy(s, "\t");
       for(l = 0; l < hmmp->a_size_3; l++) {
	 
	 m = 1;
	 while(*alphp != ';') {
	   s[m] = *alphp;
	   m++;
	   alphp++;
	 }
	 alphp++;
	 s[m] =':';
	 s[m+1] = ' ';
	 s[m+2] ='\0';
	 if(*(hmmp->emissions_3 + get_mtx_index(cur_v, l, hmmp->a_size_3)) == SILENT) {
	   ftoa(nr, 0.0, 6);
	 }
	 else {
	   ftoa(nr, *(hmmp->emissions_3 + get_mtx_index(cur_v, l, hmmp->a_size_3)), 6);
	 }
	 strcat(s, nr);
#ifdef SAVE_DEBUG
	 printf("%s\n", s);
#endif
	 strcat(s,"\n");
	 if(fputs(s, outfile) == EOF) {
	   perror("");
	 }
       }
       /* empty row */
       if(fputs("\n", outfile) == EOF) {
	 perror("");
       }
       
#ifdef SAVE_DEBUG
       printf("reached end of module vertices loop\n");
       printf("loop index i = %d\n", i);
#endif
     }
     
     if(hmmp->nr_alphabets > 3) {
       /* emission probabilities */
       if(fputs("Emission probabilities 4\n", outfile) == EOF) {
	 perror("");
       }
#ifdef SAVE_DEBUG
       printf("saving emiss probs\n");
#endif
       alphp = hmmp->alphabet_4;
       strcpy(s, "\t");
       for(l = 0; l < hmmp->a_size_4; l++) {
	 
	 m = 1;
	 while(*alphp != ';') {
	   s[m] = *alphp;
	   m++;
	   alphp++;
	 }
	 alphp++;
	 s[m] =':';
	 s[m+1] = ' ';
	 s[m+2] ='\0';
	 if(*(hmmp->emissions_4 + get_mtx_index(cur_v, l, hmmp->a_size_4)) == SILENT) {
	   ftoa(nr, 0.0, 6);
	 }
	 else {
	   ftoa(nr, *(hmmp->emissions_4 + get_mtx_index(cur_v, l, hmmp->a_size_4)), 6);
	 }
	 strcat(s, nr);
#ifdef SAVE_DEBUG
	 printf("%s\n", s);
#endif
	 strcat(s,"\n");
	 if(fputs(s, outfile) == EOF) {
	   perror("");
	 }
       }
       /* empty row */
       if(fputs("\n", outfile) == EOF) {
	 perror("");
       }
       
#ifdef SAVE_DEBUG
       printf("reached end of module vertices loop\n");
       printf("loop index i = %d\n", i);
#endif
     }
  }
  
  
  /* ------------------------------------- */
  if(fputs("-------------------------------------------------------\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
}
