#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "structs.h" /* data structures etc */
#include "funcs.h" /* function header */


//#define SAVE_DEBUG

extern int verbose;

void savemodule(FILE*, struct hmm_multi_s*, struct module_multi_s*);


savehmm(FILE *outfile, struct hmm_multi_s *hmmp)
{
  char s[1000];
  char nr[32];
  time_t t;
  int t1, t2;
  int i,j;
  struct module_s *modulep;
  struct transition_s trans, *trans_ties;
  
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
  /* alphabet */
  strcpy(s, "ALPHABET: ");
  if(fputs(strcat(s, hmmp->alphabet), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* alphabet length */
  strcpy(s, "ALPHABET LENGTH: ");
  itoa(nr, hmmp->a_size); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
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
  strcpy(s, "NR OF EMISSION PRIORFILES: ");
  itoa(nr, hmmp->nr_ed); 
  if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
  /* the PRIOR FILES */
  strcpy(s, "EMISSION PRIORFILES:");
  for(i = 0; i < hmmp->nr_ed; i++) {
    strcat(s, " ");
    strcat(s, (hmmp->emission_dirichlets)->name);
  }
  if(fputs(strcat(s, " \n"), outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
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
  for(i = 0; i < hmmp->nr_ed; i++) {
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
    modulep = *(hmmp->modules + i);
    savemodule(outfile, hmmp, modulep);
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

void savemodule(FILE *outfile, struct hmm_multi_s *hmmp, struct module_multi_s *modulep)
{
  char s[500];
  char nr[30];
  int i, j, k ,l, m;
  int nr_t, nr_et, nr_e;
  int cur_v;
  char *alphp;
  
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
  default: printf("Error: unknown type, can't save hmm\n"); exit(-1);
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
  strcpy(s, "Emission prior file: ");
  if(fputs(strcat(strcat(s, modulep->priorfile_name), "\0"), outfile) == EOF) {
    perror("");
  }
  if(fputs("\n", outfile) == EOF) {
    perror("");
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
    strcpy(s, "Emission prior scaler: ");
    ftoa(nr, *(hmmp->vertex_emiss_prior_scalers + cur_v), 6);
    strcat(s, nr);
    strcat(s, "\n");
    if(fputs(s, outfile) == EOF) {
      perror("");
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
     strcpy(s, "Nr emissions = ");
     itoa(nr, hmmp->a_size);
     if(fputs(strcat(strcat(s, nr), "\n"), outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("%s", s);
#endif
     /* transition probabilities */
     if(fputs("Transition probabilities\n", outfile) == EOF) {
       perror("");
     }
#ifdef SAVE_DEBUG
     printf("Transition probabilities\n");
     dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
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
     if(fputs("Emission probabilities\n", outfile) == EOF) {
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
  }
  
  /* ------------------------------------- */
  if(fputs("-------------------------------------------------------\n", outfile) == EOF) {
    perror("");
  }
#ifdef SAVE_DEBUG
  printf("%s", s);
#endif
}



void printxml_hmm(FILE *outfile, struct hmm_multi_s *hmmp)
{
  char s[1000];
  char nr[32];
  char nr2[32];
  time_t t;
  int t1, t2;
  int i,j;
  struct module_s *modulep;
  struct transition_s trans, *trans_ties;
  fprintf(outfile,"<hmm>\n");
  fprintf(outfile,"<name>%s</name>\n", hmmp->name);

  /* time */
  time(&t);
  fprintf(outfile,"<time-of-creation>%s</time-of-creation>\n",ctime(&t));
  fprintf(outfile,"<alphabet>%s</alphabet>\n",hmmp->alphabet);
  itoa(nr, hmmp->a_size); 
  fprintf(outfile,"<alphabet-length>%s</alphabet-length>\n",nr);
  itoa(nr, hmmp->nr_m); 
  fprintf(outfile,"<nr-of-modules>%s</nr-of-modules>\n",nr);
  itoa(nr, hmmp->nr_v); 
  fprintf(outfile,"<nr-of-vertices>%s</nr-of-vertices>\n",nr);
  itoa(nr, hmmp->nr_t); 
  fprintf(outfile,"<nr-of-transitions>%s</nr-of-transitions>\n",nr);
  itoa(nr, hmmp->nr_d); 
  fprintf(outfile,"<nr-of-distribution-groups>%s</nr-of-distribution-groups>\n",nr);
  itoa(nr, hmmp->nr_ttg); 
  fprintf(outfile,"<nr-of-transition-tie-groups>%s</nr-of-transition-tie-groups>\n",nr);
  itoa(nr, hmmp->nr_ed); 
  fprintf(outfile,"<nr-of-emission-priorfiles>%s</nr-of-emission-priorfiles>\n",nr);
  fprintf(outfile,"<emission-priorfiles>\n");

  for(i = 0; i < hmmp->nr_ed; i++) {
    fprintf(outfile,"<name>%s</name>\n", (hmmp->emission_dirichlets)->name);
  }
  fprintf(outfile,"</emission-priorfiles>\n");


  fprintf(outfile,"<nr-of-transition-priorfiles>this part is strange in the source file savehmm.c</nr-of-transition-priorfiles>\n",nr);



  //  /* nr of PRIOR FILES */
  // strcpy(s, "NR OF TRANSITION PRIORFILES: "); 
  //  if(fputs(strcat(s, "0\n"), outfile) == EOF) {
  //  perror("");
  // }
  //  /* the PRIOR FILES */
  // strcpy(s, "TRANSITION PRIORFILES:");
  // for(i = 0; i < hmmp->nr_ed; i++) {
  //   strcat(s, " ");
    //strcat(s, (hmmp->emission_dirichlets)->name);
  //  }
  // if(fputs(strcat(s, "\n"), outfile) == EOF) {
  //   perror("");
  // }

  fprintf(outfile,"<modules>\n");


  /* save the modules */
  for(i = 0; i < hmmp->nr_m; i++) {

    modulep = *(hmmp->modules + i);
    printxml_module(outfile, hmmp, modulep);


  }
 fprintf(outfile,"</modules>\n");
  

 fprintf(outfile,"<emission-distribution-groups>\n");
  /* the groups */
  j = 0;
  for(i = 0; i < hmmp->nr_d; i++) {
    itoa(nr, i+1);
    fprintf(outfile,"<group nr=\"%s\">\n",nr);


    while((t1 = *(hmmp->distrib_groups+j)) != END) {
      itoa(nr, t1);
      fprintf(outfile,"<entry>%s</entry>\n",nr);
      j++;
    }
    j++;
    fprintf(outfile,"</group>\n");
  }
 fprintf(outfile,"</emission-distribution-groups>\n");

 fprintf(outfile,"<transition-tie-groups>\n");

  /* the groups */
  j = 0;
  for(i = 0; i < hmmp->nr_ttg; i++) {
    itoa(nr, i+1);
    fprintf(outfile,"<group tie-nr=\"%s\">\n",nr);
    while(((hmmp->trans_tie_groups)+j)->from_v != END) {
      t1 = ((hmmp->trans_tie_groups)+j)->from_v;
      itoa(nr, t1);
      t1 = ((hmmp->trans_tie_groups)+j)->to_v;
      itoa(nr2, t1);
      fprintf(outfile,"<entry from-v=\"%s\" to-v=\"%s\"/>\n",nr,nr2);
      j++;
    }
    j++;
  }
 fprintf(outfile,"</transition-tie-groups>\n");
 

  fprintf(outfile,"</hmm>\n");
  fflush(outfile);
}

void printxml_module(FILE *outfile, struct hmm_multi_s *hmmp, struct module_multi_s *modulep)
{
  char s[500];
  char nr[30];
  char nr2[30];
  int i, j, k ,l, m;
  int nr_t, nr_et, nr_e;
  int cur_v;
  char *alphp;
  
  /* name */

  fprintf(outfile,"<module name=\"%s\" type=\"",modulep->name);

  /* type */
  strcpy(s, "Type: ");
  switch(modulep->type) {
  case SINGLENODE: fprintf(outfile,"Singlenode"); break;
  case CLUSTER: fprintf(outfile,"Cluster"); break;
  case HIGHWAY: fprintf(outfile,"Highway"); break;
  case U_TURN: fprintf(outfile,"U_Turn"); break;
  case SINGLELOOP: fprintf(outfile,"Singleloop");  break;
  case FORWARD_STD: fprintf(outfile,"Forward_std"); break;
  case FORWARD_ALT: fprintf(outfile,"Forward_alt"); break; 
  case PROFILE7:  fprintf(outfile,"Profile7"); break;
  case PROFILE9: fprintf(outfile,"Profile9"); break;
  default: fprintf(stderr,"Error: unknown type, can't save hmm\n"); exit(1);
  }
  fprintf(outfile,"\">\n");
  itoa(nr, modulep->nr_v); 
  fprintf(outfile,"<nr-of-vertices>%s<nr-of-vertices>\n",nr);
  fprintf(outfile,"<emission-prior-filename>%s</emission-prior-filename>\n",modulep->priorfile_name);


//  /* transition prior file */
//  strcpy(s, "Transition prior file: null\0");
//  if(fputs(s, outfile) == EOF) {
//    perror("");
//  }
//  /* empty row */
//  if(fputs("\n\n", outfile) == EOF) {
//    perror("");
//  }

  fprintf(outfile,"<vertices>\n");
  for(i = 0; i < modulep->nr_v; i++) {
    
#ifdef SAVE_DEBUG
    printf("saving vertex %d...\n", *(modulep->vertices + i));
    printf("nr_v left: %d\n", modulep->nr_v - i);
#endif



    /* Vertex nr */
    cur_v = *(modulep->vertices + i);
    itoa(nr, cur_v);
     fprintf(outfile,"<vertex nr=\"%s\" type=\"",nr);

    switch(modulep->v_type) {
    case STARTV:     fprintf(outfile,"start"); break;
    case ENDV: fprintf(outfile,"end"); break;
    case STANDARDV: fprintf(outfile,"standard"); break;
    case LOCKEDV: fprintf(outfile,"locked");  break;
    case SILENTV: fprintf(outfile,"silent"); break;
    case PROFILEV: fprintf(outfile,"%s",get_profile_vertex_type(cur_v, hmmp->silent_vertices)); break;
    default: fprintf(stderr,"Error: unknown vertex type, can't save hmm\n"); exit(1);
    }

    
    /* vertex label */

    fprintf(outfile,"\">\n<label>%c</label>\n",*(hmmp->vertex_labels + cur_v));


    /* Transition prior scaler */

    ftoa(nr, *(hmmp->vertex_trans_prior_scalers + cur_v), 6);

    fprintf(outfile,"<transition-prior-scaler>%s</transition-prior-scaler>\n",nr);
    ftoa(nr, *(hmmp->vertex_emiss_prior_scalers + cur_v), 6);
    fprintf(outfile,"<emission-prior-scaler>%s</emission-prior-scaler>\n",nr);

    /* nr transitions */

     nr_t = 0;
     for(j = 0; j < hmmp->nr_v-1; j++) {
       if(*(hmmp->transitions + get_mtx_index(cur_v, j, hmmp->nr_v)) != 0) {
	 nr_t++;
       }
     }
     itoa(nr, nr_t);
    fprintf(outfile,"<nr-of-transitions>%s</nr-of-transitions>\n",nr);


     /* nr end transitions */

     nr_et = 0;
     if(*(hmmp->transitions + get_mtx_index(cur_v, hmmp->nr_v-1, hmmp->nr_v)) != 0.0) {
       nr_et++;
     }
     itoa(nr, nr_et);

    fprintf(outfile,"<nr-end-transitions>%s</nr-end-transitions>\n",nr);

    itoa(nr, hmmp->a_size);
    fprintf(outfile,"<nr-emissions>%s</nr-emissions>\n",nr);
    fprintf(outfile,"<transition-probabilities>\n");

     /* transition probabilities */
    
     for(k = 0; k < hmmp->nr_v-1; k++) {
       if(*(hmmp->transitions + get_mtx_index(cur_v, k, hmmp->nr_v)) != 0.0) {
	 itoa(nr, k);
         fprintf(outfile,"<vertex nr=\"%s\"",nr);
	 ftoa(nr, *(hmmp->transitions + get_mtx_index(cur_v, k, hmmp->nr_v)), 6);
         fprintf(outfile," prob=\"%s\"\>\n",nr);
       }
     }
    fprintf(outfile,"</transition-probabilities>\n");


     if(*(hmmp->transitions + get_mtx_index(cur_v, hmmp->nr_v-1, hmmp->nr_v)) != 0.0) {
       itoa(nr, hmmp->nr_v - 1);
       ftoa(nr2, *(hmmp->transitions + get_mtx_index(cur_v, hmmp->nr_v-1, hmmp->nr_v)), 6);
       fprintf(outfile,"<end-transition-probabilities vertexnr=\"%s\" prob=\"%s\"/>\n",nr,nr2);
     }
    fprintf(outfile,"<emission-probabilities>\n");

     alphp = hmmp->alphabet;
     for(l = 0; l < hmmp->a_size; l++) {

    fprintf(outfile,"<letter token=\"");

       
       m = 1;
       while(*alphp != ';') {
	 fprintf(outfile,"%c",*alphp);
	 m++;
	 alphp++;
       }



       alphp++;

       if(*(hmmp->emissions + get_mtx_index(cur_v, l, hmmp->a_size)) == SILENT) {
	 ftoa(nr, 0.0, 6);
       }
       else {
	 ftoa(nr, *(hmmp->emissions + get_mtx_index(cur_v, l, hmmp->a_size)), 6);
       }
       fprintf(outfile,"\" prob=\"%s\"\>",nr);
     }

    fprintf(outfile,"</emission-probabilities>\n");

  }
  fprintf(outfile,"</module>\n");  

}
