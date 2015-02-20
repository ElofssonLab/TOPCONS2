#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nets.h"


int read_net(char* filename,network* net)
{
  char	buff[512];	/* Input string */
  char	word[100];	/* PDB file line mode */
  char  tmp[100];
  FILE  *fp;
  int temp;
  long int pos;
  double tmp_in;
  int i=0;
  int j=0;
  int k=0;
  int get_w1=0;
  int get_b1=0;
  int get_w2=0;
  int get_b2=0;
  char a;

  
  fp=fopen(filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      while(fscanf(fp,"%s",buff)!=EOF)
      {
	
	if(strcmp("nin",buff)==0)
	  {
	    fscanf(fp,"%d", &net->nin);
	    //	    printf("net->nin: %d\n", net->nin);
	  }
	  else if(strcmp("nhidden",buff)==0)
	    {
	      fscanf(fp,"%d",&net->nhidden);
	      
	    }
	  else if(strcmp("nout",buff)==0)
	    {
	      fscanf(fp,"%d",&net->nout);
	    }
	  else if(strcmp("w1",buff)==0)
	    {
	      //Store values in nhidden x nin matrix
	      for(j=0;j<net->nin;j++)
		{
		  for(i=0;i<net->nhidden;i++)
		    {
		      fscanf(fp,"%lf",&net->w1[i][j]); 
		    }
		}
	    }
	  else if(strcmp("b1",buff)==0)
	    {
	      for(i=0;i<net->nhidden;i++)
		{
		  fscanf(fp,"%lf",&net->b1[i]);
		  //printf("%lf ",net->b1[i]);
		  //		}
		}
	    }
	  else if(strcmp("w2",buff)==0)
	    {
	      for(j=0;j<net->nhidden;j++)
		{
		  for (k=0;k<net->nout;k++)
		    {
		      fscanf(fp,"%lf",&net->w2[j][k]);
		    }
		}
	    }
	  else if(strcmp("b2",buff)==0)
	    {
	      for(k=0;k<net->nout;k++)
		{
		  fscanf(fp,"%lf",&net->b2[k]);
		}
	    }
//	  //zprintf("%s\n",buff);
      }
      fclose(fp);
/* 	  printf("nin %d\n",net->nin); */
/* 	  printf("nhidden %d\n",net->nhidden); */
/* 	  printf("nout %d\n", net->nout); */
/* 	  printf("w1\n"); */
/* 	  for(j=0;j<net->nin;j++) */
/* 	    { */
/* 	      for(i=0;i<net->nhidden;i++) */
/* 		{ */
/* 		  printf("%lf ", net->w1[i][j]);  */
/* 		} */
/* 	      printf("\n"); */
/* 	    } */
/* 	  printf("b1\n"); */
/* 	  for(i=0;i<net->nhidden;i++) */
/* 	    { */
/* 	      printf("%lf ", net->b1[i]);  */
/* 	    } */
/* 	  printf("\nw2\n"); */
/* 	  for(i=0;i<net->nhidden;i++) */
/* 	    { */
/* 	      printf("%lf ", net->w2[i]);  */
/* 	    } */
/* 	  printf("\nb2\n"); */
/* 	  printf("%lf\n", net->b2);  */
    }
  else
    {
      printf("Couldn't open file %s\n",filename);
      exit(1);
    }
//  for(i=0;i<net->nhidden;i++)
//    {
//	for(j=0;j<net->nin;j++)
//	  {
//	    printf("%6.5lf ",net->w1[i][j]);
//	  }
//	printf("\n");
//    }
  
  return 0;
}

double* netfwd(double* values,network* net)
{
  int i,j,k;
  double* nodsum;
  nodsum = malloc(2000*sizeof(double));
  double* output;
  output = malloc(20*sizeof(double));

  for(i=0;i<net->nhidden;i++)
    {
      nodsum[i]=0.0;
      for(j=0;j<net->nin;j++)
	{
	  nodsum[i]+=net->w1[i][j]*values[j];
	}
      nodsum[i]+=net->b1[i];
      nodsum[i]=tanh(nodsum[i]);
    }

  for (k=0;k<net->nout;k++)
    {
      output[k] = 0.0;
      for (i=0;i<net->nhidden;i++)
	{
	  output[k]+=nodsum[i]*net->w2[i][k];
	}
      output[k]+=net->b2[k];
    }

  free(nodsum);
  return output;
}

