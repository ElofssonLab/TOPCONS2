#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nets.h"

main(int argc,char *argv[])		/* Main routine */
{
  char netfile[1000];
  char datafile[1000];
  double* values;
  double* result;
  double output = 0.0;
  network net[1];   /* Number of networks to */
  FILE* data;
  int i = 0;
  int k = 0;

  if(argc==3)
    {
      strcpy(netfile,argv[1]);
      strcpy(datafile,argv[2]);
    }
  else
    {
      printf("Usage: pred <netfile> <datafile>\n");
      exit(1);
    }
  
  read_net(netfile,&net[0]);

  values = malloc(net->nin*sizeof(double));

  data = fopen(datafile, "r");
  while ( !feof(data) )
    {
      fscanf(data,"%lf",&values[i]);
      i++;
      if ( i == net->nin ) // i already incremented
	{
	  result=netfwd(values,&net[0]);
	  for (k=0;k<net->nout;k++)
	    {
	      output = result[k];
	      output = 1.0/(1 + exp(0 - output));
	      printf("%.15lf", output);
	      if ( k == net->nout - 1 ) { printf("\n"); } 
	      else {printf(" ");}
	    }
	  free(result);
	  i = 0;
	}
    }

  free(values);
  fclose(data);

}



