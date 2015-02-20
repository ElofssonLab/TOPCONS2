#ifndef nets_h
#define nets_h
#define TRUE		1		/* Boolean definitions */
#define FALSE		0
#define PI		3.14159265	/* Useful constant */
#undef max
#define max(a,b)    ((a) > (b) ? (a) : (b))
#undef min
#define min(a,b)    ((a) < (b) ? (a) : (b))
#define SIZE (sizeof(a) / sizeof(a[0]));


typedef struct
{
  int nin;
  int nhidden;
  int nout;
  int nwts;
  double w1[20][2000]; //[2000];
  double b1[20];
  double w2[2000][20];
  double b2[20];
}network;


int read_net(char* filename,network* net);
double* netfwd(double* values,network* net);

//void free_net(network *net);

//typedef struct {
//  struct
//  {
//    double x,y,z;		/* Atomic coordinates */
//    double rms;		/* RMS deviation */
//    char residue[8];	/* PDB info for output */
//    char name[8];
//    int number;
//    int resnum;
//    int rescount;
//    int selected;
//  } atm[MAXATMS];
//  int CA_ref[MAXRES];
//  int res_ref[MAXRES];
//  double xcen,ycen,zcen;
//  int	atoms;			/* Current # of atoms */
//  int   residues;
//  char  sequence[MAXRES];
//  char	filename[1000];		/* filename to read molecule from */
//} molecule;
//
//typedef struct 
//{
//  double x,y,z;
//}my_vector;
//
//
//
//
//int    read_molecules(molecule *m);
//int    read_molecules_backbone(molecule *m);
//void   strncpy_NULL(char *dest, char *src, size_t n);
//int    get_atomtype3(char *name, char *res);
//int    get_atomtype(char *name, char *res);
//int    get_res6(char *res);
//int    get_res(char *res);
//double distance(molecule *m,int atomno1, int atomno2);
//void   print_type(int type_no, FILE *fp);
//int    get_res(char *res);
//void   print_res(int res,FILE *fp);
//double crd(molecule *m,int atomno1, int atomno2);   /*closest residue distance */
//double fatness(molecule *m);
//char   aa321(char *name);
//char*  assign_ss(molecule *m,float cutoff, float angle);
//int    calc_index(int size,int row,int col);
//int    hbond(molecule *m,int atomno1, int atomno2,float cutoff, float angle);
//char*  read_psipred(char* filename);
#endif
