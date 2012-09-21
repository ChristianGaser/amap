/*
 * modified version of spm_bwlabel.c 1790 2008-06-05 11:27:02Z spm $
 * from Jesper Andersson
 */

/****************************************************************
 **
 ** Set of routines implementing a 2D or 3D connected component
 ** labelling algorithm. Its interface is modelled on bwlabel
 ** (which is a routine in the image processing toolbox) and
 ** takes as input a binary image and (optionally) a connectednes
 ** criterion (6, 18 or 26, in 2D 6 will correspond to 4 and 18
 ** and 26 to 8). It will output an image/volume where each 
 ** connected component will have a unique label.
 **
 ** The implementation is not recursive (i.e. will no crash for
 ** large connected components) and is losely based on
 ** Thurfjell et al. 1992, A new three-dimensional connected
 ** components labeling algorithm with simultaneous object
 ** feature extraction capability. CVGIP: Graphical Models 
 ** and Image Processing 54(4):357-364.
 **
 ***************************************************************/

#include <math.h>
#include <limits.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>

/* Silly little macros. */

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif


/* Function prototypes. */

unsigned int do_initial_labelling(unsigned char  *bw,   /* Binary map */
                                  int            *dim,  /* Dimensions of bw */
                                  unsigned int   conn,  /* Connectivity criterion */
                                  unsigned int   *il,   /* Initially labelled map */
                                  unsigned int   **tt); /* Translation table */

unsigned int check_previous_slice(unsigned int *il,     /* Initial labelling map */
                                  unsigned int r,       /* row */
                                  unsigned int c,       /* column */
                                  unsigned int sl,      /* slice */
                                  int          dim[3],  /* dimensions of il */
                                  unsigned int conn,    /* Connectivity criterion */
                                  unsigned int *tt,     /* Translation table */
                                  unsigned int ttn);    /* Size of translation table */

void fill_tratab(unsigned int    *tt,      /* Translation table */
                 unsigned int    ttn,      /* Size of translation table */
                 unsigned int    *nabo,    /* Set of neighbours */
                 unsigned int    nr_set);  /* Number of neighbours in nabo */

int translate_labels(unsigned int    *il,     /* Map of initial labels. */
                        int             dim[3],  /* Dimensions of il. */
                        unsigned int    *tt,     /* Translation table. */
                        unsigned int    ttn,     /* Size of translation table. */
                        unsigned int    *l);     /* Final map of labels. */



/* Here starts actual code. */

unsigned int do_initial_labelling(unsigned char  *bw,   /* Binary map */
                                  int            *dim,  /* Dimensions of bw */
                                  unsigned int   conn,  /* Connectivity criterion */
                                  unsigned int   *il,   /* Initially labelled map */
                                  unsigned int   **tt)  /* Translation table */
{
   unsigned int     i = 0, j = 0;
   unsigned int     nabo[8];
   unsigned int     label = 1;
   unsigned int     nr_set = 0;
   unsigned int     l = 0;
   unsigned int     sl=0, r=0, c=0;
   unsigned int     ttn = 1000;

   *tt = (unsigned int *) malloc(ttn*sizeof(unsigned int));
   if(*tt == NULL) {
     fprintf(stderr,"Memory allocation error\n");
     exit(EXIT_FAILURE);
   }

   for (sl=0; sl<dim[2]; sl++)
   {
      for (c=0; c<dim[1]; c++)
      {
         for (r=0; r<dim[0]; r++)
         {
            nr_set = 0;
            if (bw[index(r,c,sl,dim)] > 0)
            {
               nabo[0] = check_previous_slice(il,r,c,sl,dim,conn,*tt,ttn);
                  
               if (nabo[0]) {nr_set += 1;}
               /*
                  For six(surface)-connectivity
               */
               if (conn >= 6)
               {
                  if (r)
                  {
                     if ((l = il[index(r-1,c,sl,dim)])) {nabo[nr_set++] = l;}
                  }
                  if (c)
                  {
                     if ((l = il[index(r,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
               }
               /*
                  For 18(edge)-connectivity
                  N.B. In current slice no difference to 26.
               */
               if (conn >= 18)
               {
                  if (c && r)
                  {
                     if ((l = il[index(r-1,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
                  if (c && (r < dim[0]-1))
                  {
                     if ((l = il[index(r+1,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
               }
               if (nr_set)
               {
                  il[index(r,c,sl,dim)] = nabo[0];
                  fill_tratab(*tt,ttn,nabo,nr_set);
               }
               else
               {
                  il[index(r,c,sl,dim)] = label;
                  if (label >= ttn) 
                  {
                      ttn += 1000; 
                      *tt = realloc(*tt,ttn*sizeof(unsigned int));
                  }
                  (*tt)[label-1] = label;
                  label++;
               }
            }
         }
      }
   }

   /*
      Finalise translation table
   */

   for (i=0; i<(label-1); i++)
   {
      j = i;
      while ((*tt)[j] != j+1)
      {
         j = (*tt)[j]-1;
      }
      (*tt)[i] = j+1;
   }
 
   
   return(label-1);
}

unsigned int check_previous_slice(unsigned int *il,     /* Initial labelling map */
                                  unsigned int r,       /* row */
                                  unsigned int c,       /* column */
                                  unsigned int sl,      /* slice */
                                  int          dim[3],  /* dimensions of il */
                                  unsigned int conn,    /* Connectivity criterion */
                                  unsigned int *tt,     /* Translation table */
                                  unsigned int ttn)     /* Size of translation table */
{
   unsigned int    l=0;
   unsigned int    nabo[9];
   unsigned int    nr_set = 0;

   if (!sl) return(0);
  
   if (conn >= 6)
   {
      if ((l = il[index(r,c,sl-1,dim)])) {nabo[nr_set++] = l;}
   }
   if (conn >= 18)
   {
      if (r) {if ((l = il[index(r-1,c,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (c) {if ((l = il[index(r,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (r < dim[0]-1) {if ((l = il[index(r+1,c,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (c < dim[1]-1) {if ((l = il[index(r,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
   }
   if (conn == 26)
   {
      if (r && c) {if ((l = il[index(r-1,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if ((r < dim[0]-1) && c) {if ((l = il[index(r+1,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (r && (c < dim[1]-1)) {if ((l = il[index(r-1,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if ((r < dim[0]-1) && (c < dim[1]-1)) {if ((l = il[index(r+1,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
   }

   if (nr_set) 
   {
      fill_tratab(tt,ttn,nabo,nr_set);
      return(nabo[0]);
   }
   else {return(0);}
}

void fill_tratab(unsigned int    *tt,      /* Translation table */
                 unsigned int    ttn,      /* Size of translation table */
                 unsigned int    *nabo,    /* Set of neighbours */
                 unsigned int    nr_set)   /* Number of neighbours in nabo */
{
   int           i = 0, j = 0, cntr = 0;
   unsigned int  tn[9];
   unsigned int  ltn = UINT_MAX;

   /*
   Find smallest terminal number in neighbourhood
   */

   for (i=0; i<nr_set; i++)
   {
      j = nabo[i];
      cntr=0;
      while (tt[j-1] != j) 
      {
         j = tt[j-1];
         cntr++;
         if (cntr>100) {printf("\nOoh no!!"); break;}
      }
      tn[i] = j;
      ltn = MIN(ltn,j);
   }
   /*
   Replace all terminal numbers in neighbourhood by the smallest one
   */
   for (i=0; i<nr_set; i++)
   {
      tt[tn[i]-1] = ltn;
   }

   return;
}

int translate_labels(unsigned int    *il,     /* Map of initial labels. */
                        int             dim[3],  /* Dimensions of il. */
                        unsigned int    *tt,     /* Translation table. */
                        unsigned int    ttn,     /* Size of translation table. */
                        unsigned int    *l)      /* Final map of labels. */
{
   int            n=0;
   int            i=0;
   unsigned int   ml=0;
   unsigned int   cl = 0;
   unsigned int   *fl = NULL;

   n = dim[0]*dim[1]*dim[2];

   for (i=0; i<ttn; i++) {ml = MAX(ml,tt[i]);}

   fl = (unsigned int *) malloc(ml*sizeof(unsigned int)); 
   if(fl == NULL) {
     fprintf(stderr,"Memory allocation error\n");
     exit(EXIT_FAILURE);
   }

   for (i=0; i<n; i++)
   {
      if (il[i])
      {
         if (!fl[tt[il[i]-1]-1])
         {
            cl ++; 
            fl[tt[il[i]-1]-1] = cl;
         }
         l[i] = fl[tt[il[i]-1]-1];
      }
   }

   free(fl);

   return((int)cl);
}


/* Gateway function with error check. */

void fill_cluster(unsigned char *bw, int dim[3])
{
   int            n, i, j, *count, max_count = -1, nl = 0;
   int            ind_max = 0;
   unsigned int   conn = 26;
   unsigned int   ttn = 0;
   unsigned int   *il = NULL;
   unsigned int   *tt = NULL;
   unsigned char   *bw_bak = NULL;
   unsigned int    *l = NULL;

   n = dim[0]*dim[1]*dim[2];

   /* Allocate memory for initial labelling map. */
   count  = (int *) malloc(65536*sizeof(int));
   l  = (unsigned int *) malloc(n*sizeof(unsigned int));
   il = (unsigned int *) malloc(n*sizeof(unsigned int));
   bw_bak = (unsigned char *) malloc(n*sizeof(unsigned char));

   memset(l, 0,n*sizeof(unsigned int));
   memset(il,0,n*sizeof(unsigned int));

   if((l == NULL) || (il == NULL) || (bw_bak == NULL)) {
     fprintf(stderr,"Memory allocation error\n");
     exit(EXIT_FAILURE);
   }

   /* rescue bw */
   for (i=0; i<n; i++) bw_bak[i] = bw[i];

   for (i=0; i<n; i++) if(bw[i]>0) bw[i] = 0; else bw[i] = 1;

   /* Do initial labelling and create translation table. */
   ttn = do_initial_labelling(bw,dim,conn,il,&tt);

   /* Translate labels to terminal numbers. */
   nl = translate_labels(il,dim,tt,ttn,l);

   for (j=0; j<ttn; j++) count[j] = 0;
   
   /* count cluster sizes */
   for (i=0; i<n; i++) if((int)l[i]>0) count[(int)l[i]]++;
   
   /* find index of largest cluster */
   for (j=0; j<ttn; j++) 
   {
     if (count[j] > max_count)
     {
       ind_max = j;
       max_count = count[j];
     }
   }

   for (i=0; i<n; i++) {
     if(l[i] == ind_max) bw[i] = 0;
     else bw[i] = 255;
   }

   free(il);
   free(tt);
   free(l);
   free(bw_bak);
   
   return;
}


void get_largest_cluster(unsigned char *bw, int dim[3])
{
   int            n, i, j, *count, max_count = -1, nl = 0;
   int            ind_max = 0;
   unsigned int   conn = 26;
   unsigned int   ttn = 0;
   unsigned int   *il = NULL;
   unsigned int   *tt = NULL;
   unsigned char   *bw_bak = NULL;
   unsigned int    *l = NULL, l_max = 0;

   n = dim[0]*dim[1]*dim[2];

   /* Allocate memory for initial labelling map. */
   count  = (int *) malloc(65536*sizeof(int));
   l  = (unsigned int *) malloc(n*sizeof(unsigned int));
   il = (unsigned int *) malloc(n*sizeof(unsigned int));
   bw_bak = (unsigned char *) malloc(n*sizeof(unsigned char));

   memset(l, 0,n*sizeof(unsigned int));
   memset(il,0,n*sizeof(unsigned int));

   if((l == NULL) || (il == NULL) || (bw_bak == NULL)) {
     fprintf(stderr,"Memory allocation error\n");
     exit(EXIT_FAILURE);
   }

   /* rescue bw */
   for (i=0; i<n; i++) bw_bak[i] = bw[i];

   for (i=0; i<n; i++) if(bw[i]>0) bw[i] = 1;

   /* Do initial labelling and create translation table. */
   ttn = do_initial_labelling(bw,dim,conn,il,&tt);

   /* Translate labels to terminal numbers. */
   nl = translate_labels(il,dim,tt,ttn,l);
   
   for (j=0; j<ttn; j++) count[j] = 0;
   
   /* count cluster sizes */
   for (i=0; i<n; i++) if(l[i]>0) count[l[i]]++;
   
   /* find index of largest cluster */
   for (j=0; j<ttn; j++) 
   {
     if (count[j] > max_count)
     {
       ind_max = j;
       max_count = count[j];
     }
   }

fprintf(stderr,"\n%d/%d %d\n\n",l_max,max_count,ind_max);   
   for (i=0; i<n; i++) {
     if(l[i] == ind_max) bw[i] = bw_bak[i];
     else bw[i] = 0;
   }

   free(il);
   free(tt);
   free(l);
   free(bw_bak);
   
   return;
}





