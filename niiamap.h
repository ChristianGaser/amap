/*#include  <bicpl.h> */

char *mask_filename = NULL;
int n_pure_classes = 3;
int n_classes = 6;
int Niters = 200;
int subsample = 8;
int iters_nu = 40;
int iters_adf[2] = {0, 0};
int pve = 1;
int correct_nu = 1;
int write_fuzzy = 0;
int write_nu = 0;
int write_label = 1;
double thresh_brainmask = 0.05;
double thresh_kmeans = 0.5;