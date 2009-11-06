/*
 * Christian Gaser
 * $Id$ 
 *
 */

char *mask_filename = NULL;
int n_pure_classes = 3;
int n_classes = 6;
int Niters = 200;
int subsample = 16;
int iters_nu = 40;
int pve = 1;
int correct_nu = 1;
int write_seg[3] = {0, 1, 0};
int write_nu = 0;
int write_label = 1;
double bias_fwhm = 50.0;
double thresh_brainmask = 0.05;
double thresh_kmeans = 0.5;
double weight_MRF = 1.0;
