#include <ParseArgv.h>
#include "mincnumap.h"
#include <float.h>
#include "optimizer3d.h"
#include "diffeo3d.h"

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#define BG 1

static ArgvInfo argTable[] = {
  {"-mask", ARGV_STRING, (char *) 1, (char *) &mask_filename, 
       "Prior brainmask."},
  {"-csf", ARGV_STRING, (char *) 1, (char *) &acsf_filename, 
       "CSF prior."},
  {"-gm", ARGV_STRING, (char *) 1, (char *) &agm_filename, 
       "GM prior."},
  {"-wm", ARGV_STRING, (char *) 1, (char *) &awm_filename, 
       "WM prior."},
  {"-nc", ARGV_INT, (char *) 1, (char *) &n_classes,
       "Number of classes."},
  {"-flips", ARGV_INT, (char *) 1, (char *) &Nflips,
       "Number of flips to end."},
  {"-iters", ARGV_INT, (char *) 1, (char *) &Niters,
       "Number of iterations to end."},
  {"-sub", ARGV_INT, (char *) 1, (char *) &subsample,
       "Subsampling for Amap approach."},
  {"-iters_nu", ARGV_INT, (char *) 1, (char *) &iters_nu,
       "Number of iterations for nu correction."},
  {"-no_nucorrect", ARGV_CONSTANT, (char *) FALSE, (char *) &correct_nu,
       "Do not use nu correction."},
  {"-mrf", ARGV_FLOAT, (char *) 1, (char *) &weight_MRF,
       "Weight of MRF prior (0..1)."},
  {"-thresh", ARGV_FLOAT, (char *) 1, (char *) &thresh_brainmask,
       "Threshold for prior brainmask (0..1)."},
  {"-thresh_kmeans", ARGV_FLOAT, (char *) 1, (char *) &thresh_kmeans,
       "Threshold for Kmeans algorithm (0..1)."},
  {"-reduceto3", ARGV_CONSTANT, (char *) TRUE, (char *) &reduceto3,
       "Reduce 6 or 5 classes to 3."},
  {"-write_fuzzy", ARGV_CONSTANT, (char *) TRUE, (char *) &write_fuzzy,
       "Write fuzzy segmentations."},
  {"-write_nu", ARGV_CONSTANT, (char *) TRUE, (char *) &write_nu,
       "Write nu corrected image."},
  {"-write_label", ARGV_CONSTANT, (char *) TRUE, (char *) &write_label,
       "Write label image (default)."},
  {"-nowrite_label", ARGV_CONSTANT, (char *) FALSE, (char *) &write_label,
       "Do not write label image."},
  {"-amap", ARGV_CONSTANT, (char *) TRUE, (char *) &use_amap,
       "Use Amap approach without nu correction (default)."},
  {"-numap", ARGV_CONSTANT, (char *) TRUE, (char *) &use_numap,
       "Use Numap approach with nu correction."},
  {"-watershed", ARGV_CONSTANT, (char *) TRUE, (char *) &use_watershed,
       "Use masking based on watershed."},
  {"-nowarp", ARGV_CONSTANT, (char *) FALSE, (char *) &warp_priors,
       "Use diffeomorphic warping to register prior maps."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


int main (int argc, char *argv[])
{
  Volume	volume, label_out, prob_out, nu_out, mask_in;
  Volume	csf_in, gm_in, wm_in;
  char      *input_filename, *output_filename;
  int       i, j, dims[3], thresh, thresh_kmeans_int;
  int		x, y, z, z_area, y_dims;
  long      area, vol;
  char		*axis_order[3] = { MIzspace, MIyspace, MIxspace };
  char		*arg_string, buffer[1024], *str_ptr;
  unsigned char *label, *prob, *mask, *marker, *init_mask, *priors;
  double	*src, *nu;
  double    val, max_src, min_src, separations[3];

  /* Get arguments */
  if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
    (void) fprintf(stderr, 
    "\nUsage: %s [options] in.mnc out.mnc\n",
                     argv[0]);
    (void) fprintf(stderr, 
      "       %s -help\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  if (n_classes < 2) {
    print_error( "Minimum number of classes is 2.\n");		
    return( 1 );
  }    

  // if Numap is selected and Amap is default
  if ((use_numap) && (use_amap)) {
    use_numap = TRUE;
    use_amap = FALSE;
  }

  input_filename  = argv[1];
  output_filename = argv[2];

  if (iters_nu <= 0)
    correct_nu = FALSE;

  if (Niters == 0)
    write_fuzzy = FALSE;

  if (use_amap)  fprintf(stderr,"Use Amap approach");
  if (use_numap) fprintf(stderr,"Use Numap approach");
  if (correct_nu)
    fprintf(stderr," with nu correction.\n");
  else {
    fprintf(stderr,".\n");
    iters_nu = -1;
  }

  if ((use_watershed) && (mask_filename == NULL) && (agm_filename == NULL)) {
    fprintf(stderr,"For watershed masking you have to define mask image.\n");
    use_watershed = FALSE;
  }

  // do not write nu corrected image if correction is not selected
  if (!correct_nu) write_nu = FALSE;

  // correct reduceto3 for 6 classes
  if (Niters == 0)
    if (reduceto3) {
      reduceto3 = 0;
      fprintf(stderr,"You need at least 1 iteration to reduce 5 or 6 classes.\n");
    }

  // correct reduceto3 for 6 classes
  if (reduceto3)
    if (n_classes == 6)
      reduceto3 = 2;

  // read data and scale it to 0..255
  strcpy(buffer, input_filename);
  str_ptr = strrchr(buffer, '.');
  if (input_volume(input_filename, MAX_VAR_DIMS, axis_order,
            NC_DOUBLE, FALSE, 0.0, 0.0, TRUE, &volume, NULL) != OK)
    return(1);
  
  // read mask and check for same size
  if (mask_filename != NULL) {
    if (input_volume(mask_filename, MAX_VAR_DIMS, axis_order,
              NC_BYTE, FALSE, 0.0, 255.0, TRUE, &mask_in, NULL) != OK)
      return(1);
    if( ! volumes_are_same_grid( volume, mask_in )) {
            print_error( "Mask file %s has different size or origin.\n", mask_filename );		
          return( 1 );    
    }
  }

  // read priors
  int n_priors = 0;
  if (acsf_filename != NULL) {
    n_priors++;
    if (input_volume(acsf_filename, MAX_VAR_DIMS, axis_order,
              NC_BYTE, FALSE, 0.0, 255.0, TRUE, &csf_in, NULL) != OK)
      return(1);
    if( ! volumes_are_same_grid( volume, csf_in )) {
            print_error( "CSF prior file %s has different size or origin.\n", mask_filename );		
          return( 1 );    
    }
  }
  if (agm_filename != NULL) {
    n_priors++;
    if (input_volume(agm_filename, MAX_VAR_DIMS, axis_order,
              NC_BYTE, FALSE, 0.0, 255.0, TRUE, &gm_in, NULL) != OK)
      return(1);
    if( ! volumes_are_same_grid( volume, gm_in )) {
            print_error( "GM prior file %s has different size or origin.\n", mask_filename );		
          return( 1 );    
    }
  }
  if (awm_filename != NULL) {
    n_priors++;
    if (input_volume(awm_filename, MAX_VAR_DIMS, axis_order,
              NC_BYTE, FALSE, 0.0, 255.0, TRUE, &wm_in, NULL) != OK)
      return(1);
    if( ! volumes_are_same_grid( volume, wm_in )) {
            print_error( "WM prior file %s has different size or origin.\n", mask_filename );		
          return( 1 );    
    }
  }

  // check that 3 priors were loaded and set number of classes to 3
  if (n_priors > 0) {
    n_classes = 3;
    if(n_priors != 3) {
      fprintf(stderr,"You need to define 3 prior images.\n");
      return(1);
    }
  }

  if ((n_priors != 3) && (warp_priors)) {
    fprintf(stderr,"For warping you have to define prior images.\n");
    warp_priors = FALSE;
  }


  label_out = copy_volume_definition(volume, NC_BYTE, FALSE, 0, n_classes);
  set_volume_real_range(label_out, 0, n_classes); 

  prob_out = copy_volume_definition(volume, NC_BYTE, FALSE, 0, 255.0);
  set_volume_real_range(prob_out, 0, 1.0); 

  get_volume_sizes(volume, dims);
  get_volume_separations(volume, separations);

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  label = (unsigned char *)malloc(sizeof(unsigned char)*vol);
  mask  = (unsigned char *)malloc(sizeof(unsigned char)*vol);
  src   = (double *)malloc(sizeof(double)*vol);
  prob  = (unsigned char *)malloc(sizeof(unsigned char)*vol*n_classes);
  if (n_priors > 0)
    priors = (unsigned char *)malloc(sizeof(unsigned char)*vol*(n_classes+1));
  if (correct_nu)
    nu    = (double *)malloc(sizeof(double)*vol);

  double mean[n_classes], mu[n_classes];
  for (i = 0; i < n_classes; i++)
    mu[i] = 0;
    
  // get sure that threshold for brainmask is zero if no mask is defined
  if ((mask_filename == NULL) && (n_priors == 0)) thresh_brainmask = 0.0;
  
  thresh = (int)round(255*thresh_brainmask);
  thresh_kmeans_int = (int)round(255*thresh_kmeans);

  min_src =  FLT_MAX;
  max_src = -FLT_MAX;
  
  // load source, mask and priors
  long vol2 = 2*vol;
  long vol3 = 3*vol;
  for (z = 0; z < dims[2]; z++) {
    z_area = z*area;
    for (y = 0; y < dims[1]; y++) {
      y_dims = y*dims[0];
      for (x = 0; x < dims[0]; x++) {
        i = z_area + y_dims + x;
        src[i] = get_volume_real_value(volume,x,y,z,0,0);
        // sometimes values are < -FLT_MAX due to overflow???
        if (src[i] > -FLT_MAX)
          min_src = MIN(src[i], min_src);
          max_src = MAX(src[i], max_src);
        if (n_priors > 0) {
          priors[       i] = get_volume_voxel_value(csf_in,x,y,z,0,0);
          priors[vol  + i] = get_volume_voxel_value(gm_in, x,y,z,0,0);
          priors[vol2 + i] = get_volume_voxel_value(wm_in, x,y,z,0,0);
          int sum = (int)priors[i] + (int)priors[vol+i] + (int)priors[vol2+i];
          // because of rounding issues values might be > 255
          sum = MIN(255, sum);
          priors[vol3 + i] = 255 - (unsigned char)sum;
        }
        if (mask_filename != NULL){
          mask[i] = get_volume_voxel_value(mask_in,x,y,z,0,0);
          if (correct_nu) nu[i] = 0.0;        
        } else mask[i] = 255;
      }
    }
  }

  // correct images with values < 0
  if (min_src < 0) {
    for (i=0; i<vol; i++)
      src[i] = src[i] - min_src;
  }
  
  // use bayes approach with priors or k-means for starting estimates
  if(n_priors > 0) {
    Bayes( src, label, priors, 100, separations, dims, iters_nu);
  }
  else {
    max_src = Kmeans( src, label, mask, 25, n_classes, BG, separations, dims, reduceto3, thresh, thresh_kmeans_int, iters_nu);
  }
    
  if (Niters > 0) {
    if (use_numap)
      Numap( src, label, prob, mean, n_classes, BG, Niters, Nflips, separations, dims);
    else
      Amap( src, label, prob, mean, n_classes, BG, Niters, Nflips, subsample, dims, weight_MRF);
  }
  
  if (warp_priors) {
    float *f     = (float *)malloc(sizeof(float)*vol3);
    float *g     = (float *)malloc(sizeof(float)*vol3);
    float *v     = (float *)malloc(sizeof(float)*vol3);
    float *flow  = (float *)malloc(sizeof(float)*vol3);
    float *flow1 = (float *)malloc(sizeof(float)*vol3);

    for (i = 0; i < vol3; i++) v[i] = 0.0;

    int code = 2;    // multinomial
    int loop = 3;
    int rform = 0;   // linear energy
    double lmreg = 0.01;
    static double param[3] = {1.0, 1.0, 1.0};

    struct dartel_prm* prm = (struct dartel_prm*)malloc(sizeof(struct dartel_prm)*10);
    // first three entrys of param are equal
    for (j = 0; j < loop; j++)
      for (i = 0; i < 3; i++)  prm[j].rparam[i] = param[i];

    // some entry are equal
    for (j = 0; j < loop; j++) {
      prm[j].rform = rform;
      prm[j].cycles = 3;
      prm[j].its = 3;
      prm[j].code = code;
      prm[j].lmreg = lmreg;
    }

    prm[0].rparam[3] = 4.0;   prm[0].rparam[4] = 2.0;    prm[0].rparam[5] = 1e-6; prm[0].k = 0; 
    prm[1].rparam[3] = 2.0;   prm[1].rparam[4] = 1.0;    prm[1].rparam[5] = 1e-6; prm[1].k = 0; 
    prm[2].rparam[3] = 1.0;   prm[2].rparam[4] = 0.5;    prm[2].rparam[5] = 1e-6; prm[2].k = 1; 
    prm[3].rparam[3] = 0.5;   prm[3].rparam[4] = 0.25;   prm[3].rparam[5] = 1e-6; prm[3].k = 2; 
    prm[4].rparam[3] = 0.25;  prm[4].rparam[4] = 0.125;  prm[4].rparam[5] = 1e-6; prm[4].k = 4; 
    prm[5].rparam[3] = 0.125; prm[5].rparam[4] = 0.0625; prm[5].rparam[5] = 1e-6; prm[5].k = 6; 

    // change order to gm/wm/csf
    for (i = 0; i < vol; i++) {
      f[i] = ((float)priors[i+vol])/255.0;
      g[i] = ((float)prob[i+vol])/255.0;
      f[i+vol] = ((float)priors[i+vol2])/255.0;
      g[i+vol] = ((float)prob[i+vol2])/255.0;
      f[i+vol2] = ((float)priors[i])/255.0;
      g[i+vol2] = ((float)prob[i])/255.0;
    }

    int size[4];
    for(i=0; i < 3; i++) size[i] = dims[i];
    
    // only use gm/wm
    size[3] = 2;
    int it = 0, it0, it1; 
    double ll[3]; 
    for (it0 = 0; it0 < loop; it0++) {
      int it_scratch = iteration_scratchsize((int *)size, prm[it0].code, prm[it0].k);
      float *scratch = (float *)malloc(sizeof(float)*it_scratch);
      for (it1 = 0; it1 < prm[it0].its; it1++) {
        it++;
        iteration(size, prm[it0].k, v, f, g, (float *)0, prm[it0].rform, prm[it0].rparam, prm[it0].lmreg, 
          prm[it0].cycles, prm[it0].its, prm[it0].code, flow, ll, scratch);              
        fprintf(stderr, "%02d:\t%g\t%g\t%g\t%g\n", it, ll[0], ll[1], ll[0]+ll[1], ll[2]);
        for (i = 0; i < vol3; i++) v[i] = flow[i];
      }
      free(scratch);
    }

    expdef(size, 6, -1, v, flow, flow1, (float *)0, (float *)0); 
    for (i = 0; i < vol; i++) {
      double buf[3], sum_buf = 0.0; 
      int j;
      sampn(size, f, 3, vol,
              (double)flow[i]-1.0, (double)flow[vol+i]-1.0, (double)flow[vol2+i]-1.0, buf);
      for(j=0; j<3; j++) sum_buf += buf[j];
      // set result to zero if sum of probabilities is too low
      if (sum_buf < 0.25) {
        label[i] = 0;
        prob[i] = 0;
        prob[i+vol] = 0;
        prob[i+vol2] = 0;
      }
    }

    free(flow1);
    free(v);
    free(f);
    free(g);
    free(flow);
  }

  if (use_watershed) {

    fprintf(stderr,"Watershed masking...\n");
  
    marker    = (unsigned char *)malloc(sizeof(unsigned char)*vol);
    init_mask = (unsigned char *)malloc(sizeof(unsigned char)*vol);

    // initialize mask image and marker
    for (i = 0; i < vol; i++) {
      if (label[i] > 0)
        init_mask[i] = (unsigned char)round(255.0*src[i]/max_src);
      else
        init_mask[i] = 0;
      if ((mask[i] >= 254) && (label[i] > 1))
        marker[i] = 1;
      else
        marker[i] = 0;
    }
  
    // set marker at border to 2
    for (z = 0; z < 2; z++) 
      for (y = 0; y < dims[1]; y++) 
        for (x = 0; x < dims[0]; x++) 
          marker[z*area + y*dims[0] + x] = 2;
    for (z = 0; z < dims[2]; z++) 
      for (y = 0; y < 2; y++) 
        for (x = 0; x < dims[0]; x++) 
          marker[z*area + y*dims[0] + x] = 2;
    for (z = 0; z < dims[2]; z++) 
      for (y = 0; y < dims[1]; y++) 
        for (x = 0; x < 2; x++) 
          marker[z*area + y*dims[0] + x] = 2;
    for (z = dims[2]-2; z < dims[2]; z++) 
      for (y = 0; y < dims[1]; y++) 
        for (x = 0; x < dims[0]; x++)
          marker[z*area + y*dims[0] + x] = 2;
    for (z = 0; z < dims[2]; z++) 
      for (y = dims[1]-2; y < dims[1]; y++) 
        for (x = 0; x < dims[0]; x++)
          marker[z*area + y*dims[0] + x] = 2;
    for (z = 0; z < dims[2]; z++) 
      for (y = 0; y < dims[1]; y++) 
        for (x = dims[0]-2; x < dims[0]; x++)
          marker[z*area + y*dims[0] + x] = 2;

    watershed3d(init_mask,marker,1,dims);
  
    // apply mask
    for (i = 0; i < vol; i++) {
      if (init_mask[i] == 0)
        label[i] = 0;
    }

    free(marker);
    free(init_mask);
  
  }

  // reduce classes to 3 if neccessary
  if (reduceto3) {
    if (n_classes == 6) {
      fprintf(stderr,"Reduce 6 classes to 3.\n");
      label = (unsigned char *)Prob6to3(prob, label, mean, BG, dims);
      n_classes = 3;
    }
    if (n_classes == 5) {
      fprintf(stderr,"Reduce 5 classes to 3.\n");
      label = (unsigned char *)Prob5to3(prob, label, mean, BG, dims);
      n_classes = 3;
    }
  }

  //  copy values to volume
  if (write_label || write_nu) {
    for (z = 0; z < dims[2]; z++) {
      z_area = z*area;
      for (y = 0; y < dims[1]; y++) {
        y_dims = y*dims[0];
        for (x = 0; x < dims[0]; x++) {
          if (write_label)
            set_volume_voxel_value(label_out,x,y,z,0,0,label[z_area + y_dims + x]);
          if (write_nu)
            set_volume_real_value(volume,x,y,z,0,0,src[z_area + y_dims + x]);
        }
      }
    }
  }

  // write labeled volume
  if (write_label) {
    if (output_volume( output_filename, NC_BYTE,
           FALSE, 0.0, 255.0, label_out, input_filename,
           (minc_output_options *) NULL ) != OK )
        return( 1 );
  }
  
  // cut extension
  strcpy(output_filename+strlen(output_filename)-4,"");

  // write nu corrected image
  if (write_nu) {
    (void) sprintf( buffer, "%s_nu.mnc",output_filename);
    if (output_volume(buffer, NC_DOUBLE,
           FALSE, 0.0, 0.0, volume, input_filename,
           (minc_output_options *) NULL ) != OK )
        return( 1 );
  }
  
  // write fuzzy segmentations for each class
  if (write_fuzzy) {
    for (i = 0; i<n_classes; i++) {
      for (z = 0; z < dims[2]; z++) {
        z_area = z*area;
        for (y = 0; y < dims[1]; y++) {
          y_dims = y*dims[0];
          for (x = 0; x < dims[0]; x++) {
            set_volume_voxel_value(prob_out,x,y,z,0,0,prob[(i*vol) + z_area + y_dims + x]);
          }
        }
      }
    
      (void) sprintf( buffer, "%s_seg%d.mnc",output_filename,i);
      if (output_volume(buffer, NC_BYTE,
           FALSE, 0.0, 255.0, prob_out, input_filename,
           (minc_output_options *) NULL ) != OK )
        return( 1 );    
    }
  }

  if (n_priors > 0) {
    free(priors);
    delete_volume(csf_in);
    delete_volume(gm_in);
    delete_volume(wm_in);
  }

  free(src);
  free(prob);
  free(label);
  free(mask);
  delete_volume(volume);
  delete_volume(prob_out);
  delete_volume(label_out);
  if (mask_filename != NULL) {
    delete_volume(mask_in);
  }

  return(0);
}




















