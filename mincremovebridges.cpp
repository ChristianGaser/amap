#include "CRemoveBridges.h"

extern "C"
{
#include  <volume_io.h>
#include  <bicpl.h>
}

void  usage(
  STRING   executable )
{  
  STRING   usage_str = "\n\
Usage: %s wm.mnc t1.mnc wm_corrected.mnc [labelvalue]\n\n\
	text.\n\n";

  print_error( usage_str, executable );
}

double Vkmeans(unsigned char *src, int nc, double *mean, int ni, int BG, int *dims)
// perform k-means algorithm give initial mean estimates    
{
  int i, j, j0, x, y, z, v;
  int count;
  long histo[256], lut[256], vol, area;
  long z_area, y_dims;
  double diff, dmin, dx, xnorm, sum, threshold;
  unsigned char *label;

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  label = (unsigned char *)malloc(sizeof(unsigned char)*vol);

  // build intensity histogram
  
  for (i = 0; i < 256; i++) histo[i] = 0;
  for (z=0;z<dims[2];z++) {
    z_area = z*area;
    for (y=0;y<dims[1];y++) {
      y_dims = y*dims[0];
      for (x=0;x<dims[0];x++) {
         v = (int)src[z_area + y_dims + x];
         if (v < BG) continue;
         if (v < 0) v = 0;
         if (v > 255) v = 255;	
         histo[v]++;
      }
    }
  }

  // use only value in histogram where cumsum is between 1..99%
  long cumsum[256];
  cumsum[0] = histo[0];
  for (i = 1; i < 256; i++) cumsum[i] = cumsum[i-1] + histo[i];
  for (i = 0; i < 256; i++) cumsum[i] = (long) round(1000.0*(double)cumsum[i]/(double)cumsum[255]);
  for (i = 0; i < 256; i++) if ((cumsum[i] <= 10) || (cumsum[i] >= 990)) histo[i] = 0;
  
  // loop through
  diff = 1e15;  count = 0;
  while (diff > 1.0 && count < ni) {

    // assign class labels
    for (i = BG; i < 256; i++) {
      dmin = 256.0 * 256.0;
      for (j = 0; j < nc; j++) {
	dx = (double) i - mean[j];
	dx *= dx;
	if (dx < dmin) {
	  lut[i] = j;
	  dmin = dx;
	}
      }
    }

    // find the new cluster centers
    diff = 0;
    for (i = 0; i < nc; i++) {
      xnorm = 0.0; sum = 0.0;
      for (j = BG; j < 256; j++)
	if (lut[j] == i) {
	  xnorm += histo[j];
	  sum +=  j * histo[j];
	}
      sum = xnorm > 0 ? sum /= xnorm: 0.0;
      dx = sum - mean[i];
      mean[i] = sum;
      dx *= dx;
      diff += dx;
    }
    count++;
  }

  // assign final labels to voxels
  for (i=BG; i<256; i++) {
    dmin = 1e15;
    j0 = 0;
    for (j = 0; j < nc; j++) {
      if (fabs((double) i - mean[j]) < dmin) {
	dmin = fabs((double)i - mean[j]);
	j0 = j;
      }
    }
    lut[i] = j0;
  }
  
  if (BG == 1) lut[0] = 0;

  // adjust for the background label
  diff = 0;
  
  for (z=0;z<dims[2];z++) {
    z_area = z*area;
    for (y=0;y<dims[1];y++) {
      y_dims = y*dims[0];
      for (x=0;x<dims[0];x++) {
         v = (int)src[z_area + y_dims + x];
         if (v >= BG) {
           if (v < 0) v = 0;
           if (v > 255) v = 255;
           label[z_area + y_dims + x] = (unsigned char)(lut[v] + BG);	
           diff += ((double)v - mean[lut[v]])*((double)v - mean[lut[v]]);
         }
         else
           label[z_area + y_dims + x] = 0;	
      }
    }
  }

  free(label);
  // return square error
  return(diff);
}

void Vtskmeans(unsigned char *src, int NI, int nclusters, int BG, int *dims, double *mu)
{
  int i, j, l, k, x, y, z;
  double e, emin, eps;

  long n[nclusters];
  double mean[nclusters];
  double var[nclusters];
  double Mu[nclusters];

  int val, nc;
  long vol, area, z_area, y_dims;

  area = dims[0]*dims[1];
  vol  = area*dims[2];
  
  int nc_initial = nclusters;

  // go through all sizes of cluster beginning with two clusters
  for (nc=2; nc<=nc_initial; nc++) {

    if (nc == 2) {
      // initialize for the two cluster case;
      n[0]=0; mean[0] = 0.0; var[0] = 0.0;


      for (z=0;z<dims[2];z++) {
        z_area = z*area;
        for (y=0;y<dims[1];y++) {
          y_dims = y*dims[0];
          for (x=0;x<dims[0];x++) {
            val = (int)src[z_area + y_dims + x];
            if (val < BG) continue;
            n[0]++;
            mean[0] += (double) val;
            var[0] += (double) val*(double) val;
          }
        }
      } 

      Mu[0] = n[0] != 0 ? mean[0]/n[0]: 0.0;
      var[0] = n[0] > 1 ? (var[0] - n[0]*Mu[0]*Mu[0])/(n[0] - 1.0) : 1.0;
      eps = 0.5*sqrt(var[0]);
    }
    else {
      // find the deviant (epsilon) for the node being divided
      eps = Mu[0];
      for (i=0; i<nc-2; i++)
        if (Mu[i+1] - Mu[i] < eps)
          eps = Mu[i+1] - Mu[i];
      if (255 - Mu[nc-2] < eps)
        eps = 255 - Mu[nc-2];
      eps = eps*0.5;
    }

    // go through low order clustering
    emin = 1e15;
    for (k=0; k<nc-1; k++) {
      for (i=nc-1; i>k+1; i--) mean[i] = Mu[i-1];
      mean[k+1] = Mu[k] + eps;  mean[k] = Mu[k] - eps;
      for (i=BG; i<k; i++) mean[i] = Mu[i];
      
      e = Vkmeans(src, nc, mean, NI, BG, dims);
            
      if (e < emin) {
        emin = e;
        for (i=0; i<nc; i++) 
          mu[i] = mean[i];
      }
    }
    for (i=0; i<nc; i++) 
      Mu[i] = mu[i]; 
  }
  
  // find the final clustering
  e = Vkmeans(src, nclusters, mu, NI, BG, dims);
  
  fprintf(stderr,"Final means: ");
  for (i=0; i<nclusters; i++) 
    fprintf(stderr,"%f ",mu[i]); 
  fprintf(stderr,"\n");    
  return;    
}

int  main(
  int   argc,
  char  *argv[] )
{
  int	x, y, z, i, dims[3], z_area, y_dims;
  long	area, vol;
  char	*axis_order[3] = { MIzspace, MIyspace, MIxspace };
  char	*arg_string;
  unsigned char	*wm, *t1, *wm_out;
  STRING	wm_filename, t1_filename, wm_out_filename;
  Volume	wm_volume, t1_volume, wm_out_volume;
  double mu[3], labelvalue;

  initialize_argument_processing( argc, argv );

  if( !get_string_argument( NULL, &wm_filename ) ||
    !get_string_argument( NULL, &t1_filename ) ||
    !get_string_argument( NULL, &wm_out_filename ) )
  {
    usage( argv[0] );
    return( 1 );
  }
  
  get_real_argument( -1.0, &labelvalue );
  
  /* read in WM image */
  if( input_volume( wm_filename, 3, File_order_dimension_names,
            NC_BYTE, FALSE, 0.0, 255.0,
            TRUE, &wm_volume, (minc_input_options *) NULL ) != OK)
  {
    print_error( "Could not read %s.\n", wm_filename );		
    return( 1 );
  }

  /* read in T1 image */
  if( input_volume( t1_filename, 3, File_order_dimension_names,
            NC_BYTE, FALSE, 0.0, 255.0,
            TRUE, &t1_volume, (minc_input_options *) NULL ) != OK)
  {
    print_error( "Could not read %s.\n", t1_filename );		
    return( 1 );
  }

  
  get_volume_sizes(wm_volume, dims);

  area = dims[0]*dims[1];
  vol  = area*dims[2];

  wm = (unsigned char *)malloc(sizeof(unsigned char)*vol);
  t1 = (unsigned char *)malloc(sizeof(unsigned char)*vol);
  wm_out = (unsigned char *)malloc(sizeof(unsigned char)*vol);

  wm_out_volume  = copy_volume(wm_volume);

  // load wm image and find bounding box
  int xmin = 255, xmax = 0, ymin = 255, ymax = 0, zmin = 255, zmax = 0; 
  for (z = 0; z < dims[2]; z++) {
    z_area = z*area;
    for (y = 0; y < dims[1]; y++) {
      y_dims = y*dims[0];
      for (x = 0; x < dims[0]; x++) {
        i = z_area + y_dims + x;
        wm[i] = get_volume_voxel_value(wm_volume,x,y,z,0,0);
        if (labelvalue != -1.0) {
          if ((double)wm[i] == labelvalue)
            wm[i] = 255;
          else wm[i] = 0;
        }
        // find bounding box of wm image
        if (wm[i] > 0) {
          xmax = (x+1 > xmax) ? x + 1 : xmax;
          ymax = (y+1 > ymax) ? y + 1 : ymax;
          zmax = (z+1 > zmax) ? z + 1 : zmax;
          xmin = (x-1 < xmin) ? x - 1 : xmin;
          ymin = (y-1 < ymin) ? y - 1 : ymin;
          zmin = (z-1 < zmin) ? z - 1 : zmin;
        }
      }
    }
  }
  fprintf(stderr,"Bounding box: x: %d-%d y: %d-%d z: %d-%d\n",xmin,xmax,ymin,ymax,zmin,zmax);
  
  // load t1 image in bounding box
  for (i = 0; i < vol; i++) t1[i] = 0;
  
  for (z = zmin; z < zmax; z++) {
    z_area = z*area;
    for (y = ymin; y < ymax; y++) {
      y_dims = y*dims[0];
      for (x = xmin; x < xmax; x++) {
        i = z_area + y_dims + x;
        t1[i] = get_volume_voxel_value(t1_volume,x,y,z,0,0);
      }
    }
  }

  // use K-Means to estimate averages for GM/WM
  Vtskmeans(t1, 25, 3, 1, dims, mu);
  
  CRemoveBridges *bro = new CRemoveBridges();

  // avgWhite,threshold,avgGray
  bro->remove(wm, t1, wm_out, dims[0],dims[1],dims[2],mu[2],((mu[2]+mu[1])/2.0),mu[1]);

  int unchanged = 225, removed = 236, added = 245, alternative = 226;
  int n_changes = 0;
  for (i = 0; i < vol; i++) {
  	if((wm_out[i] == added) || (wm_out[i] == alternative)) {
  		wm_out[i] = 255;
  		n_changes++;
  	}
  	if(wm_out[i] == unchanged)
  		wm_out[i] = 255;
  	if(wm_out[i] == removed) {
  		wm_out[i] = 0;
  		n_changes++;
  	}
  }
  fprintf(stderr,"%d voxels changed\n",n_changes);
  

  for (z = 0; z < dims[2]; z++) {
    z_area = z*area;
    for (y = 0; y < dims[1]; y++) {
      y_dims = y*dims[0];
      for (x = 0; x < dims[0]; x++) {
        i = z_area + y_dims + x;
        set_volume_voxel_value(wm_out_volume,x,y,z,0,0,wm_out[i]);
      }
    }
  }

  /* and save the modified volume */
  if (output_volume( wm_out_filename, NC_BYTE,
         FALSE, 0.0, 255.0, wm_out_volume, wm_out_filename,
         (minc_output_options *) NULL ) != OK )
      return( 1 );

  free(wm);
  free(t1);
  free(wm_out);
  
  delete_volume(wm_volume);
  delete_volume(t1_volume);
  delete_volume(wm_out_volume);

  return( 0 );
}
