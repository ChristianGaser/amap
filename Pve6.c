/*
 * Christian Gaser
 * $Id: Pve6.c 131 2009-07-24 20:59:26Z gaser $ 
 *
 */

/* This PVE calculation is a modified version from
 * the PVE software bundle:
 * Copyright (C) Jussi Tohka, Institute of Signal Processing, Tampere University of
 * Technology, 2002 - 2004.
 * P.O. Box 553, FIN-33101, Finland
 * E-mail: jussi.tohka@tut.fi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Amap.h"

void Pve6(double *src, unsigned char *prob, unsigned char *label, double *mean, double *var, int *dims, int update_label)
{
  int x,y,z,i,z_area,y_dims,ind,mxi;
  double mx, step, delta, mean2, var2, temp, val, min;
  unsigned char new_val[MAX_NC];
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
  
  step = 1.0/100.0;
    
  for (z = 1; z < dims[2]-1; z++) {
    z_area = z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x = 1; x < dims[0]-1; x++) {
        ind = z_area + y_dims + x;

        switch(label[ind]) {
        case 0: /* BG */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          break;
        case CSFLABEL+1: /* CSF */
          new_val[CSFLABEL] = 255;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          if(update_label == PVELABEL) label[ind] = ROUND(255/3);
          break;
        case GMLABEL+1: /* GM */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 255;
          new_val[WMLABEL]  = 0;
          if(update_label == PVELABEL) label[ind] = ROUND(2*255/3);
          break;
        case WMLABEL+1: /* WM */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 255;
          if(update_label == PVELABEL) label[ind] = 255;
          break;
        case BKGCSFLABEL+1: /* BKGCSF */
          min = HUGE;
          for (delta=0.0; delta <= 1.0; delta += step) {
            mean2 = (1-delta)*mean[CSFLABEL];
			      var2 = (1-delta)*(1-delta)*var[CSFLABEL];
            temp = src[ind]-mean2;
            val = ((temp*temp)/var2 +log(var2))/2.0;
            if(min>val) {
              new_val[CSFLABEL] = (unsigned char) ROUND(255.0*delta);
              new_val[GMLABEL]  = 0;
              new_val[WMLABEL]  = 0;
            }
            min = val;
          }
          if(update_label == PVELABEL) label[ind] = ROUND(255/3*delta);
          break;
        case GMCSFLABEL+1: /* GMCSF */
          min = HUGE;
          for (delta=0.0; delta <= 1.0; delta += step) {
            mean2 = delta*mean[CSFLABEL]+(1-delta)*mean[GMLABEL];
			      var2 = delta*delta*var[CSFLABEL]+(1-delta)*(1-delta)*var[GMLABEL];
            temp = src[ind]-mean2;
            val = ((temp*temp)/var2 +log(var2))/2.0;
            if(min>val) {
              new_val[CSFLABEL] = (unsigned char) ROUND(255.0*(1.0-delta));
              new_val[GMLABEL]  = (unsigned char) ROUND(255.0*delta);
              new_val[WMLABEL]  = 0;
            }
            min = val;
          }
          if(update_label == PVELABEL) label[ind] = ROUND(255/3*(1.0 + delta));
          break;
        case WMGMLABEL+1: /*WMGM */
          min = HUGE;
          for (delta=0.0; delta <= 1.0; delta += step) {
            mean2 = delta*mean[GMLABEL]+(1-delta)*mean[WMLABEL];
			      var2 = delta*delta*var[GMLABEL]+(1-delta)*(1-delta)*var[WMLABEL];
            temp = src[ind]-mean2;
            val = ((temp*temp)/var2 +log(var2))/2.0;
            if(min>val) {
              new_val[CSFLABEL] = 0;
              new_val[GMLABEL]  = (unsigned char) ROUND(255.0*(1.0-delta));
              new_val[WMLABEL]  = (unsigned char) ROUND(255.0*delta);
            }
            min = val;
          }
          if(update_label == PVELABEL) label[ind] = ROUND(255/3*(2.0 + delta));
          break;
        }

        prob[          ind] = new_val[CSFLABEL];
        prob[vol +     ind] = new_val[GMLABEL];
        prob[(2*vol) + ind] = new_val[WMLABEL];
        
        /* set old probabilities for mixed classes to zero */
        prob[(3*vol) + ind] = 0;
        prob[(4*vol) + ind] = 0;
        prob[(5*vol) + ind] = 0;
        
        /* get new label */
        if(update_label == LABEL) {
          mx = -HUGE;
          if(label[ind] > 0) {
            for (i = 0; i < 3; i++) {
              if (new_val[1+i*2] > mx) {
                mx = new_val[1+i*2];
                mxi = i;
              }
            }
            label[ind] = mxi + 1;
          }
        }
      }
    }
  }  
}

void Pve6_orig(double *src, unsigned char *prob, unsigned char *label, double *mean, int *dims, int update_label)
{
  int x,y,z,i,z_area,y_dims,ind,mxi;
  double w, mx;
  unsigned char new_val[MAX_NC];
  
  int area = dims[0]*dims[1];
  int vol = area*dims[2];
    
  for (z = 1; z < dims[2]-1; z++) {
    z_area = z*area;
    for (y = 1; y < dims[1]-1; y++) {
      y_dims = y*dims[0];
      for (x = 1; x < dims[0]-1; x++) {
        ind = z_area + y_dims + x;

        switch(label[ind]) {
        case 0: /* BG */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          break;
        case CSFLABEL+1: /* CSF */
          new_val[CSFLABEL] = 255;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          if(update_label == PVELABEL) label[ind] = ROUND(255/3);
          break;
        case GMLABEL+1: /* GM */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 255;
          new_val[WMLABEL]  = 0;
          if(update_label == PVELABEL) label[ind] = ROUND(2*255/3);
          break;
        case WMLABEL+1: /* WM */
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 255;
          if(update_label == PVELABEL) label[ind] = 255;
          break;
        case BKGCSFLABEL+1: /* BKGCSF */
          w = src[ind]/mean[CSFLABEL];
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL] = (unsigned char) ROUND(255.0*w);
          new_val[GMLABEL]  = 0;
          new_val[WMLABEL]  = 0;
          if(update_label == PVELABEL) label[ind] = ROUND(255/3*w);
          break;
        case GMCSFLABEL+1: /* GMCSF */
          w = (src[ind] - mean[CSFLABEL])/(mean[GMLABEL]-mean[CSFLABEL]);
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL] = (unsigned char) ROUND(255.0*(1-w));
          new_val[GMLABEL]  = (unsigned char) ROUND(255.0*w);
          new_val[WMLABEL]  = 0;
          if(update_label == PVELABEL) label[ind] = ROUND(255/3*(1.0 + w));
          break;
        case WMGMLABEL+1: /*WMGM */
          w = (src[ind] - mean[GMLABEL])/(mean[WMLABEL]-mean[GMLABEL]);
          if(w > 1.0) w = 1.0; if(w < 0.0) w = 0.0;
          new_val[CSFLABEL] = 0;
          new_val[GMLABEL]  = (unsigned char) ROUND(255.0*(1-w));
          new_val[WMLABEL]  = (unsigned char) ROUND(255.0*w);
          if(update_label == PVELABEL) label[ind] = ROUND(255/3*(2.0 + w));
          break;
        }

        prob[          ind] = new_val[CSFLABEL];
        prob[vol +     ind] = new_val[GMLABEL];
        prob[(2*vol) + ind] = new_val[WMLABEL];
        
        /* set old probabilities for mixed classes to zero */
        prob[(3*vol) + ind] = 0;
        prob[(4*vol) + ind] = 0;
        prob[(5*vol) + ind] = 0;
        
        /* get new label */
        if(update_label == LABEL) {
          mx = -HUGE;
          if(label[ind] > 0) {
            for (i = 0; i < 3; i++) {
              if (new_val[1+i*2] > mx) {
                mx = new_val[1+i*2];
                mxi = i;
              }
            }
            label[ind] = mxi + 1;
          }
        }
      }
    }
  }  
}
