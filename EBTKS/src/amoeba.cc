/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1996, Alex P. Zijdenbos, 
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- 
$RCSfile: amoeba.cc,v $
$Revision: 1.1.1.1 $
$Author: jason $
$Date: 2001-11-09 16:37:24 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include <config.h>
#include "amoeba.h"
#include <assert.h>
#include "trivials.h"

#define  FLIP_RATIO      1.0
#define  CONTRACT_RATIO  0.5
#define  STRETCH_RATIO   2.0
#define  for_less( i, start, end )  for( (i) = (start);  (i) < (end);  ++(i) )

/* ----------------------------- MNI Header -----------------------------------
@NAME       : numerically_close
@INPUT      : n1
              n2
              threshold_ratio
@OUTPUT     : 
@RETURNS    : TRUE if the numbers are within the threshold ratio
@DESCRIPTION: Decides if two numbers are close to each other.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

Boolean numerically_close(
    double  n1,
    double  n2,
    double  threshold_ratio )
{
    double  avg, diff;

    diff = n1 - n2;
    if( diff < 0.0 )  diff = -diff;

    if( n1 <= threshold_ratio && n1 >= -threshold_ratio &&
        n2 <= threshold_ratio && n2 >= -threshold_ratio )
    {
        return( diff <= threshold_ratio );
    }

    avg = (n1 + n2) / 2.0;

    if( avg == 0.0 )
        return( diff <= (double) threshold_ratio );

    if( avg < 0.0 )  avg = -avg;

    return( (diff / avg) <= (double) threshold_ratio );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_function_value
@INPUT      : amoeba
              parameters
@OUTPUT     : 
@RETURNS    : function value
@DESCRIPTION: Evaluates the function being minimized by amoeba, by calling
              the user function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

double  get_function_value(
    amoeba_struct  *amoeba,
    float          parameters[] )
{
    return( (*amoeba->function) ( amoeba->function_data, parameters ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_amoeba
@INPUT      : n_parameters
              initial_parameters
              parameter_deltas
              function
              function_data
              tolerance
@OUTPUT     : amoeba
@RETURNS    : 
@DESCRIPTION: Initializes the amoeba structure to minimize the function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  initialize_amoeba(
    amoeba_struct     *amoeba,
    int               n_parameters,
    double              initial_parameters[],
    double              parameter_deltas[],
    amoeba_function   function,
    void              *function_data,
    double              tolerance )
{
    int    i, j;

    amoeba->n_parameters = n_parameters;
    amoeba->function = function;
    amoeba->function_data = function_data;
    amoeba->tolerance = tolerance;
    amoeba->n_steps_no_improvement = 0;
    typedef float * FloatPtr;
    amoeba->parameters = new FloatPtr[n_parameters + 1];
    assert(amoeba->parameters);
    for (i = 0; i < n_parameters + 1; i++) {
      amoeba->parameters[i] = new float[n_parameters];
      assert(amoeba->parameters[i]);
    }
    amoeba->values = new double[n_parameters + 1];
    amoeba->sum    = new double[n_parameters];
/*
    ALLOC2D( amoeba->parameters, n_parameters+1, n_parameters );
    ALLOC( amoeba->values, n_parameters+1 );
    ALLOC( amoeba->sum, n_parameters );
*/

    for_less( j, 0, n_parameters )
        amoeba->sum[j] = 0.0;

    for_less( i, 0, n_parameters+1 )
    {
        for_less( j, 0, n_parameters )
        {
            amoeba->parameters[i][j] = (float) initial_parameters[j];
            if( i > 0 && j == i - 1 )
                amoeba->parameters[i][j] += parameter_deltas[j];
            amoeba->sum[j] += amoeba->parameters[i][j];
        }

        amoeba->values[i] = get_function_value( amoeba, amoeba->parameters[i] );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_amoeba_parameters
@INPUT      : amoeba
@OUTPUT     : parameters
@RETURNS    : function value
@DESCRIPTION: Passes back the current position of the amoeba (best value),
              and returns the function value at that point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

double  get_amoeba_parameters(
    amoeba_struct  *amoeba,
    double           parameters[] )
{
    int   i, j, low;

    low = 0;
    for_less( i, 1, amoeba->n_parameters+1 )
    {
        if( amoeba->values[i] < amoeba->values[low] )
            low = i;
    }

    for_less( j, 0, amoeba->n_parameters )
        parameters[j] = (double) amoeba->parameters[low][j];

    return( amoeba->values[low] );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : terminate_amoeba
@INPUT      : amoeba
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Frees the amoeba.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  terminate_amoeba(
    amoeba_struct  *amoeba )
{
  for (unsigned i = 0; i < amoeba->n_parameters + 1; i++)
    delete [] amoeba->parameters[i];
  delete [] amoeba->parameters;
  delete [] amoeba->values;
  delete [] amoeba->sum;
/*
  FREE2D( amoeba->parameters );
  FREE( amoeba->values );
  FREE( amoeba->sum );
*/
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : try_amoeba
@INPUT      : amoeba
              sum
              high
              fac
@OUTPUT     : 
@RETURNS    : value
@DESCRIPTION: Does a modification to the high vertex of the amoeba and
              returns the value of the new point.  If the new point is
              better (smaller value), it replaces the high vertex of the
              amoeba.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

double  try_amoeba(
    amoeba_struct  *amoeba,
    double           sum[],
    int            high,
    double           fac )
{
    int    j;
    double   y_try, fac1, fac2;
    float  *parameters;

    parameters = new float[amoeba->n_parameters];
    assert(parameters);
/*
    ALLOC( parameters, amoeba->n_parameters );
*/
    fac1 = (1.0 - fac) / amoeba->n_parameters;
    fac2 = fac - fac1;

    for_less( j, 0, amoeba->n_parameters )
        parameters[j] = sum[j] * fac1 + amoeba->parameters[high][j] * fac2;

    y_try = get_function_value( amoeba, parameters );

    if( y_try < amoeba->values[high] )
    {
        amoeba->values[high] = y_try;
        for_less( j, 0, amoeba->n_parameters )
        {
            sum[j] += parameters[j] - amoeba->parameters[high][j];
            amoeba->parameters[high][j] = parameters[j];
        }
    }

    delete [] parameters;
//    FREE( parameters );

    return( y_try );
}

#define  N_STEPS_NO_IMPROVEMENT  20

/* ----------------------------- MNI Header -----------------------------------
@NAME       : perform_amoeba
@INPUT      : amoeba
@OUTPUT     : 
@RETURNS    : TRUE if numerically significant improvement
@DESCRIPTION: Performs one iteration of an amoeba, returning true if a
              numerically significant improvement has been found recently.
              Even if it returns FALSE, you can keep calling this function,
              since it may be contracting with no improvement, but will 
              eventually shrink small enough to get an improvment.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

Boolean  perform_amoeba(
    amoeba_struct  *amoeba )
{
    int      i, j, low, high, next_high;
    double     y_try, y_save;
    Boolean  improvement_found;

    improvement_found = TRUE;

    if( amoeba->values[0] > amoeba->values[1] )
    {
        high = 0;
        next_high = 1;
    }
    else
    {
        high = 1;
        next_high = 0;
    }

    low = next_high;

    for_less( i, 2, amoeba->n_parameters+1 )
    {
        if( amoeba->values[i] < amoeba->values[low] )
            low = i;
        else if( amoeba->values[i] > amoeba->values[high] )
        {
            next_high = high;
            high = i;
        }
        else if( amoeba->values[i] > amoeba->values[next_high] )
            next_high = i;
    }

    if( numerically_close( amoeba->values[low], amoeba->values[high],
                           amoeba->tolerance ) )
    {
        if( ++amoeba->n_steps_no_improvement >= N_STEPS_NO_IMPROVEMENT )
            improvement_found = FALSE;
    }
    else
        amoeba->n_steps_no_improvement = 0;

    y_try = try_amoeba( amoeba, amoeba->sum, high, -FLIP_RATIO );

    if( y_try <= amoeba->values[low] )
        y_try = try_amoeba( amoeba, amoeba->sum, high, STRETCH_RATIO );
    else if( y_try >= amoeba->values[next_high] )
    {
        y_save = amoeba->values[high];
        y_try = try_amoeba( amoeba, amoeba->sum, high, CONTRACT_RATIO );

        if( y_try >= y_save )
        {
            for_less( i, 0, amoeba->n_parameters+1 )
            {
                if( i != low )
                {
                    for_less( j, 0, amoeba->n_parameters )
                    {
                        amoeba->parameters[i][j] = (amoeba->parameters[i][j] +
                                            amoeba->parameters[low][j]) / 2.0;
                    }

                    amoeba->values[i] = get_function_value( amoeba,
                                                  amoeba->parameters[i] );
                }
            }

            for_less( j, 0, amoeba->n_parameters )
            {
                amoeba->sum[j] = 0.0;
                for_less( i, 0, amoeba->n_parameters+1 )
                    amoeba->sum[j] += amoeba->parameters[i][j];
            }
        }
    }

    return( improvement_found );
}
