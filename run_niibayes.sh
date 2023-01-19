#!/bin/bash

if [ "$1" = "" ]; then
  echo usage:  $0 T1.nii
  exit
fi

./niibayes -tpm ~/Dropbox/Template_2.nii ~/Dropbox/Template_2.nii  -target ~/Dropbox/brainweb/wTemplate_T1.nii -bias-fwhm 80 -source $1