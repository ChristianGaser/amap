#AUTOMAKE_OPTIONS = check-news
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

ZLIBCFLAGS=-I$(top_srcdir)/zlib
NIFTICFLAGS=-I$(top_srcdir)/nifti
REGCFLAGS=-I$(top_srcdir)/reg-lib

AM_CPPFLAGS=$(ZLIBCFLAGS) $(NIFTICFLAGS) $(REGCFLAGS) 

EXTRA_DIST = \
  autogen.sh \
  $(m4_files)

dist_bin_SCRIPTS = \
  parallelize \
  parallelize_xserves \
  reg-apps/reg_full.sh

bin_PROGRAMS = \
  niibayes \
  niiamap \
  niiaverage \
  niismooth \
  niiornlm \
  niisanlm \
  niinucorrect \
  niiremovebridges \
  niiseg2mask \
  reg_aladin \
  reg_f3d \
  reg_resample \
  reg_resample \
  reg_tools

zlib_SOURCES = \
  zlib/adler32.c \
  zlib/crc32.c \
  zlib/gzio.c  \
  zlib/inffast.c \
  zlib/inftrees.c \
  zlib/trees.c \
  zlib/zutil.c \
  zlib/compress.c \
  zlib/deflate.c \
  zlib/infback.c \
  zlib/inflate.c \
  zlib/uncompr.c

nifti_SOURCES = \
  nifti/nifti1_io.c \
  nifti/znzlib.c

reg_aladin_SOURCES = \
  reg-apps/reg_aladin.cpp \
  reg-lib/_reg_resampling.cpp \
  reg-lib/_reg_affineTransformation.cpp \
  reg-lib/_reg_blockMatching.cpp \
  reg-lib/_reg_tools.cpp \
  $(zlib_SOURCES) \
  $(nifti_SOURCES) 

reg_f3d_SOURCES = \
  reg-apps/reg_f3d.cpp \
  reg-lib/_reg_resampling.cpp \
  reg-lib/_reg_bspline.cpp \
  reg-lib/_reg_bspline_comp.cpp \
  reg-lib/_reg_mutualinformation.cpp \
  reg-lib/_reg_affineTransformation.cpp \
  reg-lib/_reg_ssd.cpp \
  reg-lib/_reg_tools.cpp \
  $(zlib_SOURCES) \
  $(nifti_SOURCES)  

reg_resample_SOURCES = \
  reg-apps/reg_resample.cpp \
  reg-lib/_reg_resampling.cpp \
  reg-lib/_reg_bspline.cpp \
  reg-lib/_reg_bspline_comp.cpp \
  reg-lib/_reg_tools.cpp \
  reg-lib/_reg_affineTransformation.cpp \
  $(zlib_SOURCES) \
  $(nifti_SOURCES)  

reg_tools_SOURCES = \
  reg-apps/reg_tools.cpp \
  reg-lib/_reg_resampling.cpp \
  reg-lib/_reg_tools.cpp \
  reg-lib/_reg_affineTransformation.cpp \
  $(zlib_SOURCES) \
  $(nifti_SOURCES) 

niiaverage_SOURCES = \
  niiaverage.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  ParseArgv.c \
  niilib.c

niiremovebridges_SOURCES = \
  niiremovebridges.cpp \
  CRemoveBridges.cpp \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  ParseArgv.c \
  niilib.c
  
niismooth_SOURCES = \
  niismooth.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  vollib.c \
  ParseArgv.c \
  niilib.c

niiornlm_SOURCES = \
  niiornlm.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  ParseArgv.c \
  ornlm.c \
  niilib.c

niisanlm_SOURCES = \
  niisanlm.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  ParseArgv.c \
  sanlm_float.c \
  niilib.c
  
niiamap_SOURCES = \
  niiamap.c \
  Amap.c \
  Kmeans.c \
  Pve.c \
  MrfPrior.c \
  sanlm_float.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  niilib.c \
  vollib.c \
  ParseArgv.c

niiseg2mask_SOURCES = \
  niiseg2mask.c \
  vollib.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  niilib.c

niibayes_SOURCES = \
  niibayes.cpp \
  ptr4Dlib.cpp \
  reg-lib/_reg_resampling.cpp \
  reg-lib/_reg_affineTransformation.cpp \
  reg-lib/_reg_blockMatching.cpp \
  reg-lib/_reg_tools.cpp \
  affineRegistration.cpp \
  Bayes.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  WarpPriors.c \
  diffeo3d.c \
  optimizer3d.c \
  vollib.c \
  niilib.c

niinucorrect_SOURCES = \
  niinucorrect.c \
  Kmeans.c \
  $(nifti_SOURCES) \
  $(zlib_SOURCES) \
  niilib.c \
  vollib.c \
  ParseArgv.c

niibrainmask_SOURCES = \
  niibrainmask.c \
  bwlabel.c \
  Amap.c \
  Kmeans.c \
  Pve.c \
  MrfPrior.c \
  $(nifti_SOURCES) \
  niilib.c \
  ParseArgv.c \
  ornlm.c \
  vollib.c 

reg_aladin_headers = \
  reg-lib/_reg_affineTransformation.h \
  reg-lib/_reg_bspline.h \
  reg-lib/_reg_mutualinformation.h \
  reg-lib/_reg_ssd.h \
  reg-lib/_reg_blockMatching.h \
  reg-lib/_reg_bspline_comp.h \
  reg-lib/_reg_resampling.h \
  reg-lib/_reg_tools.h

zlib_headers = \
  zlib/crc32.h \
  zlib/deflate.h \
  zlib/inffast.h \
  zlib/inffixed.h \
  zlib/inflate.h \
  zlib/inftrees.h \
  zlib/trees.h \
  zlib/zconf.h \
  zlib/zlib.h \
  zlib/zutil.h

nifti_headers = \
  nifti/nifti1_io.h \
  nifti/nifti1_local.h \
  nifti/nifti1.h \
  nifti/znzlib.h

noinst_HEADERS = \
  $(reg_aladin_headers) \
  $(zlib_headers) \
  $(nifti_headers) \
  niilib.h \
  niiamap.h \
  niinucorrect.h \
  niibrainmask.h \
  niibayes.h \
  optimizer3d.h \
  diffeo3d.h \
  ParseArgv.h \
  Amap.h \
  CRemoveBridges.h \
  vollib.h

m4_files = m4/mni_REQUIRE_LIB.m4 \
  m4/mni_REQUIRE_MNILIBS.m4 \
  m4/mni_REQUIRE_OPENINVENTOR.m4 \
  m4/mni_cxx_have_koenig_lookup.m4 \
  m4/smr_CGAL_MAKEFILE.m4 \
  m4/smr_OPTIONAL_LIB.m4 \
  m4/smr_REQUIRED_LIB.m4 \
  m4/smr_WITH_BUILD_PATH.m4
