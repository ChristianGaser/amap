AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

ZLIBCFLAGS=-I$(top_srcdir)/zlib
NIFTICFLAGS=-I$(top_srcdir)/nifti
REGCFLAGS=-I$(top_srcdir)/reg-lib
EBTKSCFLAGS=-I$(top_srcdir)/EBTKS/include -I$(top_srcdir)/EBTKS/templates -I$(top_srcdir)/Splines

if EBTKS
  AM_CPPFLAGS=$(ZLIBCFLAGS) $(NIFTICFLAGS) $(REGCFLAGS) $(EBTKSCFLAGS) -DUSE_COMPMAT -DUSE_DBLMAT 
else
  AM_CPPFLAGS=$(ZLIBCFLAGS) $(NIFTICFLAGS) $(REGCFLAGS) 
endif

EXTRA_DIST = \
	autogen.sh \
	$(m4_files)

dist_bin_SCRIPTS = \
	parallelize \
	parallelize_xserves \
	reg-apps/reg_full.sh

bin_PROGRAMS = \
    niiamap \
	niiaverage \
	niismooth \
	niiornlm \
	niisanlm \
	niinucorrect \
	niiremovebridges \
	niiskullstrip \
    reg_aladin \
	reg_f3d \
	reg_resample \
	reg_tools
	
ebtks_SOURCES = \
	EBTKS/src/FileIO.cc \
	EBTKS/src/Histogram.cc \
	EBTKS/src/Path.cc \
	EBTKS/src/MString.cc \
	EBTKS/src/OrderedCltn.cc \
	EBTKS/src/dcomplex.cc \
	EBTKS/templates/Array.cc \
	EBTKS/templates/Matrix.cc \
	EBTKS/templates/MatrixSupport.cc \
	EBTKS/templates/SimpleArray.cc \
	EBTKS/templates/ValueMap.cc \
	EBTKS/templates/miscTemplateFunc.cc \
	EBTKS/templates/MatrixConv.cc \
	EBTKS/clapack/dcopy.c \
	EBTKS/clapack/dgemm.c \
	EBTKS/clapack/dgemv.c \
	EBTKS/clapack/dger.c \
	EBTKS/clapack/dlasyf.c \
	EBTKS/clapack/dscal.c \
	EBTKS/clapack/dswap.c \
	EBTKS/clapack/dsyr.c \
	EBTKS/clapack/dsysv.c \
	EBTKS/clapack/dsytf2.c \
	EBTKS/clapack/dsytrf.c \
	EBTKS/clapack/dsytrs.c \
	EBTKS/clapack/idamax.c \
	EBTKS/clapack/ieeeck.c \
	EBTKS/clapack/ilaenv.c \
	EBTKS/clapack/lsame.c \
	EBTKS/clapack/s_cmp.c \
	EBTKS/clapack/s_copy.c \
	EBTKS/clapack/xerbla.c

splines_SOURCES = \
    Splines/SplineSmooth.cc \
    Splines/TBSpline.cc \
    Splines/Spline.cc

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

niiremovebridges_SOURCES = \
	niiremovebridges.cpp \
	CRemoveBridges.cpp \
	$(nifti_SOURCES) \
    $(zlib_SOURCES) \
    ParseArgv.c \
	niilib.c
	
niiaverage_SOURCES = \
	niiaverage.c \
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
	
if EBTKS
    niiamap_SOURCES = \
	  niiamap.c \
	  Amap.c \
	  Kmeans.c \
	  Pve.c \
	  MrfPrior.c \
	  sanlm_float.c \
	  $(nifti_SOURCES) \
	  $(ebtks_SOURCES) \
      $(splines_SOURCES) \
	  niilib.c \
	  vollib.c \
	  ParseArgv.c

  niiskullstrip_SOURCES = \
      niiskullstrip.cpp \
	  affineRegistration.cpp \
	  Amap.c \
	  Kmeans.c \
	  MrfPrior.c \
	  Pve.c \
	  Bayes.c \
	  WarpPriors.c \
	  diffeo3d.c \
	  optimizer3d.c \
	  vollib.c \
	  $(nifti_SOURCES) \
	  $(ebtks_SOURCES) \
      $(splines_SOURCES) \
	  $(zlib_SOURCES) \
	  niilib.c \
	  ptr4Dlib.cpp \
	  reg-lib/_reg_resampling.cpp \
	  reg-lib/_reg_affineTransformation.cpp \
	  reg-lib/_reg_blockMatching.cpp \
	  reg-lib/_reg_tools.cpp

  niinucorrect_SOURCES = \
	  niinucorrect.c \
	  Kmeans.c \
	  $(nifti_SOURCES) \
	  $(ebtks_SOURCES) \
      $(splines_SOURCES) \
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
	  $(ebtks_SOURCES) \
      $(splines_SOURCES) \
	  niilib.c \
	  ParseArgv.c \
	  ornlm.c \
	  vollib.c 

else
  niiamap_SOURCES = \
	  niiamap.c \
	  Amap.c \
	  Kmeans.c \
	  Pve.c \
	  MrfPrior.c \
	  sanlm_float.c \
	  $(nifti_SOURCES) \
	  niilib.c \
	  vollib.c \
	  ParseArgv.c

niiskullstrip_SOURCES = \
	niiskullstrip.cpp \
	affineRegistration.cpp \
	Amap.c \
	Kmeans.c \
	MrfPrior.c \
	Pve.c \
	Bayes.c \
	WarpPriors.c \
	diffeo3d.c \
	optimizer3d.c \
	vollib.c \
	$(nifti_SOURCES) \
	$(zlib_SOURCES) \
	niilib.c \
	ptr4Dlib.cpp \
	reg-lib/_reg_resampling.cpp \
	reg-lib/_reg_affineTransformation.cpp \
	reg-lib/_reg_blockMatching.cpp \
	reg-lib/_reg_tools.cpp

  niinucorrect_SOURCES = \
	  niinucorrect.c \
	  Kmeans.c \
	  $(nifti_SOURCES) \
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

endif

ebtks_headers = \
    EBTKS/include/Complex.h \
    EBTKS/include/FileIO.h \
    EBTKS/include/Histogram.h \
    EBTKS/include/MPoint.h \
    EBTKS/include/MString.h \
    EBTKS/include/MTypes.h \
    EBTKS/include/Minc.h \
    EBTKS/include/OpTimer.h \
    EBTKS/include/OrderedCltn.h \
    EBTKS/include/Path.h \
    EBTKS/include/Polynomial.h \
    EBTKS/include/TrainingSet.h \
    EBTKS/include/amoeba.h \
    EBTKS/include/assert.h \
    EBTKS/include/backProp.h \
    EBTKS/include/dcomplex.h \
    EBTKS/include/fcomplex.h \
    EBTKS/include/matlabSupport.h \
    EBTKS/include/popen.h \
    EBTKS/include/trivials.h \
	EBTKS/templates/Array.h \
	EBTKS/templates/CachedArray.h \
	EBTKS/templates/Dictionary.h \
	EBTKS/templates/Matrix3D.h \
	EBTKS/templates/Matrix.h \
	EBTKS/templates/MatrixSupport.h \
	EBTKS/templates/MatrixTest.h \
	EBTKS/templates/miscTemplateFunc.h \
	EBTKS/templates/Pool.h \
	EBTKS/templates/SimpleArray.h \
	EBTKS/templates/Stack.h \
	EBTKS/templates/SimpleArraySpec.cc \
	EBTKS/templates/MatrixSpec.cc \
	EBTKS/clapack/f2c.h \
	EBTKS/clapack/blaswrap.h \
	EBTKS/templates/ValueMap.h

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
	$(ebtks_headers) \
	$(reg_aladin_headers) \
	$(zlib_headers) \
	$(nifti_headers) \
    niilib.h \
    niiamap.h \
    niinucorrect.h \
    niibrainmask.h \
    niiskullstrip.h \
    optimizer3d.h \
    diffeo3d.h \
    ParseArgv.h \
    Amap.h \
    CRemoveBridges.h \
    vollib.h \
	Splines/TBSpline.h \
	Splines/Spline.h 

m4_files = m4/mni_REQUIRE_LIB.m4 \
    m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4


