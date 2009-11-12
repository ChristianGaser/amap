AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = \
  niiamap \
	niiaverage \
	niismooth \
	niibrainmask 
	
mincremovebridges_SOURCES = \
	CRemoveBridges.cpp \
	mincremovebridges.cpp

mincmarchingcubes_SOURCES = \
	mincmarchingcubes.c
	
niibrainmask_SOURCES = \
	niibrainmask.c \
	bwlabel.c \
	Amap.c \
	Kmeans.c \
	Pve.c \
	MrfPrior.c \
	nifti1/nifti1_io.c \
	nifti1/znzlib.c \
	niilib.c \
	ParseArgv.c \
	SplineSmooth.cc \
	vollib.c 

niiaverage_SOURCES = \
	niiaverage.c \
	nifti1/nifti1_io.c \
	nifti1/znzlib.c \
	ParseArgv.c \
	niilib.c

niismooth_SOURCES = \
	niismooth.c \
	nifti1/nifti1_io.c \
	nifti1/znzlib.c \
	vollib.c \
	ParseArgv.c \
	niilib.c

niiamap_SOURCES = \
	niiamap.c \
	Amap.c \
	Kmeans.c \
	Pve.c \
	MrfPrior.c \
	nifti1/nifti1_io.c \
	nifti1/znzlib.c \
	niilib.c \
	vollib.c \
	SplineSmooth.cc \
	ParseArgv.c

m4_files = m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4
