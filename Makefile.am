AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = \
	mincnumap \
	mincmarchingcubes \
	mincremovebridges
	
mincremovebridges_SOURCES = \
	CRemoveBridges.cpp \
	mincremovebridges.cpp

mincmarchingcubes_SOURCES = \
	mincmarchingcubes.c
	
mincnumap_SOURCES = \
	mincnumap.c \
	aniso3d.c \
	Amap.c \
	Numap.c \
	Bayes.c \
	Kmeans.c \
	Pve5.c \
	MrfPrior.c \
	Watershed3d.c \
	optimizer3d.c \
	diffeo3d.c \
	SplineSmooth.cc

m4_files = m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4
