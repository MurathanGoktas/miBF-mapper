bin_PROGRAMS = biobloommimaker

biobloommimaker_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloommimaker_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir) \
	-I$(top_srcdir)/btllib/include

biobloommimaker_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl -lbtllib
	
biobloommimaker_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloommimaker_SOURCES = BioBloomMIMaker.cpp \
	MIBFGen.hpp \
	Options.cpp Options.h
