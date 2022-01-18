bin_PROGRAMS = biobloompair

biobloompair_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloompair_CPPFLAGS = -I$(top_srcdir)/MIBFPair \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir) \
	-I$(top_srcdir)/btllib/include \
	-I$(top_srcdir)/eigen-3.3.9

biobloompair_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl
	
biobloompair_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloompair_SOURCES = MIBFPairContigs_old.cpp
