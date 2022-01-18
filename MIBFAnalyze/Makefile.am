bin_PROGRAMS = mibfanalyze

mibfanalyze_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

mibfanalyze_CPPFLAGS = -I$(top_srcdir)/MIBFPair \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir) \
	-I$(top_srcdir)/btllib/include

mibfanalyze_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl
	
mibfanalyze_LDFLAGS = $(OPENMP_CXXFLAGS)

mibfanalyze_SOURCES = MIBFAnalyze.cpp
