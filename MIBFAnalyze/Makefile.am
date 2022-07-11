#bin_PROGRAMS = mibfanalyze mibfanalyzereads
bin_PROGRAMS = mibfanalyzereads

#mibfanalyze_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

#mibfanalyze_CPPFLAGS = -I$(top_srcdir)/MIBFPair \
#	-I$(top_srcdir)/Common \
#	-I$(top_srcdir) \
#	-I$(top_srcdir)/btllib/include

#mibfanalyze_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl -g
	
#mibfanalyze_LDFLAGS = $(OPENMP_CXXFLAGS)

#mibfanalyze_SOURCES = MIBFAnalyze.cpp



mibfanalyzereads_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -g

mibfanalyzereads_CPPFLAGS = -I$(top_srcdir)/MIBFPair \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir) \
	-I$(top_srcdir)/btllib/include

mibfanalyzereads_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl -lbtllib
	
mibfanalyzereads_LDFLAGS = $(OPENMP_CXXFLAGS)

mibfanalyzereads_SOURCES = MIBFAnalyzeReads.cpp

