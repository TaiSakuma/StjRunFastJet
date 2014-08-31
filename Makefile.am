# Copyright (C) 2010 Tai Sakuma <sakuma@bnl.gov>

sources	= \
	StjRunFastJet.cxx

headers = \
	StjRunFastJet.h

library_includedir=$(includedir)
library_include_HEADERS = $(headers)

lib_LTLIBRARIES = libStjRunFastJet.la

libStjRunFastJet_la_SOURCES = $(sources) $(headers) StjRunFastJetLinkDef.h
libStjRunFastJet_la_CXXFLAGS = $(ROOTCFLAGS)
libStjRunFastJet_la_LDFLAGS = -version-info 1:0:0
nodist_libStjRunFastJet_la_SOURCES = StjRunFastJetDict.cxx
BUILT_SOURCES = StjRunFastJetDict.cxx
CLEANFILES = StjRunFastJetDict.cxx StjRunFastJetDict.h
StjRunFastJetDict.cxx: $(headers) StjRunFastJetLinkDef.h
	$(ROOTCINT) -f $@ -c $(CPPFLAGS) $^
