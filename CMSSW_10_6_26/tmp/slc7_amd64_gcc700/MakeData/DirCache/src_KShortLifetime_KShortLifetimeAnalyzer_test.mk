ifeq ($(strip $(testKShortLifetimeKShortLifetimeAnalyzerTP)),)
testKShortLifetimeKShortLifetimeAnalyzerTP := self/src/KShortLifetime/KShortLifetimeAnalyzer/test
testKShortLifetimeKShortLifetimeAnalyzerTP_files := $(patsubst src/KShortLifetime/KShortLifetimeAnalyzer/test/%,%,$(foreach file,test_catch2_*.cc,$(eval xfile:=$(wildcard src/KShortLifetime/KShortLifetimeAnalyzer/test/$(file)))$(if $(xfile),$(xfile),$(warning No such file exists: src/KShortLifetime/KShortLifetimeAnalyzer/test/$(file). Please fix src/KShortLifetime/KShortLifetimeAnalyzer/test/BuildFile.))))
testKShortLifetimeKShortLifetimeAnalyzerTP_TEST_RUNNER_CMD :=  testKShortLifetimeKShortLifetimeAnalyzerTP 
testKShortLifetimeKShortLifetimeAnalyzerTP_BuildFile    := $(WORKINGDIR)/cache/bf/src/KShortLifetime/KShortLifetimeAnalyzer/test/BuildFile
testKShortLifetimeKShortLifetimeAnalyzerTP_LOC_USE := self  FWCore/TestProcessor catch2
testKShortLifetimeKShortLifetimeAnalyzerTP_PACKAGE := self/src/KShortLifetime/KShortLifetimeAnalyzer/test
ALL_PRODS += testKShortLifetimeKShortLifetimeAnalyzerTP
testKShortLifetimeKShortLifetimeAnalyzerTP_INIT_FUNC        += $$(eval $$(call Binary,testKShortLifetimeKShortLifetimeAnalyzerTP,src/KShortLifetime/KShortLifetimeAnalyzer/test,src_KShortLifetime_KShortLifetimeAnalyzer_test,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_TEST),test,$(SCRAMSTORENAME_LOGS)))
testKShortLifetimeKShortLifetimeAnalyzerTP_CLASS := TEST
else
$(eval $(call MultipleWarningMsg,testKShortLifetimeKShortLifetimeAnalyzerTP,src/KShortLifetime/KShortLifetimeAnalyzer/test))
endif
ALL_COMMONRULES += src_KShortLifetime_KShortLifetimeAnalyzer_test
src_KShortLifetime_KShortLifetimeAnalyzer_test_parent := KShortLifetime/KShortLifetimeAnalyzer
src_KShortLifetime_KShortLifetimeAnalyzer_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_KShortLifetime_KShortLifetimeAnalyzer_test,src/KShortLifetime/KShortLifetimeAnalyzer/test,TEST))
