ifeq ($(strip $(PyKShortLifetimeKShortLifetimeAnalyzer)),)
PyKShortLifetimeKShortLifetimeAnalyzer := self/src/KShortLifetime/KShortLifetimeAnalyzer/python
src_KShortLifetime_KShortLifetimeAnalyzer_python_parent := 
ALL_PYTHON_DIRS += $(patsubst src/%,%,src/KShortLifetime/KShortLifetimeAnalyzer/python)
PyKShortLifetimeKShortLifetimeAnalyzer_files := $(patsubst src/KShortLifetime/KShortLifetimeAnalyzer/python/%,%,$(wildcard $(foreach dir,src/KShortLifetime/KShortLifetimeAnalyzer/python ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
PyKShortLifetimeKShortLifetimeAnalyzer_LOC_USE := self  
PyKShortLifetimeKShortLifetimeAnalyzer_PACKAGE := self/src/KShortLifetime/KShortLifetimeAnalyzer/python
ALL_PRODS += PyKShortLifetimeKShortLifetimeAnalyzer
PyKShortLifetimeKShortLifetimeAnalyzer_INIT_FUNC        += $$(eval $$(call PythonProduct,PyKShortLifetimeKShortLifetimeAnalyzer,src/KShortLifetime/KShortLifetimeAnalyzer/python,src_KShortLifetime_KShortLifetimeAnalyzer_python,1,1,$(SCRAMSTORENAME_PYTHON),$(SCRAMSTORENAME_LIB),,))
else
$(eval $(call MultipleWarningMsg,PyKShortLifetimeKShortLifetimeAnalyzer,src/KShortLifetime/KShortLifetimeAnalyzer/python))
endif
ALL_COMMONRULES += src_KShortLifetime_KShortLifetimeAnalyzer_python
src_KShortLifetime_KShortLifetimeAnalyzer_python_INIT_FUNC += $$(eval $$(call CommonProductRules,src_KShortLifetime_KShortLifetimeAnalyzer_python,src/KShortLifetime/KShortLifetimeAnalyzer/python,PYTHON))
ALL_PACKAGES += KShortLifetime/KShortLifetimeAnalyzer
subdirs_src_KShortLifetime_KShortLifetimeAnalyzer := src_KShortLifetime_KShortLifetimeAnalyzer_python src_KShortLifetime_KShortLifetimeAnalyzer_plugins src_KShortLifetime_KShortLifetimeAnalyzer_test
ALL_SUBSYSTEMS+=KShortLifetime
subdirs_src_KShortLifetime = src_KShortLifetime_KShortLifetimeAnalyzer
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
