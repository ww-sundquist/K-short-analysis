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
