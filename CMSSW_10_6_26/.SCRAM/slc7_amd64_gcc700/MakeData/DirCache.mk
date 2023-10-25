ifeq ($(strip $(KShortLifetimeKShortLifetimeAnalyzerAuto)),)
KShortLifetimeKShortLifetimeAnalyzerAuto := self/src/KShortLifetime/KShortLifetimeAnalyzer/plugins
PLUGINS:=yes
KShortLifetimeKShortLifetimeAnalyzerAuto_files := $(patsubst src/KShortLifetime/KShortLifetimeAnalyzer/plugins/%,%,$(wildcard $(foreach dir,src/KShortLifetime/KShortLifetimeAnalyzer/plugins ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
KShortLifetimeKShortLifetimeAnalyzerAuto_BuildFile    := $(WORKINGDIR)/cache/bf/src/KShortLifetime/KShortLifetimeAnalyzer/plugins/BuildFile
KShortLifetimeKShortLifetimeAnalyzerAuto_LOC_USE := self  FWCore/Framework FWCore/PluginManager FWCore/ParameterSet DataFormats/TrackReco DataFormats/PatCandidates FWCore/ServiceRegistry CommonTools/UtilAlgos
KShortLifetimeKShortLifetimeAnalyzerAuto_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,KShortLifetimeKShortLifetimeAnalyzerAuto,KShortLifetimeKShortLifetimeAnalyzerAuto,$(SCRAMSTORENAME_LIB),src/KShortLifetime/KShortLifetimeAnalyzer/plugins))
KShortLifetimeKShortLifetimeAnalyzerAuto_PACKAGE := self/src/KShortLifetime/KShortLifetimeAnalyzer/plugins
ALL_PRODS += KShortLifetimeKShortLifetimeAnalyzerAuto
KShortLifetime/KShortLifetimeAnalyzer_forbigobj+=KShortLifetimeKShortLifetimeAnalyzerAuto
KShortLifetimeKShortLifetimeAnalyzerAuto_INIT_FUNC        += $$(eval $$(call Library,KShortLifetimeKShortLifetimeAnalyzerAuto,src/KShortLifetime/KShortLifetimeAnalyzer/plugins,src_KShortLifetime_KShortLifetimeAnalyzer_plugins,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS),edm))
KShortLifetimeKShortLifetimeAnalyzerAuto_CLASS := LIBRARY
else
$(eval $(call MultipleWarningMsg,KShortLifetimeKShortLifetimeAnalyzerAuto,src/KShortLifetime/KShortLifetimeAnalyzer/plugins))
endif
ALL_COMMONRULES += src_KShortLifetime_KShortLifetimeAnalyzer_plugins
src_KShortLifetime_KShortLifetimeAnalyzer_plugins_parent := KShortLifetime/KShortLifetimeAnalyzer
src_KShortLifetime_KShortLifetimeAnalyzer_plugins_INIT_FUNC += $$(eval $$(call CommonProductRules,src_KShortLifetime_KShortLifetimeAnalyzer_plugins,src/KShortLifetime/KShortLifetimeAnalyzer/plugins,PLUGINS))
