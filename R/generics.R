### List of generic function definitions
# TEMPLATE #
#setGeneric(
#    '',
#    def=function(object){standardGeneric('')}
#)

### Used by Complist
setGeneric(
    'intensityReport',
    def=function(object, ...){standardGeneric('intensityReport')}
)
### Used by PepID
setGeneric(
    'evalID',
    def=function(object){standardGeneric('evalID')}
)
### Used by xcmsRaw and Complist
setGeneric(
    'plotFragment',
    def=function(object, ...){standardGeneric('plotFragment')}
)
### Used by Complist and QC
setGeneric(
    'setOutlier',
    def=function(object, ...){standardGeneric('setOutlier')}
)
### Used by Complist
setGeneric(
    'getChrom',
    def=function(object, type, Sample, outlier.rm, mix.rm, retcor, clean, objectname){standardGeneric('getChrom')}
)
### Used by Complist
setGeneric(
    'setFolder',
    def=function(object, ...){standardGeneric('setFolder')}
)
### Used by Complist
setGeneric(
    'outlier',
    def=function(object, ...){standardGeneric('outlier')}
)
### Used by Complist
setGeneric(
    'modelPCA',
    def=function(object, ...){standardGeneric('modelPCA')}
)
### Used by Complist
setGeneric(
		'modelCaret',
		def=function(object, ...){standardGeneric('modelCaret')}
)
### Used by Complist
setGeneric(
		'predictCaret',
		def=function(object, ...){standardGeneric('predictCaret')}
)
### Used by Complist
setGeneric(
    'plotScore',
    def=function(object, ...){standardGeneric('plotScore')}
)
### Used by Complist
setGeneric(
    'plotStat',
    def=function(object, ...){standardGeneric('plotStat')}
)
### Used by Complist
setGeneric(
    'editSampleinfo',
    def=function(object, ...){standardGeneric('editSampleinfo')}
)
### Used by Complist
setGeneric(
    'sampleInfo',
    def=function(object, ...){standardGeneric('sampleInfo')}
)
### Used by Complist
setGeneric(
		'sampleNames',
		def=function(object, ...){standardGeneric('sampleNames')}
)
### Used by Complist
setGeneric(
    'filterAnova',
    def=function(object, ...){standardGeneric('filterAnova')}
)
### Used by Complist
setGeneric(
		'filterNZV',
		def=function(object, ...){standardGeneric('filterNZV')}
)
### Used by Complist
setGeneric(
		'filterCor',
		def=function(object, ...){standardGeneric('filterCor')}
)
### Used by Complist
setGeneric(
    'plotCont',
    def=function(object, ...){standardGeneric('plotCont')}
)
### Used by Complist
setGeneric(
    'plotEval',
    def=function(object, ...){standardGeneric('plotEval')}
)
### Used by Complist
setGeneric(
		'plotMDS',
		def=function(object, ...){standardGeneric('plotMDS')}
)
### Used by Complist
setGeneric(
    'getCont',
    def=function(object, ...){standardGeneric('getCont')}
)
### Used by Complist
setGeneric(
    'plotData',
    def=function(object, ...){standardGeneric('plotData')}
)
### Used by Complist
setGeneric(
    'plotChromatogram',
    def=function(object, Sample, retcor, type, colour, outlier.rm, mix.rm, filter, FDR, rtwin, title, nPeptides){standardGeneric('plotChromatogram')}
)
### Used by Complist
setGeneric(
		'plotSample',
		def=function(object, ...){standardGeneric('plotSample')}
)
### Used by Complist
setGeneric(
		'plotCoverage',
		def=function(object, ...){standardGeneric('plotCoverage')}
)
### Used by Complist
setGeneric(
		'plotDetection',
		def=function(object, ...){standardGeneric('plotDetection')}
)
### Used by Complist
setGeneric(
    'plotRetcor',
    def=function(object, outlier.rm, mix.rm, rtwin, type, title){standardGeneric('plotRetcor')}
)
### Used by PepID
setGeneric(
    'mergePepID',
    def=function(object1, object2){standardGeneric('mergePepID')}
)
### Used by PepID
setGeneric(
    'plotDendro',
    def=function(object, ...){standardGeneric('plotDendro')}
)
### Used by Complist
setGeneric(
    'getRemove',
    def=function(object, ...){standardGeneric('getRemove')}
)
### Used by Complist
setGeneric(
    'findOutlier',
    def=function(object, ...){standardGeneric('findOutlier')}
)
### Used by Complist
setGeneric(
    'getFeature',
    def=function(object, ...){standardGeneric('getFeature')}
)
### Used by Complist
setGeneric(
    'plotHeat',
    def=function(object, ...){standardGeneric('plotHeat')}
)
### Used by Complist
setGeneric(
		'getPepPeakIndex',
		def=function(object, ...){standardGeneric('getPepPeakIndex')}
)
### Used by Complist
setGeneric(
		'groupReport',
		def=function(object, ...){standardGeneric('groupReport')}
)
### Used by Complist
setGeneric(
		'coverageReport',
		def=function(object, ...){standardGeneric('coverageReport')}
)
### Used by Complist
setGeneric(
		'sampleReport',
		def=function(object, ...){standardGeneric('sampleReport')}
)
### Used by Complist
setGeneric(
		'chromReport',
		def=function(object, ...){standardGeneric('chromReport')}
)
### Used by Complist
setGeneric(
		'mixReport',
		def=function(object, ...){standardGeneric('mixReport')}
)
### Used by Complist
setGeneric(
		'pcaReport',
		def=function(object, ...){standardGeneric('pcaReport')}
)
### Used by Complist
setGeneric(
		'peakReport',
		def=function(object, ...){standardGeneric('peakReport')}
)
### Used by Complist
setGeneric(
    'getRaw',
    def=function(object, ...){standardGeneric('getRaw')}
)
### Used by Complist and Peplist
setGeneric(
    'getPeplist',
    def=function(object, ...){standardGeneric('getPeplist')}
)
### Used by Complist and Peplist
setGeneric(
    'getMatch',
    def=function(object, ...){standardGeneric('getMatch')}
)
### Used by Complist
setGeneric(
    'getPeakinfo',
    def=function(object, ...){standardGeneric('getPeakinfo')}
)
### Used by Complist and Peplist
setGeneric(
    'getRawID',
    def=function(object, ...){standardGeneric('getRawID')}
)
### Used by Parameters
setGeneric(
    'editPar',
    def=function(object, ...){standardGeneric('editPar')}
)
### Used by Complist
setGeneric(
    'getBestmatch',
    def=function(object, ...){standardGeneric('getBestmatch')}
)
### Used by Complist
setGeneric(
    'reAnnotate',
    def=function(object, ...){standardGeneric('reAnnotate')}
)
### Used by Complist
setGeneric(
    'reFindPeaks',
    def=function(object, ...){standardGeneric('reFindPeaks')}
)
### Used by Complist
setGeneric(
    'reGroup',
    def=function(object, ...){standardGeneric('reGroup')}
)
### Used by Complist
setGeneric(
    'reRTcorrect',
    def=function(object, ...){standardGeneric('reRTcorrect')}
)
### Used by Complist
setGeneric(
    'reMatch',
    def=function(object, ...){standardGeneric('reMatch')}
)
### Used by Complist
setGeneric(
		'reSetFDR',
		def=function(object, ...){standardGeneric('reSetFDR')}
)
### Used by xcmsRaw
setGeneric(
    'getMSnScan',
    def=function(object, ...){standardGeneric('getMSnScan')}
)
### Used by QC
setGeneric(
    'saveQC',
    def=function(object, ...){standardGeneric('saveQC')}
)
### Used by QC
setGeneric(
    'plotQC',
    def=function(object, ...){standardGeneric('plotQC')}
)
### Used by QC
setGeneric(
    'getStatQC',
    def=function(object, ...){standardGeneric('getStatQC')}
)
### Used by QC
setGeneric(
    'modelQC',
    def=function(object, ...){standardGeneric('modelQC')}
)