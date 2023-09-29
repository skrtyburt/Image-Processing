# RSEA
A collection of matlab files developed for region set enrichment analysis.

Mouse brain PET imaging must be registered to a common coordinate atlas framework and regional data extracted. Covariance matrix generation and subsequent connectomic analysis is conducted, and the output file is used as input into the region set enrichment analysis. Sets are defined as communities within a covariance as specified by an alpha level hierarchical multi-resolution consensus clustering algorithm. The functional weight of a community is compared to the weight of the same function outside the community via a 2-sample KS test. 
