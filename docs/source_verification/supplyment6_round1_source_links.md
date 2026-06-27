# Supplyment6 Round 1 Source Links

用途：记录本轮修改前核验过的原始论文、工具文档或数据库页面。正文中带有 `[L]` 的规则只能引用这些来源真正支持的内容；项目阈值仍应标记为 `[O]`。

## F1 QC / doublet / integration

- scDblFinder vignette: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
- Xi & Li doublet benchmark: https://doi.org/10.1016/j.cels.2020.11.008
- OSCA quality control / isOutlier: https://bioconductor.org/books/release/OSCA.basic/quality-control.html
- Harmony: https://www.nature.com/articles/s41592-019-0619-0
- scIB integration benchmark: https://www.nature.com/articles/s41592-021-01336-8
- inferCNV Bioconductor vignette: https://bioconductor.org/packages/release/bioc/vignettes/infercnv/inst/doc/inferCNV.html
- inferCNV running documentation / reference groups: https://github.com/broadinstitute/infercnv/wiki/Running-InferCNV
- CopyKAT paper (Nature Biotechnology 2021): https://www.nature.com/articles/s41587-020-00795-2
- CopyKAT open-access full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8122019/

## F2 signature and scoring

- GSE235046 GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235046
- tximport vignette: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
- limma-voom paper: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29
- AUCell vignette: https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
- SCENIC/AUCell original paper: https://www.nature.com/articles/nmeth.4463
- Ensembl Compara homology method: https://www.ensembl.org/info/genome/compara/homology_method.html
- Pseudobulk benchmark / pseudoreplication warning: https://www.nature.com/articles/s41467-021-25960-2

## F2 scoring method and mechanism modules (Supplyment8新增)

- UCell paper (NAR 2021): https://pubmed.ncbi.nlm.nih.gov/34285779/
- UCell Bioconductor documentation: https://bioconductor.org/packages/release/bioc/html/UCell.html
- Seurat AddModuleScore documentation: https://satijalab.org/seurat/reference/addmodulescore
- singscore paper: https://academic.oup.com/nar/article/46/11/e68/4973670
- singscore Bioconductor documentation: https://bioconductor.org/packages/release/bioc/html/singscore.html
- JASMINE source status: no locked official source in this manifest; use only as conditional sensitivity after execution-stage package/documentation verification.
- FerrDb V3 (NAR 2025): https://academic.oup.com/nar/article/54/D1/D572/8307356
- FerrDb V2 (NAR 2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC9825716/
- MSigDB GOBP_FERROPTOSIS: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_FERROPTOSIS
- MSigDB GOBP_LIPID_OXIDATION: https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_LIPID_OXIDATION
- Cuproptosis mechanism (Tsvetkov et al. Science 2022): https://www.science.org/doi/10.1126/science.abf0529
- REACTOME_PYROPTOSIS (MSigDB): https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_PYROPTOSIS.html
- GSDMD pyroptosis mechanism (Shi et al. Cell Research 2015): https://www.nature.com/articles/cr2015139
- Cell 2025 mitoxyperiosis original paper: https://pubmed.ncbi.nlm.nih.gov/41317732/

## F3 trajectory and SCENIC

- Monocle3 trajectory documentation: https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
- Slingshot paper: https://doi.org/10.1186/s12864-018-4772-0
- Trajectory benchmark: https://www.nature.com/articles/s41587-019-0071-9
- pySCENIC FAQ: https://pyscenic.readthedocs.io/en/latest/faq.html
- pySCENIC documentation: https://pyscenic.readthedocs.io/
- SCENIC protocol paper: https://www.nature.com/articles/s41596-020-0336-2

## F4 sample-aware LR

- edgeR User Guide: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
- edgeR glmTreat documentation: https://rdrr.io/bioc/edgeR/man/glmTreat.html
- limma dream/variancePartition documentation: https://www.bioconductor.org/packages/release/bioc/html/variancePartition.html
- CellChat paper (Nature Communications 2021): https://www.nature.com/articles/s41467-021-21246-9
- CellChat repository / documentation: https://github.com/jinworks/CellChat
- LIANA consensus framework paper (Nature Cell Biology 2022): https://www.nature.com/articles/s41556-022-01069-0
- LIANA documentation: https://saezlab.github.io/liana/

## F6 immunity

- BayesPrism paper: https://www.nature.com/articles/s43018-022-00356-3
- BayesPrism repository/docs: https://github.com/Danko-Lab/BayesPrism
- Ayers T-cell-inflamed GEP: https://www.jci.org/articles/view/91190
- Rooney CYT paper: https://www.cell.com/cell/fulltext/S0092-8674(15)00160-8
- TIDE paper: https://www.nature.com/articles/s41591-018-0136-1
- GSVA paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7
- GSVA Bioconductor documentation: https://bioconductor.org/packages/release/bioc/html/GSVA.html
- ESTIMATE paper (Nature Communications 2013): https://www.nature.com/articles/ncomms3612
- ESTIMATE package documentation: https://bioinformatics.mdanderson.org/estimate/
- MCP-counter paper (Genome Biology 2016): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5
- CIBERSORT paper: https://www.nature.com/articles/nmeth.3337
- xCell paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1
- quanTIseq paper: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0638-6
- TIMER2.0 paper: https://academic.oup.com/nar/article/48/W1/W509/5840495

## External validation data resources

- GSE239676 GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239676
- GSE239676 GEO supplementary FTP directory: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239676/suppl/
- SRP444325 SRA RunInfo endpoint: https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP444325
- ENA Portal API documentation: https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html

## F7/F8

- WGCNA FAQ mirror of official documentation: https://edo98811.github.io/WGCNA_official_documentation/faq.html
- WGCNA package paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559
- CellOracle paper: https://www.nature.com/articles/s41586-022-05688-9
- CellOracle documentation: https://morris-lab.github.io/CellOracle.documentation/
