## Non-negative Matrix Factorization-based Immune-TUmor MIcroenvironment Deconvolution (NITUMID)
*Daiwei Tang & Seyoung Park*

NITUMID is the R implementation of Non-negative Matrix Factorization-based Immune-TUmor MIcroenvironment Deconvolution, a statistical framework for tumor immune microenvironment deconvolution. 

During tumor development process, its interaction with the host immune system shapes a complex tumor microenvironment (TME) consisting of tumor, various infiltrating immune cells, and other non-cancerous cell types. Interests in understanding this microenvironment has been growing since the breakthrough of immunotherapy. 

One simple *in silico* method for profiling tumor immune microenvironment is via bulk tumor gene expression data, which can be seen as a mixture of all compoenent cell types (immune cells, tumor cells, etc.) gene expression profiles. 

NITUMID takes gene expression matrix (either from RNA-Seq or microarray) as input, and use our own curated list of signature genes to estimate proportions of immune and tumor cells, as well as their corresponding mRNA levels.


Current version of NITUMID only supports tumor microenvironment deconvolution for bulk melanoma data. And current framework includes 11 component cell types: Dendritic cells, CD4+ T cell, CD8+ T cell, Macrophages, B cells, Natural killer cells (NK/NKT), Monocytes, Plasma, MAST cell, Eosinophils and melanoma. NITUMID's framework is generally applicable, the only issue limiting us from extending to more tumor types is training data, and we expect to support more tumor types in near future.

To install NITUMID, you need *devtools* package first:
```{R}
install.packages("devtools")
library("devtools")
```

Then install as follow"
```{R}
install_github("tdw1221/NITUMID")
```

For detailed introduction and tutorial, please see NITUMID's vignette (https://github.com/tdw1221/NITUMID/blob/master/docs/nitumid-tutorial.html)
