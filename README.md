## Non-negative Matrix Factorization-based Immune-TUmor MIcroenvironment Deconvolution (NITUMID)


NITUMID is the R implementation of Non-negative Matrix Factorization-based Immune-TUmor MIcroenvironment Deconvolution, a statistical framework for tumor immune microenvironment deconvolution. 

During tumor development process, its interaction with the host immune system shapes a complex tumor microenvironment (TME) consisting of tumor, various infiltrating immune cells, and other non-cancerous cell types. Interests in understanding this microenvironment has been growing since the breakthrough of immunotherapy. 

One simple *in silico* method for profiling tumor immune microenvironment is via bulk tumor gene expression data, which can be seen as a mixture of all compoenent cell types (immune cells, tumor cells, etc.) gene expression profiles. Let's consider a gene expresison matrix $Y \in R^{m\times n}$, where $m$ rows are $m$ genes and $n$ columns are $n$ samples. Under the above framework, we can have the following:

$$
Y = W\cdot H$$
$$W \in R^{m\times k}$$
$$R \in R^{k \times n}
$$

Here $k$ is the number of component cell types, i.e., we assume each sample's gene expression profile (each column of $Y$) is a linear mixture of the $k$ component cell types ($k$ columns of $W$) with their corresponding weights ($n$ columns of $H$). Further, since the fractions of all component cell types should sum up to 1, we requires each column of $H$ sum up to 1. Notice that since $Y$, $W$ and $H$ here are all non-negative matrices, it becones essentially a non-nagative matrix factorization problem (NMF) under regularizaiton.

Current version of NITUMID only supports tumor microenvironment deconvolution for bulk melanoma data. And current framework includes 11 component cell types: Dendritic cells, CD4+ T cell, CD8+ T cell, Macrophages, B cells, Natural killer cells (NK/NKT), Monocytes, Plasma, MAST cell, Eosinophils and melanoma. NITUMID's framework is generally applicable, the only issue limiting us from extending to more tumor types is training data, and we expect to support more tumor types in near future.

Common microarray or RNA-Seq give measurements for >10000 genes' expression, most of which are not significantly different across cell types, so they would not help deconvolution but instead adding more computational burden. So it is reasonable to only choose a small subset of genes (a.k.a signature genes) whose expression levels help us differentiate among cell types.

The ideal signature genes are those who are exclusively expressing in one specific cell types. Also, if one gene is expression mostly in a few of our component cell types, it should also help. 

In our current version of NITUMID, we curated a list of 53 signature genes for our 11 cell types, their information can be seen in `signature_marker_melanoma` in this package:

```{R}
library("NITUMID")
head(signature_marker_melanoma)
```

`Cell` tells you whcih specific cell type is this gene a signature gene to. If a gene has no-NA `Secondary_cells` values, that means despite this gene's significant expression level in one cell type, in the `Secondary_cells`, it also has reasonable expression level in those secondary cell types.

Back to our NMF framework $Y=W\times H$, now we have $m=53$. Recall that the i-th column of $W$ represents the gene expression profile (of these 53 genes) in the i-th cell type, so for the signature genes for i-th cell type, their values in this column should be significantly larger. However, NMF by itself cannot achieve this kind of structure, so we introduce a guide matrix $A$ that helps us impose our desired strucutre for signature matrix. $A$ is a matrix with same dimension as $W$, and it is a trichotomous matrix consists of $1,0.5,0$.  An entry 1 for a specific cell type means the corresponding gene is highly expressed in that cell type, an entry 0 means that the gene is not expressed in the cell type, and an entry 0.5 means that the gene has an intermediate expression level in that cell type. 

**However, to make it compatible with our optimization setup, the $A$ for argument should be $\mathbb{1}-A$, in our manuscript, we defined the A that used in argument as \tilde{A}. Start now on, we will denote the guide matrix still as $A$ but the argument input A as $\tilde{A}$**

You can take a look at our current $\tilde{A}$ whose variable name is `A_melanoma_v5`

```{R}
dim(A_melanoma_v5)
A_melanoma_v5[1:3,]
```

Another issue with NMF is its instability, which comes as a price for flexibility. To get over this issue, for each of the NMF run, we use different parameters and choose the most stable setup. So for each run of NITUMID, on top of the original dataset, we will do a permutation and measure their results' consistency. NITUMID will return a consistency table as well as factorization results under different parameters, and we usually choose the factorization results under the most stable parameters.

For full details about how we set up the optimization proceudure, choose parameter and consistenct measurement, please see our manuscript for more details.