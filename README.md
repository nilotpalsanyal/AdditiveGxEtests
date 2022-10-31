# Additive Gene by Environment (GxE) Interaction Tests Under the Trend Effect of genotypes <br> <sub>Nilotpal Sanyal, Matthieu Rochemonteix, Summer Han</sub>

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation-of-cgen"
    id="toc-installation-of-cgen">Installation of CGEN</a>
-   <a href="#import-data" id="toc-import-data">Import Data</a>
-   <a href="#additive-gxe-tests-under-the-trend-of-effect-of-genotypes"
    id="toc-additive-gxe-tests-under-the-trend-of-effect-of-genotypes">Additive
    GxE tests under the trend of effect of genotypes</a>
    -   <a
        href="#without-g-e-independence-assumption-prospective-likelihood-indepfalse"
        id="toc-without-g-e-independence-assumption-prospective-likelihood-indepfalse">1.
        Without G-E independence assumption (Prospective likelihood)
        (<code>indep=FALSE</code>)</a>
    -   <a
        href="#under-g-e-independence-assumption-retrospective-likelihood-indeptrue"
        id="toc-under-g-e-independence-assumption-retrospective-likelihood-indeptrue">2.
        Under G-E independence assumption (Retrospective likelihood)
        (<code>indep=TRUE</code>)</a>
-   <a href="#references" id="toc-references">References</a>

# Introduction

This page serves as a tutorial for the additive gene by environment
(GxE) interaction tests under the trend effect of genotypes that have
been proposed by [Rochemonteix et al. (2020)](#ref1) and [Sanyal et al. (2021)](#ref2) and
implemented in the `additive.test` function of the R package `CGEN`.
Currently, these tests are available only for binary environmental
variables.

Specifically, under the trend effect of genotypes, we illustrate the
following tests:

-   Likelihood ratio tests (LRTs) without/with gene-environment
    independence assumption ([Rochemonteix et al. (2020)](#ref1)).

-   Wald tests without/with (UML/CML) gene-environment
    independence assumption, and a robust Wald test (EB) based on an empirical
    Bayes-type shrinkage estimator that combines estimates from the
    former Wald tests (UML and CML) ([Sanyal et al. (2021)](#ref2)).

The tests that use the gene-environment independence assumption are
based on the ‘retrospective likelihood’ function and constrained maximum 
likelihood (CML) estimates where the constraint is imposed by the independence 
assumption. When that assumption is not made, the use of the retrospective 
likelihood function along with the corresponding unconstrained maximum 
likelihood estimates (UML) is equivalent to the use of traditionally used 
‘propsective likelihood’ function and the estimates obtained therefrom.

We start by loading the package.

# Installation of CGEN

The R package `CGEN` can be installed from its Bioconductor repository
using:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CGEN")

library(CGEN)
```

# Import Data

For the illustration of the tests, we use the `Xdata2` dataset which
comes inside the `CGEN` package. This dataset contains sample covariate
and outcome data that have been taken from a lung cancer study. Let us
load the data and look at its contents.

``` r
data(Xdata2, package="CGEN") 

str(Xdata2)
#'data.frame':    11449 obs. of  8 variables:
# $ case.control: int  0 0 0 0 0 0 0 1 0 0 ...
# $ SNP         : int  0 0 0 1 1 0 2 1 1 0 ...
# $ smoking     : int  1 1 1 1 1 0 0 0 1 1 ...
# $ cov1        : num  25 25 15 35 5 25 5 15 15 5 ...
# $ cov2        : int  1 1 1 1 0 0 0 0 1 1 ...
# $ cov3        : int  4 8 3 4 4 4 8 0 8 4 ...
# $ cov4        : num  3.22 3.22 2.71 3.56 1.61 ...
# $ study       : int  4 4 3 5 2 4 2 3 3 2 ...
```

So, the data contain, for 11449 subjects, case-control type 0-1
outcomes, genotype values of a SNP (0,1,2), values of a binary (0-1)
environmental variable smoking, values of four covariates and a study
variable indicating which study subjects are taken from. We test for
additive SNP x smoking interaction in this data.

# Additive GxE tests under the trend of effect of genotypes

We use the `additive.test` function to conduct additive GxE interaction
tests under the assumption of an additive genetic model, i.e., the
presence of the trend effect of genotypes. The additive genetic model is
specified by the option `genetic.model = 0`. Further, the
gene-environment independence assumption is specified by setting
`indep=TRUE`. By default, `indep=FALSE` which corresponds to the
prospective likelihood-based analysis—which we illustrate first.

## 1. Without G-E independence assumption (Prospective likelihood) (`indep=FALSE`)

The following code performs tests for additive SNP x smoking interaction
under the trend effect of genotypes without any independence assumption. 
In the output, we get

-   LRT p value and test statistic under prospective likelihood.
-   Wald test p value and RERI estimates under prospective likelihood.

``` r
test_noindep <- additive.test(data = Xdata2, response.var = "case.control",
                snp.var = "SNP", exposure.var = "smoking", main.vars = c("cov1", "cov2",
                "cov3", "cov4", "study"), strata.var = "study", op = list(genetic.model = 0))

names(test_noindep)
#[1] "additive"   "model.info"

names(test_noindep$additive)
#[1] "pval.add.LRT" "pval.add.UML" "LRT.add"      "RERI.UML"    

# p value of the LRT
test_noindep$additive$pval.add.LRT
#[1] 0.03197455

# test statistic of the LRT
test_noindep$additive$LRT.add
#[1] 4.599861

# p value of the Wald (UML) test
test_noindep$additive$pval.add.UML
#[1] 0.03924117

# RERI, standard error and gradient estimates (UML)
test_noindep$additive$RERI.UML
#$value
#[1] -0.1530518
#
#$variance
#[1] 0.005511217
#
#$gradient
#     x11.UML     x2_1.UML x11:x2_1.UML 
#   0.5597640   -0.2059929    1.5068229
```

## 2. Under G-E independence assumption (Retrospective likelihood) (`indep=TRUE`)

Under the assumption that SNP and smoking are independent (specified by
`indep=TRUE`), the following code performs tests for additive SNP x
smoking interaction under the trend effect of genotypes. From the
output, we get

-   LRT p value and test statistic under retrospective likelihood.
-   Wald test p values and RERI estimates under retrospective
    likelihood.

``` r
test_indep <- additive.test(data = Xdata2, response.var = "case.control",
              snp.var = "SNP", exposure.var = "smoking", main.vars = c("cov1", "cov2",
              "cov3", "cov4", "study"), strata.var = "study", op = list(genetic.model = 0,
              indep = TRUE))

names(test_indep)
#[1] "additive"   "model.info"

names(test_indep$additive)
#[1] "pval.add.LRT" "pval.add.UML" "pval.add.CML" "pval.add.EB"  "LRT.add"      "RERI.UML"    
#[7] "RERI.CML"     "RERI.EB"     

# p value of the LRT
test_indep$additive$pval.add.LRT
#[1] 0.4595864

# test statistic of the LRT
test_indep$additive$LRT.add
#[1] 0.5469018

# p value of the Wald UML test
test_indep$additive$pval.add.UML
#[1] 0.03924117

# p value of the Wald CML test
test_indep$additive$pval.add.CML
#[1] 0.4661798

# p value of the Wald EB test
test_indep$additive$pval.add.EB
#[1] 0.131338

# RERI, standard error and gradient estimates for UML
test_indep$additive$RERI.UML
#$value
#[1] -0.1530518
#
#$variance
#[1] 0.005511217
#
#$gradient
#     x11.UML     x2_1.UML x11:x2_1.UML 
#   0.5597640   -0.2059929    1.5068229 
   
# RERI, standard error and gradient estimates for CML
test_indep$additive$RERI.CML
#$value
#[1] -0.036034
#
#$variance
#[1] 0.002445213
#
#$gradient
#     x11.CML     x2_1.CML x11:x2_1.CML 
#    0.557127    -0.119424     1.473737
    
# RERI, standard error and gradient estimates for EB
test_indep$additive$RERI.EB
#$value
#[1] -0.1194704
#
#$gradient
#     x11.UML     x2_1.UML x11:x2_1.UML      x11.CML     x2_1.CML x11:x2_1.CML 
#  0.62820373  -0.23117872   1.69105517  -0.06811736   0.01460142  -0.18018706 
#
#$variance
#[1] 0.006269466
```

# References

<div id="ref1">
De Rochemonteix, M., Napolioni, V., Sanyal, N., Belloy, M.E., Caporaso, N.E., Landi, M.T., Greicius, M.D., Chatterjee, N. and Han, S.S., 2021. A likelihood ratio test for gene-environment interaction based on the trend effect of genotype under an additive risk model using the gene-environment independence assumption. American journal of epidemiology, 190(1), pp.129-141. https://doi.org/10.1093/aje/kwaa132. <br><br>

<div id="ref2">
Sanyal, N., Napolioni, V., de Rochemonteix, M., Belloy, M.E., Caporaso, N.E., Landi, M.T., Greicius, M.D., Chatterjee, N. and Han, S.S., 2021. A Robust Test for Additive Gene-Environment Interaction Under the Trend Effect of Genotype Using an Empirical Bayes-Type Shrinkage Estimator. American journal of epidemiology, 190(9), pp.1948-1960. https://doi.org/10.1093/aje/kwab124.
