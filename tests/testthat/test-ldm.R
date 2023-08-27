context("Testing `ldm` function")
library(LDM)
library(testthat)

data(throat.otu.tab5)
data(throat.meta)

# test
test_that("`ldm` function provides expected results", {
    res.ldm <- ldm(formula=throat.otu.tab5 | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
                   data=throat.meta, seed=67817, n.perm.max=5000, n.cores=1)
    res_global_p <- signif(res.ldm$p.global.omni, 3)
    expect_equivalent(res_global_p, c(0.0034, 0.7010))
})
