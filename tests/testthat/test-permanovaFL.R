context("Testing `permanovaFL` function")
library(LDM)
library(testthat)

data(throat.otu.tab5)
data(throat.meta)

# test
test_that("`permanovaFL` function provides expected results", {
    res.perm <- permanovaFL(throat.otu.tab5 | (Sex+AntibioticUse) ~ SmokingStatus+PackYears, 
                            data=throat.meta, dist.method="bray", seed=82955, n.cores=1, verbose=FALSE)
    res_p <- signif(res.perm$p.permanova, 3)
    expect_equivalent(res_p, c(0.0022, 0.6450))
})
