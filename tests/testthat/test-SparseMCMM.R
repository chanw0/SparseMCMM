library(SparseMCMM)
library(testthat)


Treatment=SimulatedData$Treatment;
otu.com=SimulatedData$otu.com
outcome=SimulatedData$outcome

# test
test_that("`SparseMCMM` function provides expected results", {
  set.seed(1234)
  res=SparseMCMM(Treatment, otu.com, outcome, n.split=1, num.per=NULL)
  res_prim = res$`Esitmated Causal Effects`
  test_output = round(as.numeric(res_prim[2]), 2)
  expect_equal(test_output, 0.26)
})
