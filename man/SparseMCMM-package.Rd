\name{SparseMCMM-package}
\alias{SparseMCMM-package}
\docType{package}
\title{
SparseMCMM: Estimating and testing the microbial causal mediation effect with the high-dimensional and compositional microbiome data
}
\description{
Sparse Microbial Causal Mediation Model (SparseMCMM) is designed for the high dimensional and compositional microbiome data. SparseMCMM utilizes the linear log-contrast regression and Dirichlet regression to quantify the causal direct effect of the treatment and the causal mediation effect of the microbiome on the outcome under the counterfactual framework while addressing the compositional structure of microbiome data. Further it implements regularization techniques to handle the high-dimensional microbial mediators and identify the signature causal microbes. Furthermore, a splitting strategy (Rinaldo et al; 2019) is incorporated to account for the biases introduced by the regularization techniques employed.

SparseMCMM is particularly effective in examining the mediation effect of the microbiome within a standard three-factor (treatment, microbiome, and outcome) causal study design (Wang et al. 2020). MOverover, the analytic procedure of SparseMCMM can be harnessed to explore the influences of the microbiome on health disparities.This is depicted in an extension of the model, SparseMCMM_HD, as elucidated in Wang et al. (2023). We also discuss the differences and relevance between SparseMCMM and SparseMCMM_HD (Wang et al; 2023). It's noteworthy that the mathematical expressions of the Residual Disparity Measure (RDM), Manipulable Disparity Measure (MDM), and Overall Disparity Measure (ODM), proposed by SparseMCMM_HD, align precisely with the formulas for the Causal Direct Effect of treatment (DE), the Mediation Effect through the microbiome (ME), and the Total Effect (TE) on the outcome, as formulated in our SparseMCMM. To simplify the discussion, we will refer to these as DE, ME, and TE henceforth.

SparseMCMM consists of three components:

Component 1: Report the estimates and if applicable standard errors of DE, ME and TE respectively based on the split strategy.

Component 2: Report the point and if applicable 95\% confidence interval estimates of $ME_j$ for each microbe.

Component 3: Report the overall mediation test results: OME and CME tests.

Consequentially, SparseMCMM provides a clear and sensible causal path analysis among an exposure, compositional microbiome and outcome of interest.}
\details{
\tabular{ll}{
Package: \tab SparseMCMM\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2023-06-19\cr
License: \tab Artistic-2.0 \cr
}
}

\author{
Chan Wang and Huilin Li.

Maintainer: Chan Wang <Chan.Wang@nyulangone.org> and Huilin Li <Huilin.Li@nyulangone.org>
}
\references{
Wang C, Hu J, Blaser MJ, and Li H (2020). Estimating and testing the microbial causal mediation effect with high-dimensional and compositional microbiome data. Bioinformatics. 36(2):347-355.

Wang C, Ahn J, Tarpey T, Stella SY, Hayes RB and Li H (2023). A microbial causal mediation analytic tool for health disparity and applications in body mass index.

}
\keyword{package}
