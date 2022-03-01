## Assess how map dimensionality affects fitted and predictive error.

[Smith et al. 2004](https://www-science-org.proxy.uchicago.edu/doi/10.1126/science.1097211?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed):
p. 3 of their supplement says:
Prediction error decreased only slightly as the number of dimensions increased,indicating (along with previous evidence that the cluster structure of the 2, 3, 4, and 5dimensional maps were the same) that distortion due to compression from higher dimensions to 2 dimensions is acceptably small for this data set.

[Bedford et al. 2014](https://elifesciences.org/articles/01914):
from p. 3 of their main text:
"By fitting the BMDS model to the training dataset, we are able to predict HI titers in the test dataset and compare these predicted titers to observed titers. We find that a two-dimensional model has better predictive power than models of lower or higher dimension in all four influenza lineages (models 1–5; Table 1). We find that this 2D model performs well, yielding an average absolute predictive error of between 0.78 and 0.91 log 2 HI titers across influenza lineages (model 2; Table 1), in line with the results of Smith et al. (2004)."

**Qs in this directory**

* By how much does map error decrease in Smith et al. if we we-fit using more than 2D?
* To what extent does the map's toplogy (projected into D1, D2) change if we infer the map in higher dimensions?
* To what extent do inferred pairwise distances change if we infer the map in higher dimensions?

**1. Check data**
`check_data.Rmd`

Smith et al., 2004 did not publish supplementary data.
However, the supplement of Bedford et al. 2014 contains data from Smith et al. 2014. 
This notebook imports the data from the Bedford et al. supplement, and checks that it contains the observations used by Bedford et al. 2014 or Smith et al. 2004 to infer maps.

