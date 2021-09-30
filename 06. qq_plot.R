# QQ plot generation and lambda value estimation
# A. Load libraries
library(qqman)
library("GenABEL")

# B. QQ plot
thisqq <- qq(consoritum_metal_output$`P-value`)
# Example of meta-analysis output file = "meta_analysis_courage_05052021.tbl" (See 04. Linear regression.R)

ggsave(plot = thisqq, filename = "qq.png", units = c("in"), dpi = 600, height = 5.5, width = 5.5)

# C. Estimation of lambda
estlambda(consoritum_metal_output$`P-value`, method="median")
