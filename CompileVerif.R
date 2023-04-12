library("ctv")
ctv2html("Distributions.md")
browseURL("Distributions.html")

#see help at https://github.com/cran-task-views/ctv/blob/main/Contributing.md

check_ctv_packages("Distributions.md")




remotes::install_github("DylanDijk/CTVsuggest")
library("CTVsuggest")
tvsugg2 <- CTVsuggest(taskview = "Distributions", n = 10)
?CTVsuggest


library(RWsearch)
crandb_down()
tvdb_down()

intv <- tvdb_pkgs("Distributions")
tvsugg <- setdiff(s_crandb("distributions", "probability", mode = "and"), intv)

tvsugg2$Packages %in% tvsugg
