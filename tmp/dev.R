
devtools::build_vignettes()

devtools::document()

devtools::test()

devtools::load_all()

devtools::check()

devtools::build()

devtools::install()

browseVignettes("QCtreeCNV")

library(data.table)

f <-  QCtreeCNV:::step3
debugonce(f)

rhub::check_on_linux("../")
rhub::check_on_windows("../")
rhub::check_for_cran("../")
