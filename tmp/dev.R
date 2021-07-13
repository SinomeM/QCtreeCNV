
devtools::build_vignettes()

devtools::document()

devtools::test()

devtools::load_all()

devtools::check()

devtools::build()

library(data.table)

f <-  QCtreeCNV:::step3
debugonce(f)
