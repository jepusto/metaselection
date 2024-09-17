
dat <- read.csv2("data-raw/full_tidy.csv", stringsAsFactors = FALSE) 
mods <- names(dat)[14:31]
self_control <- dat[, c("studyid", "esid", "name", "g", "var.g", "se.g", "outcome", "comparison", mods)]
names(self_control) <- gsub(x = names(self_control), pattern = "\\.", replacement = "_")  


usethis::use_data(self_control, overwrite = TRUE)
