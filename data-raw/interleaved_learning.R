

dat <- read.csv("data-raw/Brunmair_data.csv", stringsAsFactors = FALSE) 

names(dat) <- tolower(names(dat))

mods <- names(dat)[c(3, 10:22)]
interleaved_learning <- dat[, c("name_year", "effectsizeid", "g", "vg", "studyid", "article_num", "items_num", mods)]
names(interleaved_learning)[1] <- "study" 
names(interleaved_learning)[2] <- "esid" 
names(interleaved_learning)[5] <- "sample_id" 
names(interleaved_learning)[7] <- "item_num" 

usethis::use_data(interleaved_learning, overwrite = TRUE)
