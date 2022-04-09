library(stringr)

X_ref <- as.matrix(read.csv("data/kidney_expression.csv"))
kidney_sample_data <- read.csv("data/kidney_sample_data.csv")
kidney_sample_data$external_id <- str_replace_all(kidney_sample_data$external_id,"-",".")

useful_map_cnts <- colSums(X_ref)

cnt.df <- data.frame(external_id = names(useful_map_cnts),useful_map_cnts=useful_map_cnts)
kidney_sample_data <- merge(kidney_sample_data,cnt.df)

save(kidney_sample_data,file="./data/kidney_sample_data.rda")