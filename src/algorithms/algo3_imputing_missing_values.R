
# Function to impute missing values
imputed.df.fun <- function(dataset) {
  # Imputed df in logit space
  dataset <- dataset[!rowSums(dataset > 1, na.rm = TRUE),]
  imputed.df.logit <- logit(dataset) %>% impute_EM() 
  imputed.df <- logitInv(imputed.df.logit)
  sum.imp <- c()
  for (i in 1:nrow(dataset)) {
    # For each row, sum the missing elements which are imputed and subtract them from U
    sum.imp[i] <- sum(imputed.df[i, which(is.na(dataset[i, ]))])
  }
  # Yields a set of realistic compositions, from compositions with missing values 
  imputed.df$U <- abs(imputed.df$U - sum.imp)
  return(imputed.df)
}

gen.mval.fun <- function(df, prop) {
  # inner function, select random elements and take their sum
  random_row_sum <- function(row) {
    # select 2 values
    selected_indices <- sample(length(row), floor(prop * length(row)))
    selected_values <- row[selected_indices]
    sum(selected_values)
    return(list(sum(selected_values), selected_indices))
  }
  
  # create a sum vector 
  sum.vec <- c()
  for (i in 1:nrow(df)) {
    random.row.o <- random_row_sum(df[i, ])
    sum.vec[i] <- random.row.o[[1]]
    df[i, c(random.row.o[[2]])] <- NA
  }
  output <- cbind(df, sum.vec) %>% as_tibble()
  return(output)
}
