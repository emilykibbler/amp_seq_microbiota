print_norm_check <- function(input_test){
  if (unlist(input_test)[2] < 0.05) {
  print("P value indicates not normally distributed data")
  }
  if (unlist(input_test)[2] >= 0.05) {
    print("P value indicates normally distributed data")
  }
}

normal_stats <- function(ps) {
  # Error checks
  if (class(ps) != "phyloseq") {
    print("Second argument should be a phyloseq object")
    print("Probably one that ends in something like _clean_rar")
  }
  print("Reminder:")
  print("I made this function to evaluate Chao1 and Shannon metrics")
  df <- estimate_richness(ps, measures = c("Chao1", "Shannon")) 
  # # measure evenness for each sample
  even <- df$Shannon/log(df$Chao1)
  # sd is sample data here not standard deviation
  rar_sd = as(sample_data(ps), "matrix")
  rar_sd = as.data.frame(rar_sd)
  result_df <- cbind(df, even, rar_sd)

  # # make a histogram to look at the shape of the data (bell curve? skew?). You can save this graph for your own benefit if you want.
  hist(result_df$Chao1)
  # Kurtosis
  print("Kurtosis value is:")
  print(kurtosis(result_df$Chao1))
  print("Red flag for kurtosis is absolute value >2; orange flag >1")
  print("Positive or negative value indicates direction of tail")
  # Shapiro - Shannon
  shan_shap <- shapiro.test(result_df$Shannon)
  print("Shannon diversity shapiro test:")
  print(shan_shap)
  print_norm_check(shan_shap)
  # Shapiro - Observed
  chao_shap <- shapiro.test(result_df$Chao1)
  print("Chao1 diversity shapiro test:")
  print(chao_shap)
  print_norm_check(chao_shap)
  # Shapiro - Evennes
  even_shap <- shapiro.test(result_df$even)
  print("Evenness shapiro test:")
  print(even_shap)
  print_norm_check(even_shap)
 return(result_df)
}


