# clear the environment
rm(list = ls())

# Set up the working directory where you put your results
setwd("/Users/emilyluytan/Desktop/PbPbC_Analysis/Exposure_Analysis_N2Z_MediaSeperated/Chem_347_M2/Windows_Analysis/Pb_Min_Pd_Accum")

## all csv files in a specified path
file_list <- list.files(path = "/Users/emilyluytan/Desktop/PbPbC_Analysis/Exposure_Analysis_N2Z_MediaSeperated/Chem_347_M2/Windows_Analysis/Pb_Min_Pd_Accum", pattern = ".csv")

for (i in file_list) {
  temp_data <- read.csv(i, header = TRUE)
  temp_data_new <- temp_data[,2:ncol(temp_data)]
  colnames(temp_data_new)[1] = "disease"
  
  temp_regression <- glm(disease ~ ., family = "binomial", data = temp_data_new)
  sum_temp_regression <- summary(temp_regression)
  
  df_basic <- data.frame(varible = colnames(temp_data)[2], num_sample = nrow(temp_data),
                         call = toString(sum_temp_regression$call), 
                         aic = sum_temp_regression$aic,
                         iteration = sum_temp_regression$iter)
  write.table(df_basic, file = paste(substr(i,1,nchar(i)-4),"_results.csv", sep = ""), sep = ",", col.names = TRUE, row.names = FALSE, append = TRUE)
  
  temp_regression_coefficient <- data.frame(sum_temp_regression$coefficients)
  temp_regression_coefficient$oddsRatio = exp(temp_regression_coefficient$Estimate)
  temp_regression_coefficient$CI_oddsratio_low = exp(temp_regression_coefficient$Estimate - 1.96 * temp_regression_coefficient$Std..Error)
  temp_regression_coefficient$CI_oddsratio_up = exp(temp_regression_coefficient$Estimate + 1.96 * temp_regression_coefficient$Std..Error)

  temp_regression_coefficient$CI_oddsratio_low_2 = exp(temp_regression_coefficient$Estimate - 2.58 * temp_regression_coefficient$Std..Error)
  temp_regression_coefficient$CI_oddsratio_up_2 = exp(temp_regression_coefficient$Estimate + 2.58 * temp_regression_coefficient$Std..Error)
  
  write.table(temp_regression_coefficient,file = paste(substr(i,1,nchar(i)-4),"_results.csv", sep = ""), sep = ",", col.names = TRUE, row.names = TRUE, append = TRUE)
}