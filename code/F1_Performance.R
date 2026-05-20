#test the performance of your model on synthetic generated test spots you can use this function to benchmark and get a sense of the model's performance.
#'
#' @param pred Object of class matrix containing the ground truth composition of each spot as obtained from the function syn_spot_comb_topic_fun.R, 2nd element.
#' @param real Object of class matrix with the predicted topic probability distributions for each spot.
#' @return This function returns a list with TP, TN, FP, FN and the Jensen-Shannon Divergence index.
#' @export
#' @examples
#'

test_synthetic_perf_ns <- function(pred,
                                       real) {
  # Check variables
  if (!is.matrix(pred)) stop("ERROR: pred must be a matrix object!")
  if (!is.matrix(real)) stop("ERROR: real must be a matrix object!")
  
  colnames(real) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                          x = colnames(real),
                                          perl = TRUE)
  colnames(pred) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                             x = colnames(pred),
                                             perl = TRUE)
  #load required packages
  suppressMessages(require(philentropy))
  
  ##### Get metrics between real-predicted proportions #####
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  for (i in seq_len(nrow(pred))) {
    
    # Create matrix to feed to metrics
    x <- rbind(pred[i, ],
               real[i, ])
    cat(sprintf("i: %s", i), sep = "\n")
    #### Calculate TP-TN-FP-FN ####
    for (index in colnames(pred)) {
      if (x[1, index] > 0 & x[2, index] > 0) {
        tp <- tp + 1
      } else if (x[1, index] == 0 & x[2, index] == 0) {
        tn <- tn + 1
      } else if (x[1, index] > 0 & x[2, index] == 0) {
        fp <- fp + 1
      } else if (x[1, index] == 0 & x[2, index] > 0) {
        fn <- fn + 1
      }
    }; rm(index)
    
  }; rm(i)
  
  #### Performance metrics ####
  accuracy <- round((tp + tn) / (tp + tn + fp + fn), 2)
  sensitivity <- round(tp / (tp + fn), 2)
  specificity <- round(tn / (tn + fp), 2)
  precision <- round(tp / (tp + fp), 2)
  recall <- round(tp / (tp + fn), 2)
  FPR <- round(fp / (tn + fn), 2)
  F1 <- round(2 * ((precision * recall) / (precision + recall)), 2)
  F0.5 <- round(1.25 * ((precision * recall) / ((0.25 * precision) + recall)), 2)
  
  
  cat(sprintf("The following summary statistics are obtained:
              Accuracy: %s,
              Sensitivity: %s,
              Specificity: %s,
              Precision(TPR): %s,
              recall: %s,
              F1 score: %s,
              F0.5 score: %s
              TP: %s, 
              TN: %s, 
              FP: %s, 
              FN: %s,
              FPR: %s",
              accuracy, sensitivity, specificity, precision, recall, F1, F0.5, tp, tn, fp, fn, FPR), sep = "\n")
  
  cat("raw statistics are returned in the list - TP, TN, FP, FN",
      sep = "\n")
  return(list(TP = tp, TN = tn, FP = fp, FN = fn, Accuracy = accuracy, 
              Sensitivity = sensitivity, Specificity = specificity, 
              Precision = precision, Recall = recall, F1_score = F1, F0.5_score = F0.5))
}
