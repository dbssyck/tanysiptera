#' vocal_diagnose
#'
#' This function examines vocal diagnosability based on the criterion designed by Isler, Isler & Whitney (1998).
#' \cr\cr Reference: Isler, M. L., Isler, P. R., & Whitney, B. M. (1998). Use of vocalizations to establish species limits in antbirds (Passeriformes: Thamnophilidae). The Auk, 115(3), 577-590.
#'
#' To be diagnosable, vocal parameters must meet both of the following criteria:
#' 1. No overlap of measurements between taxa
#' 2. x_a + t_a * SD_a < x_b - t_b * SD_b, where x = mean, SD = standard deviation, t = t-score at 97.5 percentile of the t distribution for n-1 degrees of freedom, and population a is the population with smaller set of measurements, population b with bigger set of measurements.
#'
#'
#'@param data Dataframe containing vocal parameters. Each row should include parameters from one vocal sample (i.e., one individual). Column 1 must contain taxa/population assignment of the sample, and subsequent columns must contain measurements of one vocal parameter each.
#'@param min.sample Default = 2. Minimum sample size allowed for each taxa/populations for Isler test to be conducted.
#'@param prose Default = FALSE. If TRUE, the result will be a dataframe with rownames spelled out as "Taxa 1 vs Taxa 2" for each pairwise comparison. If FALSE, the result will be a dataframe with Taxa 1 in first column and Taxa 2 in second column for each pairwise comparison.
#'
#'@examples
#'test_data <- read.csv(system.file("extdata", "test_data.txt", package="tanysiptera"), sep = "\t")
#'head(test_data)
#'vocal_diagnose(test_data)
#' @export


# Initiate function
vocal_diagnose <- function(data, min.sample = 2, prose = FALSE) {

# Store dataset and change name of column 1
dataset <- data
colnames(dataset)[1] <- "taxa"

# Generate sample size for each taxa
sample_size <- table(dataset$taxa)

# Remove taxa with sample size of 1, or smaller than specified min.sample value. If minimum sample size specified is 1, overwrite and make it 2.
if( min.sample==1 ) { minimum_sample_size <-2 } else { minimum_sample_size <- min.sample }
samples_checked <- sample_size [ sample_size >=  minimum_sample_size ] # generate list of taxa above minimum sample size
dataset_trimmed <-  dataset [ dataset$taxa %in% names(samples_checked), ] # subset dataset to those above minimum sample size
vocal_parameters <- subset(dataset_trimmed, select=-c(taxa)) # create dataframe with numerical values (vocal parameters only) for subsequent calculations

# Generate list of pairwise combinations
combinations <- t (combn(names (samples_checked), 2))

# Create function to test for overlap
overlap <- function (taxa1, taxa2, param, values) {
  ## Extract respective means so that population a and population b (smaller and bigger set of measurements respectively) can be determined
  vocal_means_taxa1 <- vocal_means [ vocal_means$Group.1 == taxa1, colnames(vocal_means) == param ]
  vocal_means_taxa2 <- vocal_means [ vocal_means$Group.1 == taxa2, colnames(vocal_means) == param ]
  col_num <- which (colnames(values$plus) == param) # extract column number of relevant parameter (same across values$plus and values$minus)
  row_names <- rownames(values$plus) # generate vector of rownames so that relevant rows can be easily extracted in the function

  if ( ( vocal_means_taxa1 < vocal_means_taxa2 ) &
       ( values$plus [ which (row_names == taxa1), col_num ] > values$minus [ which (row_names == taxa2), col_num ] ) )
  { 0 } else {
    if ( ( vocal_means_taxa1 > vocal_means_taxa2 ) &
         ( values$minus [ which (row_names == taxa1), col_num ] < values$plus [ which (row_names == taxa2), col_num ] ) )
    { 0 } else { 1 } # if no overlap, print 1
  }
}

# Aggregated values of each vocal parameter for each taxa
vocal_means <- aggregate ( vocal_parameters, list (dataset_trimmed$taxa), mean) # means
vocal_sds <- aggregate (vocal_parameters, list (dataset_trimmed$taxa), sd) # standard deviations
vocal_max <- aggregate ( vocal_parameters, list (dataset_trimmed$taxa), max) # maximum values
vocal_min <- aggregate ( vocal_parameters, list (dataset_trimmed$taxa), min) # minimum values
rownames(vocal_max) <- vocal_max$Group.1
rownames(vocal_min) <- vocal_min$Group.1
raw_values <- list ( "plus" = vocal_max, "minus" = vocal_min) # store max and min values into list for overlap() function

# Student's t-score at 97.5th percentile, for N-1 degrees of freedom where N is the sample size, for each taxa
vocal_tscores <- qt(0.975, df = samples_checked - 1)

# Create empty matrix to store the values required to check if criterion 2 is met
criteria_two_check_plus <- matrix (nrow = length(samples_checked), ncol = ncol(vocal_parameters))
rownames(criteria_two_check_plus) <- names (samples_checked)
colnames(criteria_two_check_plus) <- colnames (vocal_parameters)
criteria_two_check_minus <- criteria_two_check_plus

# Calcualte values for criteria two check
for (taxa in 1:nrow(criteria_two_check_plus)) {
  for (parameter in 1:ncol(criteria_two_check_plus)) {
    taxa_name <- rownames (criteria_two_check_plus) [taxa]
    param_name <- colnames (criteria_two_check_plus) [parameter]

    x <- vocal_means [ vocal_means$Group.1 == taxa_name, colnames (vocal_means) == param_name ]
    t <- vocal_tscores [ names (vocal_tscores) == taxa_name ]
    SD <- vocal_sds [ vocal_sds$Group.1 == taxa_name, colnames (vocal_sds) == param_name ]
    criteria_two_check_plus [taxa, parameter] <- x + t * SD
    criteria_two_check_minus [taxa, parameter] <- x - t * SD
  }
}

formula_values <- list ( "plus" = criteria_two_check_plus, "minus" = criteria_two_check_minus) # store max and min values into list for overlap() function

# Create empty matrix to output diagnosability and inspection matrix
diagnosability <- matrix (nrow = nrow(combinations), ncol = ncol(vocal_parameters))
colnames(diagnosability) <- colnames (vocal_parameters)
inspect <- diagnosability

# Test whether both criteria are met for all pairwsie comparisons
for (pair in 1:nrow(diagnosability)) {
  for (parameter in 1:ncol(diagnosability)) {
    taxa1 <- combinations[pair, 1]
    taxa2 <- combinations[pair, 2]
    param_name <- colnames(diagnosability)[parameter]
    # If criteria 1 and 2 are both diagnosable, print 1 in the diagnosability matrix and 3 in the inspect matrix
    if ( ( overlap(taxa1, taxa2, param_name, raw_values) == 1 ) & ( overlap (taxa1, taxa2, param_name, formula_values) == 1 ) )  {
      diagnosability [ pair, parameter ] <- 1
      inspect [ pair, parameter ] <- 3 } else {
        # If criteria 1 only, print 1 in the inspect matrix
        if ( ( overlap (taxa1, taxa2, param_name, raw_values) == 1 ) & ( overlap (taxa1, taxa2, param_name, formula_values) != 1 ) )
        { inspect [ pair, parameter ] <- 1 } else {
          # If criteria 2 only, print 2 in the inspect matrix
          if ( ( overlap (taxa1, taxa2, param_name, raw_values) != 1 ) & ( overlap (taxa1, taxa2, param_name, formula_values) == 1 ) )
            inspect [ pair, parameter ] <- 2 }
      }
  }
}

# Sum columns to count number of pairs diagnosed per parameter
param.sum <- colSums (diagnosability, na.rm = T)

# Sum rows to count number of diagnosable parameters per pair
diagnosability <- cbind (diagnosability, rowSums (diagnosability, na.rm = T))
colnames(diagnosability)[ncol(diagnosability)] <- "No. Diagnosable Parameters"

# Fix rownames accordingly
if( prose == FALSE ) {
  diagnosability <- cbind ( "Taxa 1" = combinations[,1], "Taxa 2" = combinations[,2], data.frame (diagnosability) )
  inspect <- cbind ( "Taxa 1" = combinations[,1], "Taxa 2" = combinations[,2], data.frame (inspect) )
} else {
  rownames(diagnosability) <- paste(combinations[,1], "vs", combinations[,2])
  rownames(inspect) <- paste(combinations[,1], "vs", combinations[,2])
}

# Store outputs into list
diagnose <<- list( "diagnosability" = data.frame (diagnosability),
                   "inspect" = data.frame(inspect),
                   "param.sum" = data.frame(param.sum))
} # End of function
