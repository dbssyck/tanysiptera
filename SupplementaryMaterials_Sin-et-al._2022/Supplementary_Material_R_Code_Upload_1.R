# Supplementary Material: Tanysiptera Bioacoustics
#
#
#
# Load Libraries ####
library(ggfortify)
library(MKinfer)
library(psych)
library(tanysiptera)


# Import file ####
dataset <- as.data.frame (read.csv ("\\path\\to\\file.csv", sep = ""), header = TRUE)


# Clean up dataset ####
## Aggregate parameters from same sound recordings
dataset_file_aggregated <- aggregate (dataset [, colnames(dataset) != "Accession.Number" ], list (dataset$Accession.Number), mean)
colnames (dataset_file_aggregated) [1] <- "Accession.Number" # Change first column name from Group.1 to "Accession.Number"

## Combine recordings from same sound individual
### Retrieve information on repeated individuals because they were lost during aggregate(). (Probably not the most elegant method but it works.)
levels(dataset_file_aggregated$Same.Individual) <- levels(dataset$Same.Individual)
for (i in 1:nrow(dataset_file_aggregated)) {
  rownumber <- which(dataset$Accession.Number == dataset_file_aggregated$Accession.Number[i])
  dataset_file_aggregated$Same.Individual[i] <- as.character(dataset$Same.Individual[rownumber[1]])
}

### Aggregate parameters across sound recordings of same individuals
dataset_indiv_aggregated <- aggregate (dataset_file_aggregated, list (dataset_file_aggregated$Same.Individual), mean, na.rm = T)
colnames (dataset_indiv_aggregated) [1] <- "Same.Individual" # Change first column name from Group.1 to "Same.Individual"

### Retrieve information on Accession (pick one representative) because they were lost during aggregate(). (Probably not the most elegant method but it works.)
for (i in 1:nrow(dataset_indiv_aggregated)) {
  rownumber <- which(dataset$Same.Individual == dataset_indiv_aggregated$Same.Individual[i])
  dataset_indiv_aggregated$Accession.Number[i] <- as.character(dataset$Accession.Number[rownumber[1]])
}

### Create dataset for further analysis by extracting relevant columns and recombining dataset with individuals that were not across multiple recordings
relevant_columns <- c("Accession.Number", "Species", "Taxa", "Centre.Frequency..Hz.", "Bandwidth.of.first.note..Hz.", "Duration.of.first.note..seconds.", "Shape.of.second.note", "Bandwidth.ratio.of.second.note.to.last.note","Duration.ratio.of.second.note.to.last.note", "Total.number.of.notes","Position.of.note.with.lowest.centre.frequency")
dataset_analyse <- rbind ( dataset_indiv_aggregated [, colnames(dataset_indiv_aggregated) %in% relevant_columns],
                           dataset_file_aggregated [ is.na (dataset_file_aggregated$Same.Individual), colnames(dataset_file_aggregated) %in% relevant_columns] )
for (i in 1:nrow(dataset_analyse)) {
  rownumber <- which(dataset$Accession.Number == dataset_analyse$Accession.Number[i])
  dataset_analyse$Species[i] <- as.character(dataset$Species[rownumber[1]])
  dataset_analyse$Taxa[i] <- as.character(dataset$Taxa[rownumber[1]])
} ### Retrieve information on species and taxa. Again might not be the most elegant method (?) but it works.


# ### To remove relevant taxa and rerun PCA
# taxa_remove <- c("hydrocharis",
#                  # "rosseliana",
#                  "obiensis")
#
# dataset_analyse <- dataset_analyse[ !(dataset_analyse$Taxa %in% taxa_remove), ]

vocal_parameters <- dataset_analyse[,-(1:3)]



# Check correlation ####
pearson <- cor(vocal_parameters, method = "pearson")



# PCA ####
## Check suitability using Bartlett's Test
cortest.bartlett( cor(vocal_parameters), n = 75)


## Run PCA
pca_scores <- prcomp (vocal_parameters, scale = T)
pca_summary <- summary(pca_scores)
pca_summary$sdev ^ 2 # obtain eigenvalue
pca_scores$rotation # get loadings


## Scree plot (Adapted from https://www.analyticsvidhya.com/blog/2016/03/pca-practical-guide-principal-component-analysis-python/)
pr_var <- pca_scores$sdev^2 # calculate variance
prop_varex <- pr_var/sum(pr_var) # proportion of variance explained
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b") # cumulative screen plot (same as above, but inverted)


## Combine PCA scores and relevant information
post_pca <- data.frame(pca_scores$x, Taxa =  paste (dataset_analyse$Taxa), stringsAsFactors=FALSE)


## Plot to visualise
autoplot(pca_scores, data = post_pca, colour = 'Taxa')
autoplot(pca_scores, data = post_pca, colour = 'Taxa', x = 2, y = 3)


# Add cluster information ####
post_pca$cluster = NA
post_pca[grep("carolinae", post_pca[,"Taxa"]), "cluster"]  <- "carolinae"
post_pca[grep("ellioti", post_pca[,"Taxa"]), "cluster"]  <- "ellioti"
post_pca[grep("riedelii", post_pca[,"Taxa"]), "cluster"]  <- "riedelii"
post_pca[grep("hydrocharis", post_pca[,"Taxa"]), "cluster"]  <- "hydrocharis"
post_pca[grep("meyeri|galatea|minor|vulcani", post_pca[,"Taxa"]), "cluster"]  <- "galatea"
post_pca[grep("acis|boanensis|nais", post_pca[,"Taxa"]), "cluster"]  <- "nais"
post_pca[grep("obiensis", post_pca[,"Taxa"]), "cluster"]  <- "obiensis"
post_pca[grep("rosseliana", post_pca[,"Taxa"]), "cluster"]  <- "rosseliana"
post_pca[grep("emiliae", post_pca[,"Taxa"]), "cluster"]  <- "doris"
post_pca[grep("doris", post_pca[,"Taxa"]), "cluster"]  <- "doris"
post_pca[grep("brunhildae|sabrina|margarethae|browningi", post_pca[,"Taxa"]), "cluster"]  <- "margarethae"



# ANOVA ####
## Remove clusters with sample size < 3
run_anova <- post_pca [ post_pca$cluster %in% names (subset (table(post_pca$cluster), table(post_pca$cluster) >2)), ]
# write.csv(ggdata_cluster, paste ("C:\\Users\\Keita\\Documents\\ZZZ_Keita_Files\\University\\D_Common-Paradise-kingfisher\\5_Results\\ggdata_cluster.csv", sep = ""), row.names=TRUE)

## Conduct one way anova and Tukey test for each PC
one_way_anova_PC1 <- aov(formula = PC1 ~ cluster, data = run_anova)
tukey_PC1 <- TukeyHSD(one_way_anova_PC1)
one_way_anova_PC2 <- aov(formula = PC2 ~ cluster, data = run_anova)
tukey_PC2 <- TukeyHSD(one_way_anova_PC2)
one_way_anova_PC3 <- aov(formula = PC3 ~ cluster, data = run_anova)
tukey_PC3 <- TukeyHSD(one_way_anova_PC3)




# Diagnosability using nonparametric bootstrap ####
sample_size <- table(dataset_analyse$Taxa) # Generate sample size of each taxa
list_of_taxa_to_compare <- names (sample_size [sample_size > 2]) # Generate list of taxa to compare
combinations <- t (combn(list_of_taxa_to_compare, 2)) # generate list of pairwise combinations

run_bootstrap <- dataset_analyse [, c("Duration.ratio.of.second.note.to.last.note", "Bandwidth.ratio.of.second.note.to.last.note","Taxa") ] # generate dataframe of vocal parameters to run bootstrap analysis

DBM_pairwise <- matrix (nrow = nrow(combinations), ncol = (ncol(run_bootstrap) - 1) ) # create empty matrix to store diagnosability output
rownames(DBM_pairwise) <- paste(combinations[,1], "vs", combinations[,2])
colnames(DBM_pairwise) <- colnames(run_bootstrap)[-3]

## Run nonparametric diagnosability test for all pairwise comparisons. (Probably not the most elegant method but it works.)
for (pair in 1:nrow(DBM_pairwise)) {
  for (param in 1:ncol(DBM_pairwise)) {
    taxa1 <- combinations[pair, 1]
    taxa2 <- combinations[pair, 2]
    param <- colnames(DBM_pairwise)[param]

    values_taxa1 <- run_bootstrap [ run_bootstrap$Taxa == taxa1, param ]
    values_taxa2 <- run_bootstrap [ run_bootstrap$Taxa == taxa2, param ]

    bootstrap <- boot.t.test(values_taxa1, values_taxa2,
                             alternative = c("two.sided"),
                             paired = FALSE,
                             var.equal = FALSE,
                             conf.level = 0.975,
                             R = 10000)

    if ( bootstrap$p.value < 0.025 )
    { DBM_pairwise[pair,param] <- 1 } # if diagnosable, print 1
  }
}



# Isler Criterion
vocal_diagnose(dataset_analyse[,-c(1,2)])


