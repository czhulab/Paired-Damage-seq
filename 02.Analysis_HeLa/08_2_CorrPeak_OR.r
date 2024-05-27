library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
theme_set(theme_cowplot())

GetFisherTable <- function(filePath) {
  # Read the first 4 lines from the file
  lines <- readLines(filePath, n = 4)
  
  # Extract the digits from these lines
  numbers <- sapply(lines, function(line) {
    as.numeric(gsub(".*: ", "", line))
  })
  
  # Assign the extracted numbers to variables for clarity
  numQueryIntervals <- numbers[1]
  numDbIntervals <- numbers[2]
  numOverlaps <- numbers[3]
  numPossibleIntervals <- numbers[4]
  
  # Calculate the remaining cells of the contingency table
  in_a_not_in_b <- numQueryIntervals - numOverlaps
  not_in_a_in_b <- numDbIntervals - numOverlaps
  not_in_a_not_in_b <- numPossibleIntervals - numQueryIntervals - numDbIntervals + numOverlaps
  
  # Construct the contingency table
  contingencyTable <- matrix(c(numOverlaps, in_a_not_in_b, not_in_a_in_b, not_in_a_not_in_b),
                             nrow = 2, byrow = TRUE,
                             dimnames = list(c("in -a", "not in -a"),
                                             c("in -b", "not in -b")))
  
  return(contingencyTable)
}



GetFisherResult <- function(fisher_out){
    p_val <- fisher_out$p.value
    lower <- fisher_out$conf.int[1]
    upper <- fisher_out$conf.int[2]
    OR <- as.numeric(fisher_out$estimate)
    return(c(p_val, lower, upper, OR))
}



Comparisons <- c('enhancer', 'H3K27me3', 'H3K27ac', "H3K9me3", "H3K4me1", "H3K36me3", "H3K4me3")

###################################################################################################
###################################################################################################
########################################### Union Peak ############################################
UnionPeak_out <- list()
UnionPeakResultdir <- "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/TPA/DS/AllPeak/"
UnionFiles <- list.files(UnionPeakResultdir, "UnionPeak_", full.names = TRUE)

for (h in Comparisons){
    file <- grep(h, UnionFiles, value = TRUE)
    table <- GetFisherTable(file)
    results <- fisher.test(table) %>% GetFisherResult()
    UnionPeak_out[[length(UnionPeak_out) + 1]] <- results
}

UnionResults <- do.call(rbind, UnionPeak_out)
rownames(UnionResults) <- Comparisons
colnames(UnionResults) <- c("p_val", "CI_low", "CI_high", "OR")
UnionResults <- as.data.frame(UnionResults)

###################################################################################################
###################################################################################################
######################################### Top     Peak ############################################
TopPeak_out <- list()
TopPeakResultdir <- "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/TPA/DS/FDR/Positive/"
TopPeakFiles <- list.files(TopPeakResultdir, full.names = TRUE)

for (h in Comparisons){
    file <- grep(h, TopPeakFiles, value = TRUE)
    table <- GetFisherTable(file)
    results <- fisher.test(table) %>% GetFisherResult()
    TopPeak_out[[length(TopPeak_out) + 1]] <- results
}

TopPeakResults <- do.call(rbind, TopPeak_out)
rownames(TopPeakResults) <- Comparisons
colnames(TopPeakResults) <- c("p_val", "CI_low", "CI_high", "OR")
TopPeakResults <- as.data.frame(TopPeakResults)

###################################################################################################
###################################################################################################
##################################### bottom      Peak ############################################
BotPeak_out <- list()
BotPeakResultdir <- "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/TPA/DS/FDR/Negative/"
BotPeakFiles <- list.files(BotPeakResultdir, full.names = TRUE)

for (h in Comparisons){
    file <- grep(h, BotPeakFiles, value = TRUE)
    table <- GetFisherTable(file)
    results <- fisher.test(table) %>% GetFisherResult()
    BotPeak_out[[length(BotPeak_out) + 1]] <- results
}

BotPeakResults <- do.call(rbind, BotPeak_out)
rownames(BotPeakResults) <- Comparisons
colnames(BotPeakResults) <- c("p_val", "CI_low", "CI_high", "OR")
BotPeakResults <- as.data.frame(BotPeakResults)
###################################################################################################
###################################################################################################
##################################### data transformation #########################################
TopPeakResults$Comparison <- rownames(TopPeakResults)
TopPeakResults$Dataset <- "Positive Linked Peaks"

UnionResults$Comparison <- rownames(UnionResults)
UnionResults$Dataset <- "Union Peak"

BotPeakResults$Comparison <- rownames(BotPeakResults)
BotPeakResults$Dataset <- "Negative Linked Peaks"
# Combine the two datasets
CombinedResults <- rbind(TopPeakResults, UnionResults, BotPeakResults) 

# Convert to long format for ggplot
CombinedResults$Comparison <- factor(CombinedResults$Comparison, levels = unique(CombinedResults$Comparison))

###################################################################################################
###################################################################################################
############################################ Plotting #############################################
# Define colors
dotCOLS = c("#a6d8f0","#f9b282", "#54e87d")
barCOLS = c("#008fd5","#de6b35", "#22ab48")


# Plot
p <- ggplot(CombinedResults, aes(x = Comparison, y = OR, color = Dataset, fill = Dataset)) + 
  geom_linerange(aes(ymin = CI_low, ymax = CI_high), size=1, position=position_dodge(width = 0.5)) + # Corrected line
  #geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  scale_fill_manual(values = barCOLS) +
  geom_hline(yintercept=1, lty=2) +
  scale_color_manual(values = dotCOLS) +
  coord_flip() +
  labs(x = "", y = "Odds Ratio and Confidence Interval", fill = "Dataset", color = "Dataset")

# Display the plot
pdf("/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/Joint/TPA/DS/Plots/FDR_forest_DS.pdf")
p
dev.off()








