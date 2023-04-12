library(tidyverse)
library(sqldf)
source("UtilityFunctions.R")

# Read in the initial Yeast Sequence Data Set from UCI, and give it functional names
df1 <- read.csv("yeast.csv", header = FALSE, sep = "", dec = ".")
colnames(df1) <- c("Sequence", "MCG", "GVH", "ALM", "MIT", "ERL", "POX", "VAC", "NUC", "Class")

# 1
# clean the given data set of all extraneous data that doesn't represent a full DNA sequence
# Then, remove all duplicates
# This shaves the total size of the data set down nearly 200 instances
df1_clean <- subset(df1, Class != "EXC" & Class != "VAC" & Class != "POX" & Class != "ERL")
df1_clean <- unique(df1_clean)

# 2
# make two subsets for data that represent the two subgroups of test organelles
# These two subsets will be used in the future manipulations
df2 <- subset(df1_clean, Class == "MIT" | Class == "NUC")
df3 <- subset(df1_clean, Class == "ME1" | Class == "ME2" | Class == "ME3" | Class == "CYT")

# 3
# subset df2 is partitioned into a subset of mitochondrial sequences
# The new subset has a couple attributes measured and displayed placed in a new data frame
df2_mit <- subset(df2, Class == "MIT")
df2_mit_dist <- data.frame(df2_mit$Sequence, df2_mit$MCG, df2_mit$GVH)
colnames(df2_mit_dist) <- c("MIT Sequence", "MCG", "GVH")

# 4
# the new data frame df2_mit_dist is operated on
# Then, the results of the operations are compared against the grandparent set (df2)
mcg_min <- min(df2_mit_dist$MCG)
gvh_min <- min(df2_mit_dist$GVH)
mcgtot_min <- min(df2$MCG)
gvhtot_min <- min(df2$GVH)
mins_names <- c("Min MCG (MIT)", "Min MCG (Total)", "Min GVH (MIT)", "Min GVH (Total)")
dfmins <- data.frame(mcg_min, mcgtot_min, gvh_min, gvhtot_min)
colnames(dfmins) <- mins_names
mcg_max <- max(df2_mit_dist$MCG)
gvh_max <- max(df2_mit_dist$GVH)
mcgtot_max <- max(df2$MCG)
gvhtot_max <- max(df2$GVH)
maxes_names <- c("Max MCG (MIT)", "Max MCG (Total)", "Max GVH (MIT)", "Max GVH (Total)")
dfmaxes <- data.frame(mcg_max, mcgtot_max, gvh_max, gvhtot_max)
colnames(dfmaxes) <- maxes_names
mcg_avg <- trunc(mean(df2_mit_dist$MCG)*10^2)/10^2
gvh_avg <- trunc(mean(df2_mit_dist$GVH)*10^2)/10^2
mcgtot_avg <- trunc(mean(df2$MCG)*10^2)/10^2
gvhtot_avg <- trunc(mean(df2$GVH)*10^2)/10^2
avgs_names <- c("Avg MCG (MIT)", "Avg MCG (Total)", "Avg GVH (MIT)", "Avg GVH (Total)")
dfavgs <- data.frame(mcg_avg, mcgtot_avg, gvh_avg, gvhtot_avg)
colnames(dfavgs) <- avgs_names
df2_mit_dist_final <- cbind(dfmins, dfmaxes, dfavgs)

# 5
# we return to normal df2 and scrub out all Yeast sample names with a B in them
# From the resultant set, we bin the data according to the ALM score, in 6 bins
df2_SEQ <- df2[ , "Sequence"]
df2_temp <- grep("B", df2_SEQ, value = TRUE)
df2_SEQ_clean <- data.frame(df2_temp)
df2_temp_clean <- sqldf("select Sequence from 'df2' except select * from 'df2_SEQ_clean'")
df2_noB <- merge(df2, df2_temp_clean, by = "Sequence")
df2_noB_bins <- cut(df2_noB$ALM, breaks = 6, labels = c("First", "Second", "Third", "Fourth", "Fifth", "Sixth"))
bins_table <- table(df2_noB_bins)

# 6
# the df2 subset is partitioned by its NUC score (not the NUC class value, but the nuclear protein distribution ratio)
# The new subset is then sampled 3 times and placed into a new data frame
df2_nuc <- subset(df2, NUC == 0.22)
df2_nuc_s1 <- data.frame(df2_nuc[sample(nrow(df2_nuc), size = 50, replace = TRUE), ])
df2_nuc_s2 <- data.frame(df2_nuc[sample(nrow(df2_nuc), size = 50, replace = TRUE), ])
df2_nuc_s3 <- data.frame(df2_nuc[sample(nrow(df2_nuc), size = 50, replace = TRUE), ])
df2_nuc_temp <- merge(df2_nuc_s1, df2_nuc_s2, by = "Sequence")
df2_nuc_samples <- merge(df2_nuc_temp, df2_nuc_s3, by = "Sequence")

# 7
# the data has the duplicate attributes removed for legibility.
# Then, a simple pie chart is created to express the ratio between NUC sequences and MIT sequences (after being counted with SQL)
# Ideally, the ratio will always favor the NUC group as larger because it is over-represented in the original set
# For consistency, this chart is made using the ggplot2 library, even though it would technically be faster to make a simple pie chart
df2_nuc_samples_clean <- subset(df2_nuc_samples, select = -c(MCG.x, GVH.x, ALM.x, MIT.x, ERL.x, POX.x, VAC.x, NUC.x, Class.x, MCG.y, 
                                                             GVH.y, ALM.y, MIT.y, ERL.y, POX.y, VAC.y, NUC.y, Class.y))
count_nuc <- sqldf("select count(Class) from 'df2_nuc_samples_clean' where Class = 'NUC'")
count_mit <- sqldf("select count(Class) from 'df2_nuc_samples_clean' where Class = 'MIT'")
df2_pie <- data.frame(Class = c("NUC", "MIT"), Presence = c(as.integer(count_nuc), as.integer(count_mit)))
# The given algorithm for converting a bar graph to a pie chart was borrowed from the R Graph Gallery
# Link to instructions: https://r-graph-gallery.com/piechart-ggplot2.html
# For the plot generation, refer to the block at the bottom of the script

# 8
# we finally examine the df3 subset for membranous sequences
# Remove all the cytoskeletal sequences, then make a graph to represent the three remaining classes (sub-classes of cellular membrane)
df3_mem <- subset(df3, Class != "CYT")
me1_avg <- trunc(as.numeric(sqldf("select avg(MIT) from 'df3_mem' where Class = 'ME1'"))*10^2)/10^2
me2_avg <- trunc(as.numeric(sqldf("select avg(MIT) from 'df3_mem' where Class = 'ME2'"))*10^2)/10^2
me3_avg <- trunc(as.numeric(sqldf("select avg(MIT) from 'df3_mem' where Class = 'ME3'"))*10^2)/10^2
df3_mem_avgs <- data.frame(Class = c("ME1", "ME2", "ME3"), Average = c(me1_avg, me2_avg, me3_avg))
# For the plot generation, refer to the block at the bottom of the script

# 9
# use the df3_mem subset to caculate a sum of POX values
# Use a conditional code block to determine whether this value is useful for comparison or not
# If not, simply remove that attribute from the data set
pox_sum <- sum(df3_mem$POX)
if (pox_sum >= 5) {
  cyt_sum <- as.numeric(sqldf("select sum(POX) from 'df3' where Class = 'CYT'"))
  if (pox_sum > cyt_sum) {
    print("Cytoskeletal POX < Membrane POX")
    print(" ")
  } else {
    print("Cytoskeletal POX > Membrane POX")
    print(" ")
  }
} else {
  df3_mem_clean <<- subset(df3_mem, select = -c(POX))
}

# 10
# a subset of df3 with only cytoskeletal data will be created and then have its attributes changed
# The POX, ERL, and Class name (redundant) will be removed, and then that resultant set will be scrubbed for specific MIT/NUC values
# After the cleaning is done, the result set will be exported using my export function (included with my submission) at the bottom of the file
# If the export function doesn't work, please check the file path at the top of this source file, or un-comment the line at the bottom of this block
df3_cyt <- subset(df3, Class == "CYT")
df3_cyt_clean <- subset(df3_cyt, select = -c(POX, ERL, Class))
df3_cyt_clean <- subset(df3_cyt_clean, MIT <= 0.2 & NUC >= 0.22)
# write.table(df3_cyt_clean, file = "Cytoskeletal_Emergency.csv", append = FALSE, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

# The following block is a user interface, intended to be used to print out the tables in blocks according to their numbering
# Follow the instructions in an interactive environment/terminal to see the answer to each question
# In the event that you don't want this, simply copy-paste from the list of conditionals which answer you want to see into a console
userInput = " "
userFilename = " "
cat(paste("Please select the question you'd like to see:", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", sep = "\n"))
userInput <- readLines("stdin", 1)
userInput <- as.integer(userInput)
if (userInput == 1) {
  df1_clean
} else if (userInput == 2) {
  df2
  df3
} else if (userInput == 3) {
  df2_mit_dist
} else if (userInput == 4) {
  df2_mit_dist_final
} else if (userInput == 5) {
  bins_table
} else if (userInput == 6) {
  df2_nuc_samples
} else if (userInput == 7) {
  # The below is the pie chart for script 7
  ggplot(df2_pie, aes(x = "", y = Presence, fill = Class)) + geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) +
    ggtitle("Ratio of NUC Sequences to MIT Sequences among 3 Random Samples") + xlab(" ") + ylab(" ")
} else if (userInput == 8) {
  # The below is the bar graph for script 8
  ggplot(df3_mem_avgs, aes(x = "", y = Average, fill = Class)) + geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    ggtitle("Average MIT Distribution Score for ME1, ME2, and ME3 Membrane Sequence Samples") + xlab("Membrane Category") + 
    ylab("Average MIT Ratio") + scale_fill_hue(c = 50)
} else if (userInput == 9) {
  df3_mem_clean
} else if (userInput == 10) {
  cat(paste("Enter the name for your export file:", " ", sep = "\n"))
  userFilename <- readLines("stdin", 1)
  # The below is the export for script 10. This uses my export function included in UtilityFunctions.R
  # When used in a terminal, this will export the CSV for script 10 as a filename of the user's choice
  export(df3_cyt_clean, userFilename)
} else {
  paste("Sorry! Select from the assignment only.", " ", sep = "\n")
}
