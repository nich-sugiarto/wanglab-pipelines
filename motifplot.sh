if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  echo Where the meta file is tab-delimited, and formated as follows:
  echo locationOfMetaCSV	DESeqResultTable
  exit 1
fi

mkdir -p goPlots
mkdir -p PBS
mkdir -p log

folder=$(cd "$(dirname "$0")";pwd)  # Save current folder as a variable

libraryText=$1

while IFS=$'\t' read -r -a varArray; do
	metascape=${varArray[0]}  # Location for metascape table
	de=${varArray[1]}

	mName=$(basename "$metascape" ".csv")

	cat >${folder}/PBS/${mName}'_GO.R' <<EOF
library(ggplot2)
library(dplyr)
library(tidyverse)
library(janitor)
library(viridis)
library(stringr)
library(ggrepel)
library(RColorBrewer)

#2/16/23
#top 100 motifs
#label/annotate significant (below E-50)
#fix labels (substrings)
#Change color of annotated motif name dots (firebrick3 (scale_color_manual))

Table <- read.table("knownResults.txt", comment.char="", sep="\t", header = TRUE) %>%
	clean_names() %>%
	unique() %>%
	arrange(log_p_value)

# Split the motif_names by the parentheses. 
# The backslashes in front of the "(" are necessary so that it's read as plain text
# and not some special character
# The back slash in front of the backslash is so that when you cat this in as an R script,
# It reads in the backslash as it is and not its own special character
# Tbh, just ask me later
# Returns a list of lists
cleanNames <- str_split(Table\$motif_name, pattern = "\\\\(")  

# Found this one online. It takes the first sub element from the list of lists
cleanNames <- sapply(cleanNames, "[[", 1)

# Overrides the old names with the new names
Table\$motif_name <- cleanNames

# Version you can paste right into the terminal (minus the comments out ofc)
# cleanNames <- str_split(Table$motif_name, pattern = "\\(")
# cleanNames <- sapply(cleanNames, "[[", 1)
# Table$motif_name <- cleanNames

# Version where you can have it all on one line:
# Table\$motif_name <- sapply(str_split(Table\$motif_name, pattern = "\\\\("), "[[", 1)

Table_100 <- head(Table,100)

motifLvls <- as.list(Table_100$motif_name)

Plot <- ggplot(data=Table_100, aes(x=factor(motif_name, level = unique(motifLvls)), y=log_p_value)) + 
	geom_point() +
	scale_y_reverse() +
	theme_bw() +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank()) + 
	ylab("Log_P_Value") +
	ggtitle("H514 AH_K27A Ranked P Values")

ggsave("testPlot.pdf")

filtered_table <- Table_100 %>%
	filter(log_p_value<(-50)) %>%
	arrange(log_p_value)

labeled_plot <- Plot +
	geom_label_repel(data=filtered_table,
	mapping = aes(x=factor(motif_name, level = unique(motifLvls)), y=log_p_value, 
	label = motif_name), min.segment.length = 0.0000001, size = 2)
	
ggsave("labeledTestPlot.pdf")

filtered_table$motif_name <- list_motif_names


#SEPARATE (for metascape plot)
# motif names left axis, one for each dot (split string by /, make list, take first element)
## list = str_split(filtered_table$motif_name, "(")
## sapply(list,"[[,1)

# Log(pvalue) - color <-- pAdj? 
# Size of dot - % of sequences

