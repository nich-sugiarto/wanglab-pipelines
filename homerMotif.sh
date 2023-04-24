#!/bin/bash

# USAGE: 
# This pipeline takes in one positional argument:
# 	$1 - target folder

# This script looks at all the bed files located in the user's desired folder,
# and runs HOMER's findMotifs program offshore on the discovery HPC. 

# Requires that HOMER is available in you PATH variable

mkdir -p PBS
mkdir -p log

folder=$(pwd)

# If a name is not provided
if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target folder
  exit 1
fi

cd $1

for file in *.bed; do
	base=$(basename "$file" ".bed")
	mkdir -p ${folder}/homer_sizegiven/${base}

  cat >${folder}/PBS/${base}_motif_sizegiven'.R' <<EOF
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
cleanNames <- str_split(Table\$motif_name, pattern = \\\\(")  

# Found this one online. It takes the first sub element from the list of lists
cleanNames <- sapply(cleanNames, "[[", 1)

# Overrides the old names with the new names
Table\$motif_name <- cleanNames

# Version you can paste right into the terminal (minus the comments out ofc)
# cleanNames <- str_split(Table\$motif_name, pattern = "\\(")
# cleanNames <- sapply(cleanNames, "[[", 1)
# Table\$motif_name <- cleanNames

# Version where you can have it all on one line:
# Table\$motif_name <- sapply(str_split(Table\$motif_name, pattern = "\\\\("), "[[", 1)   

Table_100 <- head(Table,100)

motifLvls <- as.list(Table_100\$motif_name)

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
	label = motif_name), min.segment.length = 0.0000001, size = 2) +
	geom_point(data=filtered_table,aes(x=motif_name, y=log_p_value, color="red"))

ggsave("labeledTestPlot.pdf")
EOF

	cat >${folder}/PBS/${base}_motif_sizegiven'.pbs' <<EOF
#!/bin/bash -l

#SBATCH --job-name=${base}_motifs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=16:00:00
#SBATCH -o ${folder}/log/${base}_motifs_%j.txt -e ${folder}/log/${base}_motifs_%j.err.txt
# Change to job working directory
cd ${folder}
#---------------------End of header---------------------#

findMotifsGenome.pl $1/${base}.bed hg38 homer_sizegiven/${base}/ -size given

Rscript ${folder}/PBS/${base}_motif_sizegiven'.R'

EOF

	sbatch ${folder}/PBS/${base}_motif_sizegiven.pbs
done
