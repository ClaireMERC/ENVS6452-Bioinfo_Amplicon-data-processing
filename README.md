# ENVS6452-Bioinfo_Amplicon-data-processing
Graph the relative abundance of bacterial species in samples from sequenced data sets of Illumina-sequenced paired-end files.

# Outline of work flow
From sequence data sets of Illumina-sequenced paired-end fastq files, this script makes an amplicon sequence variant (ASV) table with abundance, assign taxonomy to the output sequences, and produces different types of graph to represent the relative abundance. 

# Getting ready
Download data set and unzip
If not down yet, download dada2 (https://benjjneb.github.io/dada2/dada-installation.html) 

```{r}
library(dada2); packageVersion("dada2") #Invoke libraries

setwd("~/workingR/Assignment3_ASVclass") #set up your working directory
path <- setwd("~/workingR/Assignment3_ASVclass") #And make an object "path" to make it easier to point at your working directory when needed

list.files(path) #Print the files in your working directory to see if the data are successfully loaded and how they look like
```

Read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files (those manipulations MAY HAVE TO BE MODIFIED if your filename format is different).
```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE)) #Make an object with Forwards amplicons
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE)) #Make an object with Reverser amplicons
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #Extract sample names
```

# 1. Quality control and prep
Inspect reads quality profiles to determine truncation lenght.
Your reads must still overlap after truncation in order to merge them later.
Here we are working with 2x300 V4 sequence data and  the V4 region is around 390bp, so the forward and reverse reads should almost completely overlap and our trimming can be completely guided by the quality scores.

```{r}
plotQualityProfile(fnFs[1:2]) # visualizing the quality profiles of the forward reads
```
Based on the quality profile of your pumice sample only (not the negative control), choose trimming position where the mean quality score (green line) drops below 30 (here at 250)

```{r}
plotQualityProfile(fnRs[1:2]) #visualize the quality profile of the reverse reads
```
The reverse reads are usually of worse quality. DADA2 incorporates quality information into its error model which makes the algorithm robust to lower quality sequence, but trimming as the average qualities crash will improve the algorithm’s sensitivity to rare sequence variants.
As for forward read, choose trimming position where the mean quality score of your sample drops below 30 (here at 210)
```

Assign the filenames for the filtered fastq.gz files, place filtered files in filtered/ subdirectory, and make an objects with filtered Forward and Reverse amplicons.
```{r}
#Forward
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
#Reverse
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))  
names(filtRs) <- sample.names
```

Filter and trim:
Here we use standard filtering parameters: 
maxN=0 meaning DADA2 requires no Ns, which is a very strict parameter
truncQ=2, rm.phix=TRUE 
maxEE=2 sets the maximum number of “expected errors” allowed in a read (here, error allowed = 2 for F, 2 for R), which is a better filter than simply averaging quality scores.
```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,210),  #truncLen=c(260,210) numbers decided from quality trim=c(F,R) 
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, 
              multithread=FALSE) # set multithead as TRUE for Mac, FALSE on Windows

head(out) #look at the result
```

# 2. Evaluate errors rates (and trim again)
will take several minutes for the computer to process
```{r}
errF <- learnErrors(filtFs, multithread=FALSE) #Learn the error rates for the Forward reads:
errR <- learnErrors(filtRs, multithread=FALSE) #Learn the error rates for the Reverse reads:

plotErrors(errF, nominalQ=TRUE) #visualize the estimated error rates
```
The red line shows the error rates expected under the nominal definition of the Q-score. Here the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected. Everything looks reasonable and we proceed with confidence.

Sample Inference: We now apply the core sample inference algorithm to the filtered and trimmed sequence data
```{r}
dadaFs <- dada(filtFs, err=errF, multithread=FALSE) #for the Forward reads
dadaRs <- dada(filtRs, err=errR, multithread=FALSE) #for the Reverse reads

dadaFs[[1]] #Inspect the returned dada-class object
```

# 3. Merge paired end reads
forming contigs (full denoised sequence) by joining the forward and reverse reads together.
By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments).
```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

head(mergers) # Inspect the merger data frames
```
The mergers object is a list of data frames from each sample. Each data frame contains the merged sequence, its abundance, and the indices of the forward and reverse sequence variants that were merged. Paired reads that did not exactly overlap were removed by mergePairs, further reducing spurious output.
Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?


# 4. Evaluate unique sequences 
(eg: ASVs)

Construct sequence table: 
```{r}
seqtab <- makeSequenceTable(mergers) #Construct an amplicon sequence variant table (ASV) table
dim(seqtab)

table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths
```
The sequence table is a matrix with rows corresponding to the samples, and columns corresponding to the sequence variants.
The lengths of the merged sequences should all fall within the expected range for your amplicon. Sequences that are much longer or shorter than expected may be the result of non-specific priming. You can remove non-target-length sequences from your sequence table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256])

Remove chimeras:
Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.
The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on factors including experimental procedures and sample complexity.
```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) #multithread=TRUE on Mac and Windows
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
```
Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though). If most of your reads were removed as chimeric, upstream processing may need to be revisited. In almost all cases this is caused by primer sequences with ambiguous nucleotides that were not removed prior to beginning the DADA2 pipeline.

Track reads through the pipeline:
As a final check of the process, look at the number of reads that made it through each step in the pipeline
```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

head(track) #Look at the output
```
Good if the majority of the raw reads are kept. Outside of filtering, there should be no step in which a majority of reads are lost. If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon. If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification.


# 5. Assign taxonomy
The DADA2 package provides a native implementation of the naive Bayesian classifier method for this purpose. The assignTaxonomy function takes as input a set of sequences to be classified and a training set of reference sequences with known taxonomy, and outputs taxonomic assignments with at least minBoot bootstrap confidence.
Formatted training fastas and the Silva reference database can be found on this website: https://melindahiggins2000.github.io/N741bigdata/lesson78.html 

Download the silva_nr_v132_train_set.fa.gz and place it in your working directory
Download the species-assignment training fastas file silva_species_assignment_v132.fa.gz and place it in your working directory

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "~/workingR/miseqsopdata/MiSeq_SOP/silva_nr_v132_train_set.fa.gz", multithread=FALSE)

#Make species level assignments (based on exact matching between ASVs and sequenced reference strains) 
taxa <- addSpecies(taxa, "~/workingR/miseqsopdata/MiSeq_SOP/silva_species_assignment_v132.fa.gz")

#Remove sequence rownames for display only
taxa.print <- taxa
rownames(taxa.print) <- NULL

head(taxa.print) #inspect the taxonomic assignments
```
If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) and see if this fixes the assignments. If using DECIPHER for taxonomy, try IdTaxa (..., strand="both").

Save your objects as csv file (to save them permanently on your computer)
```{r}
write.csv(taxa.print, "~/workingR/Assignment3_ASVclass/taxa_print.csv")
write.csv(seqtab.nochim, file = "~/workingR/Assignment3_ASVclass/seqtab.csv")
write.csv(track, file= "~/workingR/Assignment3_ASVclass/track.csv")
```
It is now safe to quit R or clean your R environment as everthing you need for the following steps had been saved as csv files. 


# 6. Formate your output data

Starting from saved csv files, import them as objects:
```{r}
taxa <- read.csv(file = "~/workingR/Assignment3_ASVclass/taxa_print.csv")
seqtab.nochim <- read.csv(file = "~/workingR/Assignment3_ASVclass/seqtab.csv", header = FALSE)
track <- read.csv(file = "~/workingR/Assignment3_ASVclass/track.csv")
```

Transpose the seqtab.nochim data to look at it as a column (flipping rows into coloumns)
```{r}
flipped_seqtab.nochim<- as.data.frame(t(seqtab.nochim))
View(flipped_seqtab.nochim)
```
In the flipped csv, the header becomes random letters and moves the sample names to a new row (first row) this is something you need to edit to make it easier and more organized to use with downstream analysis tools like phyloseq
Copy the first row into the header and then delete the first row (in two steps to double check either version)
```{r}
colnames(flipped_seqtab.nochim) <- flipped_seqtab.nochim[1,] #copy the first row 
View(flipped_seqtab.nochim)
flipped_seqtab.nochim <- flipped_seqtab.nochim[-1,] #then delete the first row
View(flipped_seqtab.nochim)
```
Next change the names of the sequences to "ASVs" so it is more digestable than the nucleotide sequence itself. In flipped_seqtab.nochim.csv files -> take the first column, first row name and paste csv for all first column rows >> changes all the sequences in the first column to ASV1,2,3,4.. etc instead, so it is easier to use later on with data 
```{r}
rownames(flipped_seqtab.nochim) <- paste0("ASV", 1:nrow(flipped_seqtab.nochim))
 
flipped_seqtab.nochim_forself <- flipped_seqtab.nochim[,-1] #remove the sequences column: save as flipped_seqtab.nochim_forself
```

Save the files in case it is useful later
```{r}
#saves your flipped_seqtab.nochim file with your taxa data as one data sheet
write.csv(flipped_seqtab.nochim, file = "~/workingR/Assignment3_ASVclass/flipped_seqtab.nochim.csv")
write.csv(flipped_seqtab.nochim_forself, file ="~/workingR/Assignment3_ASVclass/flipped_seqtab.nochim_forself.csv")

#then save your flipped seqtab.nochim file with your taxa data as one data sheet 
OTUabund<-cbind(flipped_seqtab.nochim,taxa)
write.csv(OTUabund,file="C:/Users/Claire/Desktop/workingR/Assignment3_ASVclass/OTUabund.csv")

#Look at your data to make sure everything worked
taxa<-taxa[-1]
View(taxa)
```

# 7. Make the graph(s)

If not done yet, download the packages needed to evaluate data and make graph
To install phyloseq, you might have to look on bioconductor.org > packages
```{r}
#Invoke those libraries
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(phyloseq)
library(Biostrings)
```

Formate your data
```{r}
taxmat <- as.matrix(taxa) #Call your phyloseq formatted object "taxmat"

#Make the OTU table in a way the format that phyloseq needs. Call this object "otumat".
otumat <-flipped_seqtab.nochim[,-1]
view(otumat)

#Convert all the files to matrix format so they can be used in phyloseq
class(otumat)
class(taxmat)
otumat <- as.matrix(otumat)
taxmat <-as.matrix(taxmat)

#inspect and confirm that they are in matrix format and ready to be used in phyloseq
rownames(otumat) <- paste0("ASV", 1:nrow(otumat))
rownames(taxmat) <- paste0("ASV", 1:nrow(otumat))
 
#make sure that R recognizes that the OTU data is numeric, not character data
class(otumat)<-"numeric"
```

Now that you have matrices, use phyloseq to continue analysis
```{r}
OTU = otu_table(otumat, taxa_are_rows = TRUE) #tell phyloseq where your ASVs files are
TAX = tax_table(taxmat) #tell phyloseq where your taxa files are

#tell phyloseq to put it all together ( sample names, ASVs and taxa) 
physeq = phyloseq(OTU, TAX)
physeq
sample_names(physeq)
samplenames<-sample_names(physeq)
```

Make the plot(s) you'd like :

# option 1 : bar graph of relative abundance of PHYLA
```{r}
#Plotting bar using the plot_bar function of phyloseq
p<-plot_bar(physeq, fill = "Phylum")

#Remove the lines between each ASV using the stacking of ggplot's geom_bar function
pstacked<- p + geom_bar(aes(fill=Phylum), stat="identity", position="stack")

#Merging the ASVs of each phyla together to make it easier to interpret the plot:
#First, use the "tax_glom" of phyloseq to glom together taxa based on the column of your choosing. Here we are glomming together the taxa in the phylum column to see all of the phylum abundances.
ps_phylum <- tax_glom(physeq, "Phylum")
 
#Plot the ps_phylum graph using plot-bar function to look at the differences between the separate taxa with each ASV being represented and glomming the ASVs together
plot_bar(ps_phylum, fill = "Phylum")

#Looking at relative abundance is more helpful and insightful than absolute abundance. Make a transformation of sample counts and our glommed taxa from previous step and put in calculation of ASV over total ASVs. Then we psmelt it so it is easier the plot and factor it into a table of the relative abundance of the phyla:
 
#after you glom together your taxa by phylum, make a table of relative abundance by talling up each taxa, and dividing by the total taxa (eg: what percentage of the total is each phylum in each sample)
ps_phylum_relabun <- transform_sample_counts(ps_phylum, function(ASV) ASV/sum(ASV))

#Then use psmelt to melt away the phyloseq formatting and make it easier for plotting
taxa_abundance_table_phylum <- psmelt(ps_phylum_relabun)

#factor the values of Phylum
taxa_abundance_table_phylum$Phylum<-factor(taxa_abundance_table_phylum$Phylum)

#save the abundance table for easy access to look at again later
write.csv(taxa_abundance_table_phylum, file ="~/workingR/Assignment3_ASVclass/RelativeAbundance_Phylum.csv")

#label with graph and axis titles
title = "Relative Abundance of Phyla in Pumice Rock in the South Pacific "
p_realabun<-plot_bar(ps_phylum_relabun, fill = "Phylum", title=title) + ylab("Relative Abundance (%)") +

#Place the sample name parallel to the x axis  
theme(axis.text.x = element_text(angle = 0))

#take the previous relative abundance graph and stack the phyla to remove the lines in between each phyla to make more seamless 
p_abun_stacked<- p_realabun + geom_bar(aes(fill=Phylum), stat="identity", position="stack")
p_abun_stacked
```
#export these graph (plot window > export...)


# Option 2: bar graph of relative abundance of ORDERS
```{r}
#glom together taxa in the order column to graph them together
ps_order <- tax_glom(physeq, "Order")

#plot this to see what you get after glomming the taxa of the orders
plot_bar(ps_order, fill ="Order")

#make a transformation of sample counts and our glommed taxa from previous step and put in calculation of ASV over total ASVs
ps_order_relabun <- transform_sample_counts(ps_order, function(ASV) ASV/sum(ASV))

#psmelt it so it is easier to plot and factor it into a table of the relative abundance of the orders
taxa_abundance_table_order <- psmelt(ps_order_relabun)

#factor the values of Orders
taxa_abundance_table_order$Order<-factor(taxa_abundance_table_order$Order)

#save the abundance table for easy access to look at again later
write.csv(taxa_abundance_table_order, file = "~/workingR/Assignment3_ASVclass/RelativeAbundance_Order.csv")
 
#Using the abundance table made using the phyloseq analysis, plot this relative abundance of orders table as a bar graph and add a title and axis labels 
title = "Relative Abundance of Orders in Pumice Rock Samples in South Pacific"
o_realabun <- plot_bar(ps_order_relabun, fill = "Order", title=title) + ylab("Relative Abundance (%)") +

#Place the sample name parallel to the x axis  
theme(axis.text.x = element_text(angle = 0)) +

#Change the size of the graph
theme(aspect.ratio = 4/1)

#Use the stacked ability in geom_bar function to remove the lines between each order to make the plot more seamless
o_abun_stacked<- o_realabun + geom_bar(aes(fill=Order), stat="identity", position="stack")
o_abun_stacked
```

# Option 3: POINT GRAPH of relative abundance of orders
```{r}
library(scales) #invoke one more library

#glom together taxa in the order column to graph them together
ps_order <- tax_glom(physeq, "Order")

#make a transformation of sample counts and our glommed taxa from previous step and put in calculation of ASV over total ASVs
ps_order_relabun <- transform_sample_counts(ps_order, function(ASV) ASV/sum(ASV))

#psmelt it so it is easier to plot and factor it into a table of the relative abundance of the orders
taxa_abundance_table_order <- psmelt(ps_order_relabun)
taxa_abundance_table_order$Order<-factor(taxa_abundance_table_order$Order)

#save the abundance table for easy access to look at again later
write.csv(taxa_abundance_table_order, file ="~/workingR/Assignment3_ASVclass/RelativeAbundanceOrder.csv")
 
#Using the abundance table made using the phyloseq analysis, plot this relative abundance of orders table as a point graph
taxa_abundance_table_order$Abundance[taxa_abundance_table_order$Abundance == 0] <- NA #make null value = no point
point_plot <- ggplot(taxa_abundance_table_order,aes(Sample,Order)) + #make graph
geom_point(aes(size=Abundance, fill=Order),shape=21,color="black") + #shape point
theme(panel.grid.major=element_line(linetype=1,color="grey"), panel.background = element_blank()) + #format background
labs(title="Relative Abundance of Orders in Pumice Rock Samples in South Pacific") + ylab("Order") + xlab("Sample") + #add titles
scale_size(name = "Relative abundance") + #legend
theme(aspect.ratio = 6/1) #adapt graph dimensions

point_plot #show graph
```

# Option 4: PIE CHART of relative abundance of order
```{r}
#glom together taxa in the order column to graph them together
ps_order <- tax_glom(physeq, "Order")

#make a transformation of sample counts and our glommed taxa from previous step and put in calculation of ASV over total ASVs
ps_order_relabun <- transform_sample_counts(ps_order, function(ASV) ASV/sum(ASV))

#psmelt it so it is easier to plot and factor it into a table of the relative abundance of the orders
taxa_abundance_table_order <- psmelt(ps_order_relabun)
taxa_abundance_table_order$Order<-factor(taxa_abundance_table_order$Order)

#save the abundance table for easy access to look at again later
write.csv(taxa_abundance_table_order, file ="~/workingR/Assignment3_ASVclass/RelativeAbundanceOrder.csv")
 
#Using the abundance table made using the phyloseq analysis, plot this relative abundance of orders table as a pie chart
pie_chart <- ggplot(taxa_abundance_table_order, aes(x="", y=Abundance, fill=Order)) +
  geom_bar(stat="identity", width=1) + #make your graph
  coord_polar("y", start=0) + #set y values positions
facet_wrap(facets=vars(Sample)) + #facet to make on pie chart per sample
theme_void() + #erase background
labs(title="Relative Abundance of Orders in Pumice Rock Samples in South Pacific") +
theme(plot.title = element_text(hjust = 0.5))+ #titles 
theme(legend.position = "bottom",legend.title.position="top") #format legend

pie_chart #show graph
```
