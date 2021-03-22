##### RUNS IN SILICO DIGESTIONS OF MUS MUSCULUS REFERENCE GENOME AND CREATES PLOTS OF FRAGMENT DISTRIBUTION ALONG EACH CHROMOSOME. ##### 
#clear workspace
rm(list = ls())

## load the required libraries ##
library("SimRAD")
library("karyoploteR")
library("parallel")
#working directory
#setwd("~/Documents/RADseq_simulations/M.musculus_SimRAD/data/")

## Specify restriction enzyme recognition sites
#Restriction enzyme 1 - sbfI 
cs_5p1 <- "CCTGCA"
cs_3p1 <- "GG"
#Restriction enzyme 2 - mseI 
cs_5p2 <- "T"
cs_3p2 <- "TAA"

#create a list of input reference sequences for ALL M.musculus chromosomes
ReferenceFilesList <- list.files(pattern = "*.fa")

### LOOP THROUGH EACH CHR FASTA FILE FOR IN SILICO DIGESTION ###

#function for in silico digestion
in_silico_digestion <-  function(x) {
  for (i in 1:length(x)) {
    
    # read in the reference sequence
    assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".ref"), 
           ref.DNAseq(FASTA.file = x[i], subselect.contigs = FALSE))
    #Mm_chr <- ref.DNAseq(FASTA.file = x[i], subselect.contigs = FALSE)
    
    #in silico digestion
    assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".dig"), 
           insilico.digest(get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".ref")), 
                           cs_5p1, cs_3p1, cs_5p2, cs_3p2)
           )
    #Mm_chr.dig <- insilico.digest(Mm_chr, cs_5p1, cs_3p1, cs_5p2, cs_3p2)
    
    # simulate library construction (assigning a name based on the chromosome)
    assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".lib"), 
           adapt.select(get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".dig")), 
                        type = "AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2))#, envir = .GlobalEnv)
    
    #create a list of all the objects from the digest
    assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i]),
           list("reference" = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".ref")),
                #"digest" = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".dig")),
                "library" = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".lib"))
                )
           )
  }
  #return the list of all objects from the digestion
  return(get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i])))
}

#run the digestion for each chromosome in parallel 
system.time({libraries <- mclapply(ReferenceFilesList, in_silico_digestion, mc.cores = detectCores())})

#extract and name all the libraries and references by chromosome
for (i in 1:length(libraries)){
  assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".ref"), libraries[[i]]$reference)
  assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".lib"), libraries[[i]]$library)
}
rm(libraries)

### RETRIEVE THE FRAGMENT POSITIONS WITHIN EACH REFERENCE CHROMOSOME IN PARALLEL
#create a list of library names

#create a function to retrieve fragment positions
fragment_positions <- function(x) {
  for (i in 1:length(x)) {
    
    #create temporary empty vectors to write to
    start.tmp <- c()
    end.tmp <- c()
    
    #match the fragments to the reference chromosome
    for (j in 1:length(as.character(get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".lib"))))){
      if (width(get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".lib"))[j]) <= 20000){
        index <- matchPattern(as.character(get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".lib")))[j],
                              get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i], ".ref")))
      
        start.tmp[j] <- start(index)
        end.tmp[j] <- end(index)  
      }
      
    }
    
    #assign the start and end positions to the designated vector
    assign(paste0("start_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i]), start.tmp)
    assign(paste0("end_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i]), end.tmp)  
  }
  
  #create a list of start and end positions from each chromosome
  assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i]),
         list("start" = get(paste0("start_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i])),
              "end" = get(paste0("end_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i]))
         )
  )
  return(get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", x))[i])))
  
}
#retrieve fragment positions for each chromosome in parallel
system.time({positions <- mclapply(ReferenceFilesList, fragment_positions, mc.cores = detectCores())})

### PLOTTING ###
#create IRanges objects of start and end positions, and width of fragments for each chromosome
for (i in 1:length(positions)){
  positions[[i]]$start <- na.omit(positions[[i]]$start)
  positions[[i]]$end <- na.omit(positions[[i]]$end)
  assign(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"), 
         IRanges(start = positions[[i]]$start, end = positions[[i]]$end))
}
rm(positions)

#plot an ideogram
kp <- plotKaryotype(genome = "mm10")

#plot the fragments on the ideogram to see their distribution along each chromosome
for (i in 1:length(ReferenceFilesList)){
  kpPlotRegions(kp, data = GRanges(seqnames = Rle(paste0(sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i])), 
                                   ranges = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))), col = "red")
}

#plot the density of fragment distribution along each chromosome
for (i in 1:length(ReferenceFilesList)){
  kpPlotDensity(kp, data = GRanges(seqnames = Rle(paste0(sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i])), 
                                   ranges = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))))
}

# #plot the coverage of fragment distribution along each chromosome
for (i in 1:length(ReferenceFilesList)){
  kpPlotDensity(kp, data = GRanges(seqnames = Rle(paste0(sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i])),
                                   ranges = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))), window.size = 1)
}

#create a dataframe of fragment length distributions
length_distribution <- data.frame("Chromosome" = c(rep("chr1", length(Mm_chr1.lib)), rep("chr2", length(Mm_chr2.lib)), rep("chr3", length(Mm_chr3.lib)), 
                                                   rep("chr4", length(Mm_chr4.lib)), rep("chr5", length(Mm_chr5.lib)), rep("chr6", length(Mm_chr6.lib)), 
                                                   rep("chr7", length(Mm_chr7.lib)), rep("chr8", length(Mm_chr8.lib)), rep("chr9", length(Mm_chr9.lib)), 
                                                   rep("chr10", length(Mm_chr10.lib)), rep("chr11", length(Mm_chr11.lib)), rep("chr12", length(Mm_chr12.lib)), 
                                                   rep("chr13", length(Mm_chr13.lib)), rep("chr14", length(Mm_chr14.lib)), rep("chr15", length(Mm_chr15.lib)), 
                                                   rep("chr16", length(Mm_chr16.lib)), rep("chr17", length(Mm_chr17.lib)), rep("chr18", length(Mm_chr18.lib)), 
                                                   rep("chr19", length(Mm_chr19.lib)), rep("X", length(Mm_chrX.lib)), rep("Y", length(Mm_chrY.lib))),
                                  "fragment_length" = width(c(Mm_chr1.lib, Mm_chr2.lib, Mm_chr3.lib, Mm_chr4.lib, Mm_chr5.lib, Mm_chr6.lib, Mm_chr7.lib, 
                                                              Mm_chr8.lib, Mm_chr9.lib, Mm_chr10.lib, Mm_chr11.lib, Mm_chr12.lib, Mm_chr13.lib, Mm_chr14.lib, 
                                                              Mm_chr15.lib, Mm_chr16.lib, Mm_chr17.lib, Mm_chr18.lib, Mm_chr19.lib, Mm_chrX.lib, Mm_chrY.lib)))

#### EXPLORE THE DATA ####

#plot fragment size distributions
plot(length_distribution)#per chromosome
hist(length_distribution$fragment_length, breaks = length(length_distribution$fragment_length)/20)#for all sequences

#calculate mean and quantiles of fragment lengths for each chromosome (ALL FRAGMENTS)
tapply(length_distribution$fragment_length, length_distribution$Chromosome, mean)#mean
tapply(length_distribution$fragment_length, length_distribution$Chromosome, quantile)#quantiles

#calculate mean fragment lengths for <2000 bp... ONLY 743 FRAGMENTS GREATER THAN 2000bp = 0.03% of total fragments
length_distribution_short <- length_distribution[which(length_distribution$fragment_length < 2000), ]
plot(length_distribution_short)#plot distribution
tapply(length_distribution_short$fragment_length, length_distribution_short$Chromosome, mean)#mean
tapply(length_distribution$fragment_length, length_distribution$Chromosome, quantile)#quantiles

#### SIZE SELECTION ####
#function to give a histogram of fragment length distribution for the whole genome (not per chromosome)
sizeSelect <- function (fragmentLength, min.size, max.size, graph = TRUE, verbose = TRUE) 
{
  ssel <- fragmentLength[(fragmentLength) < max.size & (fragmentLength) > 
                      min.size]
  if (verbose == TRUE) {
    cat(length(ssel), " fragments between ", min.size, " and ", 
        max.size, " bp \n", sep = "")
  }
  if (graph == TRUE) {
    bk <- hist((fragmentLength), breaks = length(fragmentLength)/20, 
               plot = FALSE)$breaks
    hist((fragmentLength), border = "grey75", col = "grey75", 
         breaks = bk, main = "", xlab = "Locus size (bp)", 
         ylab = "Number of loci", xlim = c(0,2500))
    hist((ssel), border = "red", col = "red", add = TRUE, 
         breaks = bk)
    text(mean(c(min.size, max.size)), max(hist((ssel), breaks = bk, plot = FALSE)$counts), pos = 4, 
         labels = paste(length(ssel), " loci between ", min.size, " and ", max.size, " bp", sep = ""), 
         col = "red", cex = 1.0, font = 2)
  }
  return(ssel)
}

##plot desired fragment length disrtibution across whole genome
#select min and max fragment sizes
max_fragment_size <- 600
min_fragment_size <- 300

#how many fragments within a size range (in bp)?
nrow(length_distribution[length_distribution$fragment_length <max_fragment_size & length_distribution$fragment_length >min_fragment_size,])

ssel <- sizeSelect(length_distribution$fragment_length, min.size = min_fragment_size, max.size = max_fragment_size)

#calculate the percentage of the genome targeted based on the number of size selected loci
(sum(ssel)/sum(nchar(c(Mm_chr1.ref, Mm_chr2.ref, Mm_chr3.ref, Mm_chr4.ref, Mm_chr5.lib, Mm_chr6.ref, Mm_chr7.ref, Mm_chr8.ref, Mm_chr9.ref, Mm_chr10.ref, 
            Mm_chr11.ref, Mm_chr12.ref, Mm_chr13.ref, Mm_chr14.ref, Mm_chr15.ref, Mm_chr16.ref, Mm_chr17.ref, Mm_chr18.ref, Mm_chr19.ref, 
            Mm_chrX.ref, Mm_chrY.ref))))*100

#plot the distribution of the size selected loci 
kp <- plotKaryotype(genome = "mm10")

for (i in 1:length(ReferenceFilesList)){
  kpPlotRegions(kp, data = GRanges(seqnames = Rle(paste0(sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i])), 
                                   ranges = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))
                                   [get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))@width > min_fragment_size & #insert min size here
                                       get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))@width < max_fragment_size]), #insert max size here
                avoid.overlapping = F , col = "red", r1 = 0.1)
}

kp <- plotKaryotype(genome = "mm10")
for (i in 1:length(ReferenceFilesList)){
  kpPlotDensity(kp, data = GRanges(seqnames = Rle(paste0(sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i])), 
                                   ranges = get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))
                                   [get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))@width > min_fragment_size & #insert min size here
                                       get(paste0("Mm_", sub(".fa", "", sub("mm_ref_GRCm38.p4_", "", ReferenceFilesList))[i], ".pos"))@width < max_fragment_size]))#insert max size here
}
