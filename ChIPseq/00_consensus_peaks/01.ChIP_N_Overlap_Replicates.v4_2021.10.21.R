# script takes a list of BED files (replicates), 
#    and filters their elements to remove those not overlapping with at least N other replicates.
#
# then, the consensus peak set is constructed, using reduce() command in GenomicRanges.
# 
# 

# 
# initial BED files and output BED files are in 0-based coordinates
# internal calculations with GRanges are in 1-based coordinates


options(scipen=999) # prevent using scientific notation when writing the files

# helper function to read bed files into R and remove extra columns
read_bed <- function(filename, header =F, sep ="\t") {
  df <- read.table(filename, sep="\t", header = F, quote = "", stringsAsFactors = F )
  df.short <- df[,1:3] # only "chr", "start" and "end" columns
  names(df.short) <- c("chr","start","end")
  return (df.short)
}

# function: from a GRanges object, create a BED-formatted data.frame
df_from_gr <- function(gr ) {
  df <- data.frame(chr = seqnames(gr), start = start(gr) -1, end = end(gr))
  # -1 :  1-based to 0-based coord conversion 
  return (df)
}

# function: 
# from list of GRanges, constuct a list of dataframes with extra column = N replicates overlapping w peak)
n_overl_replicates <- function(abc.gr) {
  abc.new = vector("list", length(abc.gr))
  
  for (n in (1:length(abc.gr)))  {
    #a.tmp = data.frame(abc.gr[n])
    a.tmp = df_from_gr(abc.gr[[n]]) 
    a.tmp$overl = 0
    a.tmp
    others.gr = abc.gr[-n]
    for (i.gr in others.gr) { 
      
      count.tmp = countOverlaps(abc.gr[[n]] , i.gr)
      #print(count.tmp)
      count.tmp = ifelse(count.tmp >= 1, 1, 0)
      # EXPLANATION:
      # x = c(2,3,1,0,1,3)
      # y = ifelse(x >= 1, 1, 0)
      # y
      # [1] 1 1 1 0 1 1
      
      #print(count.tmp)
      a.tmp$overl = a.tmp$overl + as.vector(count.tmp)
      #print ("----")
    }
    #print(a.tmp)
    abc.new[[n]] <- a.tmp
    
  }
  return(abc.new) # list of data.frames for each replicate, with last column = N overlapping replicates
}

# function: write BED file (data.frame) as a file to disk
write_bed <- function(bed, filename, col.names =F, sep ="\t") {
  bed.short <- bed[,1:3]  
  write.table(bed.short, filename, sep="\t", col.names = F, row.names = F, quote = F )
}

# supplemental function: copy to clipoard, to paste into Excel etc
# use: wc() , or wc(command or object)
wc <- function(x = .Last.value) {
  clipr::write_clip(x)
}


###### START 

# choose the necessary group of datasets (cell type)
# CHOOSE ONE:
listFiles = list.files(pattern = "(.*)GLU.bed$")
listFiles = list.files(pattern = "(.*)SOX.bed$")
listFiles = list.files(pattern = "(.*)OLIG.bed$")

listFiles = list.files(pattern = "(.*)NeuN.bed$")
listFiles = list.files(pattern = "(.*)Olig.bed$")
listFiles = list.files(pattern = "(.*)Astro.bed$")
listFiles = list.files(pattern = "(.*)MG.bed$")

# check files
listFiles

# create list of dataframes
list_beds <- lapply(listFiles, read_bed)
# create list of GRange objects
list_gr <- lapply(list_beds, function(x) with(x, GRanges(chr, IRanges(start+1, end), strand = "*")) ) # +1 : 0- to 1-based coords

# OPTIONAL:
# print number of peaks in replicates
lapply(list_beds, nrow)

# for each peak in each replicate dataset, calculate N replicates with which it overlaps.
out = n_overl_replicates(list_gr)

# OPTIONAL
# depending on cell type, save the list with a new name
out.glu <- out
out.gaba <- out
out.olig <- out
out.mg = out
out.neur = out
out.astro = out
out.olig = out

#########################
# Calculate N of peaks in __reduced__ datasets (one GRange merged from all replicates) for each threshold values
# for a selected threshold Nthresh (minimum number of replicates a peak should overlap with), write down the consensus table of peaks

# CHANGE AS NEEDED
Nthresh = 4
#Nthresh = 2

out.full <- out

maxNoverl = length(out.full) - 1
print(maxNoverl)

filt <- vector("list", length(list_gr))

for (thresh in (0:maxNoverl ) ) {
  for (i in (1:length(out.full))) {
    tmp = out.full[[i]]
    tmp = tmp[tmp$overl >= thresh , ]
    filt[[i]] <- tmp
  }
  #print(filt[[i]])
  list_df_from_filt <- lapply(filt, function(x) with(x, GRanges(chr, IRanges(start+1, end), strand = "*" , overlap = overl)) ) 
  
  list_gr_from_filt = GRangesList(list_df_from_filt)
  unlist_filt = unlist(list_gr_from_filt)
  sort_filt = sort(unlist(list_gr_from_filt))  
  reduced_filt = reduce(sort_filt)    

  dim_filt = length(reduced_filt)
  print(thresh +1) # print the N of overlapping replicates
  print (dim_filt) # print the number of peaks in the consensus peak set for this N
  
  # for the selected Nthresh of overlaps, write down the GRanges object for consensus peaks 
  if (thresh == Nthresh) {  
    gr_consens = reduced_filt
  }
}

# convert the consensus peak set from GRanges to data.frame
gr_consens.df = df_from_gr(gr_consens)
dim(gr_consens.df)

# save file to disk
# (subsitute the name as needed)
write_bed(gr_consens.df, "OLIG.overl5.bed")

# Calculate the lengths of peaks
gr_consens.df$lengs = gr_consens.df$end - gr_consens.df$start +1
# copy data.frame to clipboard
wc(gr_consens.df)

##################

#########################
# OPTIONAL:
# create list of _filtered_ dataframes, applying a threshold  
filt <- vector("list", length(list_gr))
# substitute as needed:
thresh <- 4 # overlapping with at least 4 other replicates


for (i in (1:length(out))) {
  tmp = out[[i]]
  tmp = tmp[tmp$overl >= thresh , ]
  filt[[i]] <- tmp
}
# calculate the number of peaks in filtered datasets
lapply(filt, nrow)

# write the bed files to disk.
# each file represents the peaks from a particular replicate which have >=N overlaps with other replicates.
celltype = "ASTR"
celltype = "MG"
celltype = "MG"
celltype = "OLIG"
celltype = "NEUR"

# list of file names:
list_output_files <- paste0(celltype, 1:length(filt), ".bed")
list_output_files
# write files to disk
for (i in 1:length(filt)) {
  write_bed(filt[[i]], list_output_files[i])
}
#########################