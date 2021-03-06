
# This script generates a PRELIMINARY training set using the best quality detected and deconvoluted peak groups as targets and producing their corresponding decoys. 
# For decoy generation, the peak groups are associated by pairs and used as targets. 
# For each pair of targets, A and B, a pair of decoys is generated by keeping the same precursor and its properties and swapping the m/z values of 60% of the fragments.
# A list with generated targets and decoys is saved in Skyline format.
#
# INPUT:
#   - PeakID...csv (aligned matrix with deconvoluted MS2s)
# OUTPUT:
#   - SkylineTargetList_preliminary-target-decoy_[sample]_[date].csv
#   - TargetDecoy-preliminary-figures_[sample]_[date].pdf

basePathData = "/Users/xxx/Documents" # <--- update this path!

setwd(file.path(basePathData,"data/Asper"))
alignmentTable = "xMSDIAL/PeakID_0_20221211731.txt"

#setwd(file.path(basePathData,"data/Pput"))
#alignmentTable = "xMSDIAL/PeakID_0_202111201516.txt"

#setwd(file.path(basePathData,"data/Rhodo"))
#alignmentTable = "xMSDIAL/PeakID_0_2022126180.txt"


sampleString = basename(getwd())
minSignalToNoise = 15 
minMzDistanceMSMS = 0.4
maxCountFragments = 16 
ExplicitRTwindow = 2 # For Skyline target list

pdfFileName = paste("TargetDecoy-preliminary-figures_", sampleString, "_", Sys.Date(), ".pdf", sep='')

dat = read.csv(file = alignmentTable, sep = '\t', stringsAsFactors = FALSE, skip = 3)
colnames(dat) = gsub(" |/|\\(|\\)|\\%", "", dat[1,])
dat = dat[-1,]

dat = dat[-which(dat$AverageCCS == "-1"), ] # remove MSDIAL results with CCS = -1
dat = dat[which(dat$MSMSspectrum != ""), ] # remove features without MS/MS
dat$SNaverage = as.numeric(dat$SNaverage)
dat = dat[which(dat$SNaverage >= minSignalToNoise), ] # remove low S/N features

dat$AlignmentID = as.numeric(dat$AlignmentID)
dat$AverageMz = as.numeric(dat$AverageMz)
dat$AverageRtmin = as.numeric(dat$AverageRtmin)
dat$AverageCCS = as.numeric(dat$AverageCCS)

#hist(log10(dat$SNaverage), 100)

#__________________________________________________
# Remove redundant peaks from spectrum keeping the most intense within tolerance (e.g., from deficient peak centroding):
RemoveRedundantPeaks = function(spectrum, minMzDistanceMSMS)
{
  spectrum = spectrum[with(spectrum, order(-ProductHeight)), ]
  spectrum$index = 1:length(spectrum[,1])
  spectrum$cluster = 0
  k = 1
  while(k <= length(spectrum[,1]))
  {
    if(spectrum$cluster[k] > 0) # continue if it is already clustered
    {
      k = k + 1
      next
    }
    
    mzx = spectrum$ProductMz[k]
    indexes = which(abs(spectrum$ProductMz - mzx) <= minMzDistanceMSMS
                    & spectrum$cluster == 0)
    spectrum$cluster[indexes] = k
    k = k + 1
  }
  spectrum = spectrum[which(spectrum$index == spectrum$cluster),]
  spectrum$cluster = NULL
  spectrum$index = NULL
  return(spectrum)
}

#__________________________________________________
# Parse fragment list:
# -------------------------------------------------
dat$MSMSspectrum = gsub(":", " ", dat$MSMSspectrum)
dat$MS1isotopicspectrum = gsub(":", " ", dat$MS1isotopicspectrum)
ms2 = NULL
for(k in 1:length(dat[,1]))
{
  print(paste("Processing feature:", k))
  tokens = unlist(strsplit(dat$MSMSspectrum[k], " "))  
  if(length(tokens) == 0)
  {
    # this is the MS1 detected without fragments
    next
  }
  # values separated by space: m/z intensity, e.g.:
  #   133.01199 2438 168.98599 281
  msx = NULL
  precHeight = unlist(strsplit(dat$MS1isotopicspectrum[k], " "))[2]
  for(j in seq(1, length(tokens), 2))
  {
    x = data.frame(AlignmentID = dat$AlignmentID[k],
                   PrecursorMz = dat$AverageMz[k],
                   PrecursorHeight = precHeight,
                   ProductMz = as.numeric(tokens[j]),
                   ProductHeight = as.numeric(tokens[j + 1]))
    msx = rbind(msx, x)
  }
  msx = RemoveRedundantPeaks(msx, minMzDistanceMSMS)
  ms2 = rbind(ms2, msx)
}
x = NULL
msx = NULL
ms2$MS1isotopes = NULL
ms2$MSMSspectrum = NULL

ms2$PrecursorHeight = as.numeric(ms2$PrecursorHeight)

# Keep best fragments: --------------- 

# remove fragments with too large intensity
ms2 = ms2[which(ms2$ProductHeight <= (ms2$PrecursorHeight * 1.3)), ]

# remove fragments with too small intensity
ms2 = ms2[which(ms2$ProductHeight >= (ms2$PrecursorHeight * 0.01)), ]

# remove fragments with close m/z to precursor m/z:
ms2 = ms2[which(abs(ms2$ProductMz - ms2$PrecursorMz) > 0.2), ]

# Compute fragment intensity rank:
ms2 = ms2[with(ms2, order(AlignmentID, -ProductHeight)), ]
ms2$intensityRank = ave(ms2$ProductHeight, ms2$AlignmentID, FUN=seq_along)

# Keep only top n fragments:
ms2 = ms2[which(ms2$intensityRank <= maxCountFragments), ]

# compute number of fragments:
ms2$CountFragments = ave(ms2$intensityRank, ms2$AlignmentID, FUN=max)

# Keep features with at least 3 fragments:
ms2 = ms2[which(ms2$CountFragments >= 3),]

dat = dat[dat$AlignmentID %in% unique(ms2$AlignmentID), ]

#__________________________________________________
# Generate decoys:
# -------------------------------------------------
decoys = NULL
mztoldecoys = 50
rtExclusionWindow = 3
ms2$hasdecoy = 0

ms2 = merge(ms2, dat[,c("AlignmentID", "AverageRtmin", "AverageCCS", "Spectrumreferencefilename")], by="AlignmentID")
ms2 = ms2[with(ms2, order(-ms2$PrecursorHeight, AlignmentID, ProductMz)), ] # sort by Precursor Height to start by the high intensity features
ms2$ProductName = ""
for(featID in unique(ms2$AlignmentID))
{
  indexesMs2 = which(ms2$AlignmentID == featID & ms2$hasdecoy == 0)
  original = ms2[indexesMs2, ]
  original$indexesMs2 =indexesMs2
  
  if(length(original$AlignmentID) < 1)
    next
  
  indexes = which(ms2$Spectrumreferencefilename == original$Spectrumreferencefilename[1]
                  & abs(ms2$PrecursorMz - original$PrecursorMz[1]) <= mztoldecoys
                  & abs(ms2$AverageRtmin - original$AverageRtmin[1]) > rtExclusionWindow
                  & ms2$CountFragments == original$CountFragments[1]
                  & ms2$AlignmentID != featID
                  & ms2$hasdecoy == 0)
  
  if(length(indexes) < 1)
    next
  
  candidates = ms2[indexes, ]
  candidates = candidates[!duplicated(candidates$AlignmentID), ]
  
  # select the candidate feature with largest rt difference:
  candidates$rtdiff = abs(candidates$AverageRtmin - original$AverageRtmin[1])
  candidates = candidates[with(candidates, order(rtdiff)), ]
  candfeatID = candidates$AlignmentID[1]
  indexesMs2 = which(ms2$AlignmentID == candfeatID)
  paired = ms2[indexesMs2, ]
  paired$indexesMs2 = indexesMs2
  #------- Flag any fragments with close mz  values:
  original$closeMz = 0
  paired$closeMz = 0
  # For each original.mz, find the index with closest smaller mz value in paired.mz:
  closestSmallerIndex = findInterval(original$ProductMz, paired$ProductMz, rightmost.closed = FALSE, all.inside = TRUE)
  indexesCloseMz = which(abs(original$ProductMz - paired$ProductMz[closestSmallerIndex]) <= 4)
  if(length(indexesCloseMz) > 0)
  {
    original$closeMz[indexesCloseMz] = 1
    paired$closeMz[closestSmallerIndex[indexesCloseMz]] = 1
  }
  # check the next element too:
  closestSmallerIndex = closestSmallerIndex + 1
  indexesCloseMz = which(abs(original$ProductMz - paired$ProductMz[closestSmallerIndex]) <= 4)
  if(length(indexesCloseMz) > 0)
  {
    original$closeMz[indexesCloseMz] = 1
    paired$closeMz[closestSmallerIndex[indexesCloseMz]] = 1
  }
  original = original[with(original, order(closeMz, intensityRank)), ]
  paired = paired[with(paired, order(closeMz, intensityRank)), ]
  # Regenerate product name to ensure rank correspondence of swapped fragments:
  original$ProductName = paste("raw", original$AlignmentID, "-f", (1:nrow(original)), sep='')
  paired$ProductName = paste("raw", paired$AlignmentID, "-f", (1:nrow(paired)), sep='')
  #--------
  swapIndexes = 1:(ceiling(length(original[,1])*0.6)) # get indexes top n most intense fragments for swapping (excluding the ones with close mz values)
  if(max(original$closeMz[swapIndexes]) == 1 || max(paired$closeMz[swapIndexes]) == 1)
  {
    print("    Features cannot be paired, they have too many fragments with close m/z values:")
    print(paste("        AlignmentID =", original$AlignmentID[1], ", paired AlignmentID =", paired$AlignmentID[1]))
    next
  }
  original$closeMz = NULL
  paired$closeMz = NULL
  
  decoy1 = rbind(paired[swapIndexes,], original[-(swapIndexes),])
  decoy1$AlignmentID = original$AlignmentID[1]
  decoy1$PrecursorHeight = original$PrecursorHeight[1]
  decoy1$PrecursorMz = original$PrecursorMz[1]
  decoy1$AverageRtmin = original$AverageRtmin[1]
  decoy1$AverageCCS = original$AverageCCS[1]
  
  decoy2 = rbind(original[swapIndexes,], paired[-swapIndexes,])
  decoy2$AlignmentID = paired$AlignmentID[1]
  decoy2$PrecursorHeight = paired$PrecursorHeight[1]
  decoy2$PrecursorMz = paired$PrecursorMz[1]
  decoy2$AverageRtmin = paired$AverageRtmin[1]
  decoy2$AverageCCS = paired$AverageCCS[1]
  
  decoys = rbind(decoys, decoy1, decoy2)
  
  #indexes = which(ms2$AlignmentID == original$AlignmentID[1] | ms2$AlignmentID == paired$AlignmentID[1])
  ms2$ProductName[original$indexesMs2] = original$ProductName
  ms2$ProductName[paired$indexesMs2] = paired$ProductName
  ms2$hasdecoy[original$indexesMs2] = 1
  ms2$hasdecoy[paired$indexesMs2] = 1
}
decoy1 = NULL
decoy2 = NULL
candidates = NULL
original = NULL
paired = NULL
decoys$indexesMs2 = NULL

# keep only features with decoys:
ms2 = ms2[which(ms2$hasdecoy == 1),]
ms2$hasdecoy = NULL
decoys$hasdecoy = NULL

decoys = decoys[which(decoys$AlignmentID %in% unique(ms2$AlignmentID)),]

# update feature to use molecule name: 
#     raw: targets
#     xxx: decoys
ms2$PrecursorName = paste("raw", ms2$AlignmentID, sep = '')
ms2$Label = "target"
decoys$PrecursorName = paste("xxx", decoys$AlignmentID, sep = '')
decoys$Label = "decoy"

ms2 = rbind(ms2, decoys)
# sort
ms2 = ms2[with(ms2, order(PrecursorName, ProductMz)), ]

# Format transition list for Skyline:
newtargets = data.frame(ms2$PrecursorName, ms2[,c("PrecursorName", "PrecursorMz", "ProductName", "ProductMz", "AverageRtmin", "AverageCCS")])

colnames(newtargets) = c("MoleculeGroup", "PrecursorName", "PrecursorMz", "ProductName", "ProductMz", "Explicit Retention Time", "Collisional Cross Section (sq A)")
newtargets$'Explicit Retention Time Window' = ExplicitRTwindow
newtargets$PrecursorCharge = -1
newtargets$ProductCharge = -1

write.table(newtargets, paste("SkylineTargetList_preliminary-target-decoy_", sampleString, "_", Sys.Date(), ".csv", sep=''),
            col.names = TRUE, row.names = FALSE, sep = ",")

newtargets$Label = ms2$Label

# Generate global statistics figures:
library(ggplot2)

pdf(pdfFileName, paper="a4r", useDingbats=FALSE) # create .pdf file
    # histogram by m/z:
    p = ggplot(newtargets[!duplicated(newtargets$PrecursorName),], aes(PrecursorMz, fill=Label))
    p = p + ylab("Number of Precursor")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))
    
    p = ggplot(newtargets, aes(ProductMz, fill=Label))
    p = p + ylab("Number of Fragments")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))
    
    # histogram by abundance:
    p = ggplot(ms2[!duplicated(ms2$PrecursorName),], aes(log10(PrecursorHeight), fill=Label))
    p = p + ylab("Number of Precursor")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))
    
    p = ggplot(ms2, aes(log10(ProductHeight), fill=Label))
    p = p + ylab("Number of Fragments")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))

    # histogram by number of fragments:
    p = ggplot(ms2[!duplicated(ms2$PrecursorName),], aes(factor(CountFragments), fill=Label))
    p = p + ylab("Number of Precursor")
    p = p + xlab("Number of Fragments per peak group")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    plot(p + geom_bar(position="dodge", alpha = 0.8, colour="black"))
        
    # histogram by RT:
    p = ggplot(newtargets[!duplicated(newtargets$PrecursorName),], aes(`Explicit Retention Time`, fill=Label))
    p = p + ylab("Number of Precursor")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))
    
    # histogram by CCS:
    p = ggplot(newtargets[!duplicated(newtargets$PrecursorName),], aes(`Collisional Cross Section (sq A)`, fill=Label))
    p = p + ylab("Number of Precursor")
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))
    
    p = ggplot(newtargets[!duplicated(newtargets$PrecursorName),], aes(`Explicit Retention Time`, PrecursorMz, colour=Label))
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    p = p + facet_wrap( ~ Label, ncol=1)
    plot(p + geom_point(alpha = 0.3, size=3))
    
    p = ggplot(newtargets[!duplicated(newtargets$PrecursorName),], aes(`Collisional Cross Section (sq A)`, PrecursorMz, colour=Label))
    p = p + theme_bw()
    p = p + theme(text=element_text(size = 15), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
    p = p + facet_wrap( ~ Label, ncol=1)
    plot(p + geom_point(alpha = 0.3, size=3))

dev.off() # close .pdf file

table(ms2[!duplicated(ms2$PrecursorName),"Label"])
