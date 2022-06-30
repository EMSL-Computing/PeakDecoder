
# This script uses the extracted precursor and fragment ion signals XIC metrics of the preliminry training set exported from Skyline.
# The training set is filtered to keep the best quality peak groups.
# A support vector machine is trained using multiple scores calculated from the XIC metrics of the training set: 
# the cosine similarity between the expected and XIC intensities, and the mean and standard deviation of each precursor and its fragments for retention time (RT), LC-FWHM and mass error metrics. These scores are used as ML features which measure co-elution and similarity to the expected values. 
# After scoring the training set, the true and false positives are used to estimate a false-discovery rate (FDR).
# The trained SVM is saved.
#
# INPUT:
#   - Skyline-Report_preliminary-target-decoy.csv
# OUTPUT:
#   - PeakDecoder.Rda (trained scoring model)
#   - PeakDecoder-TrainingOutput_[sample]_[date].txt
#   - PeakDecoder-FDR-thresholds_[sample]_[date].txt
#   - TargetDecoy-filtered-figures_[sample]_[date].pdf
#   - TargetDecoy-score-and-metrics_[sample]_[date].csv
#   - ScoringTraining-performance-figures_[sample]_[date].pdf

basePathData = "/Users/xxx/Documents" # <--- update this path!

setwd(file.path(basePathData,"data/Asper"))
inputFileTraining = "Skyline-Report_preliminary-target-decoy_Asper.csv"
minPrecursorSignalToNoise = 5 
minCountFragmentsGoodMetrics = 2

#setwd(file.path(basePathData,"data/Pput"))
#inputFileTraining = "Skyline-Report_preliminary-target-decoy_Pput.csv"
#minPrecursorSignalToNoise = 5 
#minCountFragmentsGoodMetrics = 2

#setwd(file.path(basePathData,"data/Rhodo"))
#inputFileTraining = "Skyline-Report_preliminary-target-decoy_Rhodo_reducedFile.csv"
#minPrecursorSignalToNoise = 100
#minCountFragmentsGoodMetrics = 4

sampleString = basename(getwd())
pdfFileName = paste("TargetDecoy-filtered-figures_", sampleString, "_", Sys.Date(), ".pdf", sep='')
pdfFileNameMetrics = paste("ScoringTraining-performance-figures_", sampleString, "_", Sys.Date(), ".pdf", sep='')
mztol = 15 # ppm

# Load training data:
dat = read.csv(file = inputFileTraining, sep = ',', stringsAsFactors = FALSE)
dat$Area = as.numeric(dat$Area)
dat$Area[is.na(dat$Area)] = 0
dat$Height = as.numeric(dat$Height)
dat$Height[is.na(dat$Height)] = 0
dat$Mass.Error.PPM = as.numeric(dat$Mass.Error.PPM)
dat$Mass.Error.PPM[is.na(dat$Mass.Error.PPM)] = max(dat$Mass.Error.PPM, na.rm = TRUE)
dat$Fwhm = as.numeric(dat$Fwhm)
dat$Fwhm[is.na(dat$Fwhm)] = max(dat$Fwhm, na.rm = TRUE)
dat$Background = as.numeric(dat$Background)
dat$Background[is.na(dat$Background)] = max(dat$Background, na.rm = TRUE)
dat$Retention.Time = as.numeric(dat$Retention.Time)
dat$Retention.Time[is.na(dat$Retention.Time)] = max(dat$Retention.Time, na.rm = TRUE)
dat$Collisional.Cross.Section = as.numeric(dat$Collisional.Cross.Section)

colnames(dat)[colnames(dat) == "Peptide"] = "PrecursorName"
colnames(dat)[colnames(dat) == "Fragment.Ion"] = "ProductName"

# Calculate signal to noise:
dat$SignalToNoise = (dat$Area + dat$Background) / (dat$Background + 1)

# Remove precursor as rows and include precursor values as columns:
precs = dat[which(dat$ProductName == "precursor"), ]
dat = dat[- which(dat$ProductName == "precursor"), ]
dat = merge(dat, data.frame(PrecursorName=precs$PrecursorName,
                            Replicate=precs$Replicate,
                            Precursor.RT=precs$Retention.Time, 
                            Precursor.Fwhm = precs$Fwhm, 
                            Precursor.Height = precs$Height,
                            Precursor.SignalToNoise = precs$SignalToNoise), 
                by=c("PrecursorName", "Replicate"), all.x = TRUE)
precs = NULL

# Compute RT difference to precursor: 
dat$DIA.RTdiff = dat$Precursor.RT - dat$Retention.Time

# Compute FWHM difference to precursor:
dat$DIA.FWHMdiff = dat$Precursor.Fwhm - dat$Fwhm

dat$Label = "target"
dat$Label[grep("xxx", dat$PrecursorName)] = "decoy"
dat$PrecursorName = sub("raw|xxx", "", dat$PrecursorName)
dat$PrecursorName.Replicate = paste(dat$PrecursorName, dat$Replicate, sep = '.')
dat$originalFragmentRank = sub(".+f(.+)$", "\\1", dat$ProductName)

targets = dat[which(dat$Label == "target"), ]
decoys = dat[which(dat$Label == "decoy"), ]
dat = NULL

# Remove fragments without height (height will be used as library intensity)
targets = targets[which(targets$Height > 0),]
# Remove precursors without intensity or low S/N:
targets = targets[which(targets$Precursor.Height > 0),]
targets = targets[which(targets$Precursor.SignalToNoise > minPrecursorSignalToNoise),]

# Assign library intensity and rank:
targets$LibraryIntensity = targets$Height + targets$Background
targets = targets[with(targets, order(PrecursorName.Replicate, -LibraryIntensity)), ]
targets$LibraryIntensityRank = ave(targets$LibraryIntensity, targets$PrecursorName.Replicate, FUN=seq_along)
# compute number of fragments:
targets$CountFragments = ave(targets$LibraryIntensityRank, targets$PrecursorName.Replicate, FUN=max)

# Count bad metrics per fragment to keep only high quality targets:
targets$countBadMetrics = 0
indexes = which(targets$Area <= 0)
targets$countBadMetrics[indexes] = targets$countBadMetrics[indexes] + 1

# fragments with too small intensity, flag fragment if intensity < 1% of precursor height
indexes = which(targets$Height < targets$Precursor.Height * 0.01)
targets$countBadMetrics[indexes] = targets$countBadMetrics[indexes] + 1

indexes = which(abs(targets$Mass.Error.PPM) > mztol)
targets$countBadMetrics[indexes] = targets$countBadMetrics[indexes] + 1

# RT difference to precursor, flag fragment if RT difference is larger than tolerance: 
#indexes = which(abs(targets$DIA.RTdiff) > targets$Precursor.Fwhm * 0.1)
indexes = which(abs(targets$DIA.RTdiff) > 0.1)
targets$countBadMetrics[indexes] = targets$countBadMetrics[indexes] + 1

# FWHM difference to precursor, flag fragment if FWHM is larger than the precursor FWHM: 
indexes = which(abs(targets$DIA.FWHMdiff) > targets$Precursor.Fwhm * 2)
targets$countBadMetrics[indexes] = targets$countBadMetrics[indexes] + 1

hist(targets$countBadMetrics)

# Count fragments with bad metrics:
targets = targets[with(targets, order(PrecursorName.Replicate, countBadMetrics)), ]
targets$badMetricsRank = ave(targets$countBadMetrics, targets$PrecursorName.Replicate, FUN=seq_along)

targets$countFragmentsBadMetrics = 0
targets$countFragmentsBadMetrics[which(targets$countBadMetrics > 0)] = 1
targets$countFragmentsBadMetrics = ave(targets$countFragmentsBadMetrics, targets$PrecursorName.Replicate, FUN=sum)
# keep groups that have at least 2 good fragments:
targets$countFragmentsGoodMetrics = targets$CountFragments - targets$countFragmentsBadMetrics
targets = targets[which(targets$countFragmentsGoodMetrics >= minCountFragmentsGoodMetrics), ]
# remove fragments with bad metrics and ranked more than double the rank of number of good fragments:
#targets = targets[which(targets$countBadMetrics == 0 |
#  (targets$countBadMetrics > 0 & targets$badMetricsRank <= targets$countFragmentsGoodMetrics * 2)), ] # <--- For Asper and Pput
targets = targets[which(targets$countBadMetrics <= 1),] # <--- For Rhodo
                          
targets$countBadMetrics = NULL
targets$countFragmentsGoodMetrics = NULL
targets$countFragmentsBadMetrics = NULL
targets$badMetricsRank = NULL
targets$CountFragments = NULL

# Map library info to decoys:
decoys$Replicate.ProductName = paste(decoys$Replicate, decoys$ProductName, sep='.')
targets$Replicate.ProductName = paste(targets$Replicate, targets$ProductName, sep='.')
decoys = merge(targets[,c("Replicate.ProductName", "LibraryIntensity", "LibraryIntensityRank")], 
               decoys, 
               by="Replicate.ProductName")

# Check training set to balance targets and decoys distributions by number of fragments and m/z :
targets$checked = 0
decoys$checked = 0
dat = NULL
set.seed(123) # set seed to get reproducible sequence of random numbers
for(k in 1:length(targets[,1]))
{
  if(targets$checked[k] == 1)
    next
  
  # Get decoys related to this target:
  indexesDecoys = which(decoys$PrecursorName.Replicate == targets$PrecursorName.Replicate[k])
  ids = as.numeric(unique(sub("raw(.+)\\-.+", "\\1", decoys$ProductName[indexesDecoys])))
  # Get fragments from the paired targets and decoys:
  indexesTargets = which(targets$Replicate == targets$Replicate[k] & targets$PrecursorName %in% ids)
  indexesDecoys = which(decoys$Replicate == targets$Replicate[k] & decoys$PrecursorName %in% ids)
  x = rbind(targets[indexesTargets,], decoys[indexesDecoys,])
  # Sort and count by originalFragmentRank
  x = x[with(x, order(originalFragmentRank, Label)), ]
  x$CountTargetDecoy = ave(x$originalFragmentRank, x$originalFragmentRank, FUN=seq_along)
  # Select the originalFragmentRanks found 4 times because it is consistent (2 times in targets and 2 times in decoys)
  # this will keep the same number of fragments for the 2 targets and 2 decoys
  selected = x$originalFragmentRank[which(x$CountTargetDecoy == 4)]
  if(length(selected) >= 3)
  {
    x = x[which(x$originalFragmentRank %in% selected),]
    # Re-swap some decoy transitions with their target one to increase the similarity and avoid too perfect decoys:
    if(length(selected) >= 4)
    {
      x = x[with(x, order(Label, originalFragmentRank)), ]
      sizeOfReplacement = sample(1:round(length(selected)/4), 1) # to replace 1 or more
      rankToFix = sample(selected[1:round(length(selected)/2)], size=sizeOfReplacement)
      indexes = which(x$originalFragmentRank %in% rankToFix)
      decfix = x[indexes[1:round(length(indexes)/2)],]
      tarfix = x[indexes[(round(length(indexes)/2)+1):length(indexes)],]
      x = x[-indexes,]
      decfix = tarfix
      decfix$Label = "decoy"
      decfix$Protein = paste("xxx", decfix$PrecursorName, sep = '')
      dat = rbind(dat, tarfix, decfix)
    }
    dat = rbind(dat, x)
  }
  
  targets$checked[indexesTargets] = 1
}

dat$PrecursorName = dat$Protein # to add raw|xxx indicating target|decoy
dat$PrecursorName.Replicate = paste(dat$PrecursorName, dat$Replicate, sep = '.')
table(dat$Label[!duplicated(dat$PrecursorName.Replicate)])

# Update intensity ranks:
dat = dat[with(dat, order(PrecursorName.Replicate, -LibraryIntensity)), ]
dat$LibraryIntensityRank = ave(dat$LibraryIntensity, dat$PrecursorName.Replicate, FUN=seq_along)
dat = dat[with(dat, order(PrecursorName.Replicate, -Area)), ]
dat$IntensityRank = ave(dat$Area, dat$PrecursorName.Replicate, FUN=seq_along)
dat$CountFragments = ave(dat$IntensityRank, dat$PrecursorName.Replicate, FUN=max)

# _____________________________________________________________________

# Generate global statistics figures:
library(ggplot2)
dat$PrecursorName.Replicate = paste(dat$Protein, dat$Replicate, sep = '.')

pdf(pdfFileName, paper="a4r", useDingbats=FALSE) # create .pdf file
# histogram by m/z:
p = ggplot(dat[!duplicated(dat$PrecursorName.Replicate),], aes(Precursor.Mz, fill=Label))
p = p + ylab("Number of Precursor")
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))

p = ggplot(dat, aes(Product.Mz, fill=Label))
p = p + ylab("Number of Fragments")
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))

# histogram by abundance:
p = ggplot(dat, aes(log10(LibraryIntensity), fill=Label))
p = p + ylab("Number of Fragments")
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))

# histogram by number of fragments:
p = ggplot(dat[!duplicated(dat$PrecursorName.Replicate),], aes(factor(CountFragments), fill=Label))
p = p + ylab("Number of Precursor")
p = p + xlab("Number of Fragments per peak group")
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
plot(p + geom_bar(position="dodge", alpha = 0.8, colour="black"))

# histogram by RT:
p = ggplot(dat[!duplicated(dat$PrecursorName.Replicate),], aes(Retention.Time, fill=Label))
p = p + ylab("Number of Precursor")
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))

# histogram by CCS:
p = ggplot(dat[!duplicated(dat$PrecursorName.Replicate),], aes(Collisional.Cross.Section, fill=Label))
p = p + ylab("Number of Precursor")
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))

p = ggplot(dat[!duplicated(dat$PrecursorName.Replicate),], aes(Retention.Time, Precursor.Mz, colour=Label))
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
p = p + facet_wrap( ~ Label, ncol=1)
plot(p + geom_point(alpha = 0.3, size=3))

p = ggplot(dat[!duplicated(dat$PrecursorName.Replicate),], aes(Collisional.Cross.Section, Precursor.Mz, colour=Label))
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
p = p + facet_wrap( ~ Label, ncol=1)
plot(p + geom_point(alpha = 0.3, size=3))

dev.off() # close .pdf file

# --------------------------------------------------------------------------------
# Function to compute Cosine similarity
cosineSimimilarity = function(x,y)
{
  return (x %*% y / sqrt(x%*%x * y%*%y))
}
# --------------------------------------------------------------------------------
# Compute descriptors (machine learning features): 
computeDescriptors = function(df)
{
  descpt = NULL
  for(prec in unique(df$PrecursorName.Replicate))
  {
    tb = df[which(df$PrecursorName.Replicate == prec),]
    
    # Library descriptors:
    x = tb[1,c("PrecursorName", "PrecursorName.Replicate", "Label")]
    x$DIA.cosSim = cosineSimimilarity(tb$Area, tb$LibraryIntensity)
    
    # LC peak shape descriptors:
    x$DIA.RTdiffSd = sd(tb$Precursor.RT[1] - tb$Retention.Time)
    x$DIA.RTdiffMean = mean(tb$Precursor.RT[1] - tb$Retention.Time)
    x$DIA.FWHMdiffSd = sd(tb$Precursor.Fwhm[1] - tb$Fwhm)
    x$DIA.FWHMdiffMean = mean(tb$Precursor.Fwhm[1] - tb$Fwhm)
    
    # Mass error descriptor:
    x$DIA.MassErrorSd = sd(tb$Mass.Error.PPM)
    x$DIA.MassErrorMean = mean(tb$Mass.Error.PPM)
    
    descpt = rbind(descpt, x)
  }
  return(descpt)
}

# --------------------------------------------------------------------------------
# Train SVM model:

descriptors = computeDescriptors(dat)
descriptors$DIA.cosSim[is.na(descriptors$DIA.cosSim)] = 0
descriptors$Label = factor(descriptors$Label)

#install.packages("e1071")
library("e1071")

set.seed(123) # set seed to get reproducible cross validation
model = svm(Label ~ ., data = descriptors[,-(1:2)], scale = TRUE, probability = TRUE, cross=10) 

save(model, file =  "PeakDecoder.Rda")

sink(paste("PeakDecoder-TrainingOutput_", sampleString, "_", Sys.Date(), ".txt", sep=''))
summary(model)

# Score training data:
x = descriptors[,-(1:2)]
x$label = NULL
pred = predict(model, x, probability = TRUE)

# Save number of examples per class:
print("Number of examples per class:")
table(descriptors$Label)
print("")

# Check accuracy:
print("Confusion matrix:")
table(pred, descriptors$Label)

sink() # close output file


#attr(pred, "decision.values")[1:4,]
prob = attr(pred, "probabilities")
results = cbind(descriptors, prob)
results$decoy = NULL
colnames(results)[colnames(results) == "target"] = "PeakDecoderScore"

# Save results and metrics for training set:
resultstemp = results
resultstemp$featTarget = as.numeric(sub("raw|xxx(.+)", "\\1", resultstemp$PrecursorName))
resultstemp = merge(resultstemp, dat[,c("PrecursorName.Replicate", "Total.Area", "Precursor.Mz", "Precursor.SignalToNoise", "CountFragments")], 
                by = "PrecursorName.Replicate")
resultstemp = resultstemp[!duplicated(resultstemp$PrecursorName.Replicate),]
resultstemp$Total.Area = as.numeric(resultstemp$Total.Area)
write.csv(resultstemp, file = paste("TargetDecoy-score-and-metrics_", sampleString, "_", Sys.Date(), ".csv", sep=''), row.names = FALSE)
resultstemp = NULL

# --------------------------------------------------------------------------------
# Compute FDR cosine sim:
results = results[with(results, order(-DIA.cosSim)), ] # sort by score DESC
results$FDR = 0
fp = 0
tp = 0
for(k in 1:length(results[,1]))
{
  if(results$Label[k] == "decoy")
    fp = fp + 1
  if(results$Label[k] == "target")
    tp = tp + 1
  
  if(fp > 0)
    results$FDR[k] = fp/(fp + tp)
}
results$ScoreType = "DIA.cosSim"
fdrdata = results[which(results$Label == "target"),c("FDR", "ScoreType")]

# Compute FDR combined score:
results = results[with(results, order(-PeakDecoderScore)), ] # sort by score DESC
results$FDR = 0
fp = 0
tp = 0
for(k in 1:length(results[,1]))
{
  if(results$Label[k] == "decoy")
    fp = fp + 1
  if(results$Label[k] == "target")
    tp = tp + 1
  
  if(fp > 0)
  {
    results$FDR[k] = fp/(fp + tp) # Käll et al. 2007 said that this is actually the FDR, the false positive rate compares against TN
    #results$simpleFDR[k] = fp/tp #  SimpleFDR from Käll et al. 2007  <--- This is not consistent
  }
}
results$ScoreType = "PeakDecoder"
fdrdata = rbind(fdrdata, results[which(results$Label == "target"),c("FDR", "ScoreType")])
results$ScoreType = NULL

# Save FDR thresholds:
thresholdsFDR = results[which(results$Label == "target"), c("PeakDecoderScore", "FDR")]
thresholdsFDR$FDR = round(thresholdsFDR$FDR, digits = 3)
thresholdsFDR$PeakDecoderScore = round(thresholdsFDR$PeakDecoderScore, digits = 2)
thresholdsFDR = thresholdsFDR[!duplicated(thresholdsFDR$PeakDecoderScore), ]

write.csv(thresholdsFDR, file = paste("PeakDecoder-FDR-thresholds_", sampleString, "_", Sys.Date(), ".csv", sep=''), row.names = FALSE)

# --------------------------------------------------------------------------------
# Make figures:

pdf(pdfFileNameMetrics, paper="a4r", useDingbats=FALSE) # create .pdf file

# Make FDR figure:
p = ggplot(fdrdata, aes(FDR, color=ScoreType))
p = p + ylab("True positive rate")
p = p + xlab("False discovery rate")
p = p + theme_bw()
p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
              axis.text.x = element_text(angle = 40, hjust = 1))
p = p + coord_cartesian(xlim = c(0,1))
plot(p + stat_ecdf(size=1))

# Make SVM figure:
p = ggplot(data = NULL, aes(x=1:length(model$accuracies), y=model$accuracies))
p = p + theme_bw()
p = p + labs(title=paste("SVM model final accuracy:", format(model$tot.accuracy, nsmall = 0, digits=6, scientific = FALSE)))
p = p + xlab("Cross validation fold")
p = p + ylab("Accuracy")
p = p + theme(text=element_text(size = 20), axis.title = element_text(size = 20))
p = p + scale_x_continuous(breaks = 1:10)
plot(p + geom_abline(intercept = model$tot.accuracy, slope = 0, color="red", size = 0.5)
     + geom_point(colour="black")
     + geom_text(aes(x=1:length(model$accuracies), y=model$accuracies+0.07,
                     label = format(model$accuracies, nsmall = 0, digits=6, scientific = FALSE)), 
                 color="blue", position=position_dodge(.2), hjust=.5, angle = 80))

# Make histogram figures:
for(colindex in 4:(length(colnames(results))-1))
{
  p = ggplot(results, aes_string(colnames(results)[colindex], fill="Label"))
  p = p + ylab("Number of peak groups")
  p = p + theme_bw()
  p = p + theme(text=element_text(size = 30), axis.title = element_text(size = 30),
                axis.text.x = element_text(angle = 40, hjust = 1))
  plot(p + geom_histogram(position="dodge", alpha = 0.8, colour="black", bins = 30))
}

dev.off() # close .pdf file

