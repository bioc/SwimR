createFrequencyMatrix <-
function(inputPath, outputPath, method="Extrema", Threshold=0.6, DeltaPeakDt=1.6,
MinFrameBtwnMax=4, MinDelta=2.5, longPeriod=5, AvWindowSize=10, fps=15,
ZP_Length=100, WindowSize=30, MaxCompWin=2, minTime=0, maxTime=600){
	
	
peakdet <- function(v, delta, x){
		
maxtab <- vector()
mintab <- vector()
		
v <- as.vector(v)
nargin <- length(as.list(match.call()))-1
		
if(nargin < 3){
    x <- c(1:length(v))
}else{
    x <- as.vector(x)
    if(length(v)!=length(x)){
        stop("Input vectors v and x must have same length!")
    }
}
if(length(delta)>1){
    stop("Input argument DELTA must be a scalar!")
}else{	
    if(delta <= 0){
        stop("Input argument DELTA must be positive!")
    }
}
		
mn <- Inf
mx <- -Inf
mnpos <- NaN
mxpos <- NaN
		
lookformax <- 1
for(i in c(1:length(v))){
    this <- v[i]
    if(this > mx){
        mx <- this
        mxpos <- x[i]
    }
    if(this < mn){
        mn <- this
        mnpos <- x[i]
    }
			
    if(lookformax){
        if(this < (mx-delta)){
            maxtab <- rbind(maxtab, c(mxpos, mx))
            mn <- this
            mnpos <- x[i]
            lookformax <- 0
        }
    }else{
        if(this > (mn+delta)){
            mintab <- rbind(mintab, c(mnpos, mn))
            mx <- this
            mxpos <- x[i]
            lookformax = 1
        }
    }
}
tabresult <- list(maxtab=maxtab, mintab=mintab)
return(tabresult)
}

GetPeaks <- function(TimeSer, fps){
# Identify peaks and troughs in the time series.  The first and last points are considered peaks or troughs.
# Sometimes the peaks can be flat for a while, so we have to deal with that nasty situation.
Peaks <- TimeSer
PeakTime <- TimeSer
NumPeaks <- 1
LastDiff <- 1
TSlen <- nrow(TimeSer)
time_IDX <- c(1:max(nrow(TimeSer),ncol(TimeSer)))
		
for(Pt in c(2:(TSlen-1))){
    if(TimeSer[Pt,1] == TimeSer[(Pt+1), 1]){
        next
    }
    else{
        if((sign(TimeSer[Pt, 1] - TimeSer[LastDiff, 1]) + sign(TimeSer[(Pt+1), 1] - TimeSer[Pt, 1])) == 0){
            NumPeaks <- NumPeaks + 1
            Peaks[NumPeaks, 1] <- TimeSer[Pt, 1]
            PeakTime[NumPeaks, 1] <- time_IDX[Pt]
        }
    }
			
    LastDiff <- Pt
}
Peaks <- c(Peaks[c(1:NumPeaks), 1], TimeSer[TSlen, 1])
PeakTime <- c(PeakTime[c(1:NumPeaks), 1],time_IDX[TSlen])
Peaks <- array(Peaks,dim=c(length(Peaks), 1))
PeakTime <- array(PeakTime,dim=c(length(PeakTime), 1))
NumPeaks <- NumPeaks + 1
res <- list(Peaks=Peaks,PeakTime=PeakTime, NumPeaks=NumPeaks)
return(res)
}
				   
RTfilt <- function(Peaks, PeakTime,Thresh){
# Use the racetrack-filter algorithm to eliminate small cycles from the list of peaks .
# Initialize the algorithm.
				   
NumPeaks <- nrow(Peaks)
FiltPeaks <- array(0,dim=c(NumPeaks, 1))
FiltPeakTime <- array(0,dim=c(NumPeaks, 1))
S1 <- Peaks[1, 1]
S2 <- Peaks[2, 1]
S1_TIME <- PeakTime[1, 1]
S2_TIME <- PeakTime[2, 1]
		
Ind <- 2
NumKept <- 0
while(Ind < NumPeaks){
    Ind <- Ind + 1
    S3 <- Peaks[Ind, 1]
    S3_TIME <- PeakTime[Ind, 1]
			
# Calculate the absolute difference between data points.
				   
    Diff12 <- abs(S1 - S2)
    Diff23 <- abs(S2 - S3)
    Diff31 <- abs(S3 - S1)
			
# Must find the greatest difference.  If difference is greater then the threshold
# value, continue with the algorithm; otherwise go back and get another data point.
				   
    if (Diff12 < Thresh){
# See if Diff23 is greatest.
				   
        if((Diff23 >= Diff12) && (Diff23 >= Diff31)){
            S1 <- S2
            S2 <- S3
            S1_TIME <- S2_TIME
            S2_TIME <- S3_TIME
						
            if(Diff23 >= Thresh)
                break
            }else{
                if((Diff31 >= Diff12) && (Diff31 >= Diff23)){
                    S2_TIME <- S3_TIME
                    if(Diff31 >= Thresh)
				        break
                }
            }
		}else{
			Ind <- Ind - 1
			break
	 }
}
				   
# If Diff32 is greater then the threshold, write s1 to output file
# and move s1, s2, and s3 forward in the data file by one.
# If Diff32 is less than the threshold, points must be discarded.
		
while(Ind < NumPeaks){
    Ind <- Ind + 1
    S3 <- Peaks[Ind, 1]
    S3_TIME <- PeakTime[Ind, 1]
				
    Diff12 <- abs(S1-S2)
    Diff23 <- abs(S2-S3)
				
    if(Diff23 >= Thresh){
        NumKept <- NumKept + 1
        FiltPeaks[NumKept, 1] <- S1
        FiltPeakTime[NumKept, 1] <- S1_TIME
				   
        S1 <- S2
        S2 <- S3
        S1_TIME <- S2_TIME
        S2_TIME <- S3_TIME
    }else{
        Ind <- Ind + 1
        if(Ind > NumPeaks){
            break
        }
        S3 <- Peaks[Ind,1]
        S3_TIME <- PeakTime[Ind,1]
        if(abs(S1-S3) > Diff12){
            S2 = S3
        }
    }
}
			
FiltPeaks[NumKept+1, 1] <- S1
FiltPeaks[NumKept+2, 1] <- S2
FiltPeakTime[NumKept+1, 1] <- S1_TIME
FiltPeakTime[NumKept+2, 1] <- S2_TIME
						
# Eliminate the unused elements.
		    
FiltPeaks <- FiltPeaks[1:(NumKept+2),]
FiltPeakTime <- FiltPeakTime[1:(NumKept+2),]
FiltPeaks <- array(FiltPeaks,dim=c(length(FiltPeaks), 1))
FiltPeakTime <- array(FiltPeakTime,dim=c(length(FiltPeakTime), 1))
				   
re <- list(FiltPeaks=FiltPeaks, FiltPeakTime=FiltPeakTime)
return(re)
}
				   
				   
				   
FrequencyAnalysis <- function(filename, inputPath, outputPath, method, 
							  Threshold, DeltaPeakDt, MinFrameBtwnMax, MinDelta, longPeriod, 
							  AvWindowSize, fps, ZP_Length, WindowSize, MaxCompWin){
		
		
filedir <- file.path(inputPath, filename)
Data <- read.csv(filedir, header=FALSE,stringsAsFactors=FALSE)
		
#Parse
X <- array(Data[, 1], dim=c(nrow(Data), 1))
Y <- array(Data[, 2], dim=c(nrow(Data), 1))
Hip <- array(Data[, 3], dim=c(nrow(Data), 1))
Neck <- array(Data[, 4], dim=c(nrow(Data), 1))
Knee <- array(Data[, 5], dim=c(nrow(Data), 1))
		
#Calculate some related information
length <- nrow(X)
		
#Calculate time
t <- seq(1/fps, length/fps, 1/fps)
		
##FFT Analysis PreSteps
		
#Display to outpute
#cat("Processing...\n")
		
N <- WindowSize + ZP_Length
freq <- c(0:(N-1))*fps/N
		
#Low-pass filter raw data at .5 * Nyquist = 3.75 Hz Cutoff
fir1_result <- fir1(3, 0.5)
fir1_result <- as.Arma(fir1_result)
b <- fir1_result$b
a <- fir1_result$a
HipLP <- filter(b, a, Hip)
		
##Calculate frequencies for Hip using FFT
		
#Display to output
#cat("Calculating Hip Frequencies...\n")
		
#Initialize MaxHipFFT
MaxHipFFT <- array(0, dim=c(1, length))
for(i in c((WindowSize/2):(length(HipLP)-(WindowSize/2)))){
    range <- c((i-(WindowSize/2-1)):(i+(WindowSize/2)))
    HipLP_FFT <- abs(fft(c((HipLP[range]-mean(HipLP[range]))*hamming(WindowSize),rep(0, ZP_Length))))
    HipLP_FFT <- signif(HipLP_FFT, digits=8)
    Fidx <- which(HipLP_FFT==max(HipLP_FFT))
#if(length(Fidx)==1 && Fidx!=1)
#cat("I:",i,"Fidx:",length(Fidx),"\n")
    Fidx <- Fidx[1]
    MaxHipFFT[1, i] <- freq[Fidx]
}
		
#Fill in what is missing
MaxHipFFT[1, 1:(WindowSize/2-1)] <- MaxHipFFT[1, WindowSize/2]
MaxHipFFT[1, i:length(HipLP)] <- MaxHipFFT[1, i-1]
		
#Low pass filter the result
fir1_result <- fir1(75, 0.01)
fir1_result <- as.Arma(fir1_result)
b <- fir1_result$b
a <- fir1_result$a
MaxHipFFT_LP <- filter(b, a, MaxHipFFT)
		
#Move very long period cycles to zero, update 2013/01/10
idx <- which(MaxHipFFT_LP < 1/longPeriod)
MaxHipFFT_LP[idx] <- 0
		
## Find frequencies through ExtremaCount
rm(freq)
	
#this is a flag which ensures a max cannot be found until 
#after the hip angle has gone under pi at least once since the last max
ReadyForMax <- 1  
		
		
#initialize variables
count <- 1
MaxLoc <- array(1, dim=c(1, (length-MaxCompWin)))
MaxFreq <- array(1, dim=c(1, (length-MaxCompWin)))
		
#loop through data, find maxima and count frames between to maxima to determine frequency
		
for( i in c((1+MaxCompWin):(length-MaxCompWin))){
    if(ReadyForMax == 1 && HipLP[i] > (3.14+Threshold) && MaxLoc[1,count] < (i-MinFrameBtwnMax)){
        mTemp <- max(HipLP[(i-MaxCompWin):(i+MaxCompWin)])
        comp <- HipLP[i]
        if(mTemp==comp){
            count <- count+1
            MaxLoc[1, count] <- i
            Per <- t[MaxLoc[1, count]] - t[MaxLoc[1, (count-1)]]
            MaxFreq[1,count] <- 1/Per
            ReadyForMax <- 0
        }
    }else{
        if(HipLP[i] < 3.14)
            ReadyForMax <- 1
    }
}
		
# Construct a filled in version of the array
		
freqFull <- array(0, dim=c(1, length))
for( i in c(1:(count-1))){
    freqFull[1, (MaxLoc[1, i]:(MaxLoc[1, i+1]-1))] <- MaxFreq[1, i+1]
}
		
freqFull[1, (MaxLoc[1, count]:length)] <- MaxFreq[1, count]
		
#Move very long period cycles to zero, update 01/10/2013
idx <- which(freqFull < 1/longPeriod)
freqFull[idx] <- 0
		
#Create SmoothedVersion
fir1_result <- fir1(60, 0.01)
fir1_result <- as.Arma(fir1_result)
b <- fir1_result$b
a <- fir1_result$a
SmFF <- filter(b, a, freqFull)
		
		
##Peak det method
tabresult <- peakdet(HipLP, DeltaPeakDt)
maxtab <- tabresult$maxtab
mintab <- tabresult$mintab
		
maxLength <- nrow(maxtab)
MaxFreqTab <- array(0, dim=c(1, maxLength))
for(i in c(2:maxLength)){
    Per <- t[maxtab[i, 1]] - t[maxtab[(i-1), 1]]
    MaxFreqTab[1, i] <- 1/Per
}
		
#Construct a filled in version of the array
		
freqFullTab <- array(0, dim=c(1, length))
for(i in c(1:(maxLength-1))){
    freqFullTab[1, (maxtab[i, 1]:(maxtab[(i+1), 1]-1))] <- MaxFreqTab[1, i+1]
}
		
freqFullTab[1, (maxtab[maxLength, 1]:length)] <- MaxFreqTab[1, maxLength]
		
#Move very long period cycles to zero, update 01/10/2013
idx <- which(freqFullTab < 1/longPeriod)
freqFullTab[idx] <- 0
		
#Create SmoothedVersion
fir1_result <- fir1(60, 0.01)
fir1_result <- as.Arma(fir1_result)
b <- fir1_result$b
a <- fir1_result$a
SmFFTab <- filter(b, a, freqFullTab)
		
#Racetrack + getPeaks Method (RT_GP) PF, update 01/10/2013
		
re <- GetPeaks(Hip, fps)
Peaks <- re$Peaks
PeakTime <- re$PeakTime
NumPeaks <- re$NumPeaks
		
re <- RTfilt(Peaks, PeakTime, MinDelta)
FiltPeaks <- re$FiltPeaks
FiltPeakTime <- re$FiltPeakTime
		
NumPeaks <- nrow(FiltPeaks)
		
FiltPeakTime <- round(FiltPeakTime)
FiltPeakTime[FiltPeakTime[, 1]<1, 1] <- 1
		
				   
maxtab <- FiltPeakTime
		
		
maxLength <- nrow(FiltPeakTime)
				   
MaxFreqTab <- array(0, dim=c(1, maxLength))
for(i in c(2:maxLength)){
    Per <- t[maxtab[i, 1]] - t[maxtab[(i-1), 1]]
    MaxFreqTab[1, i] <- 0.5/Per
}
				   
freqFullTab <- array(0,dim=c(1, length))
for(i in c(1:(maxLength-1))){
    freqFullTab[1, (maxtab[i, 1]:(maxtab[(i+1), 1]-1))] <- MaxFreqTab[1, i+1]
}
				   
freqFullTab[1, (maxtab[maxLength,1]:length)] <- MaxFreqTab[1, maxLength]
		
idx <- which(freqFullTab < 1/longPeriod)
freqFullTab[idx] <- 0
		
#Create SmoothedVersion
fir1_result <- fir1(60, 0.01)
fir1_result <- as.Arma(fir1_result)
b <- fir1_result$b
a <- fir1_result$a
SmFF_RT_GP <- filter(b, a, freqFullTab)
		
## Compute time averaged versions of each array (FFT, ExCnt, PeakDet)
		
#AvCount <- floor(length(t)/(AvWindowSize*fps))
#AvLength <- AvWindowSize*fps
		
#FFT_AV <- vector()
#EC_AV <- vector()
#PD_AV <- vector()
#t_AV <- vector()
#for(i in c(1:AvCount)){
    #idx <- c(((i-1)*AvLength+1):(i*AvLength))
    #FFT_AV <- c(FFT_AV,mean(MaxHipFFT_LP[idx]))
    #EC_AV <- c(EC_AV,mean(SmFF[idx]))
    #PD_AV <- c(PD_AV,mean(SmFFTab[idx]))
    #t_AV <- c(t_AV,mean(t[idx]))
#}
		
## Plot Frequencies on one chart
#Display to output
#cat("Plotting...\n")
		
FigureFileName <- substr(filename, 1, nchar(filename)-4)
FigureOutPutDir <- file.path(outputPath, "FreqOutputFiles", FigureFileName,sep="")
dir.create(FigureOutPutDir, showWarnings=F)
		
FigureFileFulName <- paste(FigureOutPutDir, FigureFileName, "Fig.jpg",sep="")

jpeg(FigureFileFulName, height=800, width=800, res=100)
par("ask"=FALSE)
y_m <- max(c(MaxHipFFT_LP,SmFF, SmFFTab, SmFF_RT_GP)) + 0.2
plot(t/60, MaxHipFFT_LP, col="red", type="l", xlab="Time (min)", ylab="Frequency (Hz)", main="Hip Frequencies", ylim=c(0,y_m))
lines(t/60, SmFF, col="blue")
lines(t/60, SmFFTab, col="green")
lines(t/60, SmFF_RT_GP, col="magenta")

legend("topright", "(x,y)", legend=c("Fast Fourier Transform", "Extrema Count", "Peak Delta", "RT+GP"),
	   col=c("red", "blue", "green", "magenta"), lty=c(1, 1, 1, 1), lwd=c(1, 1, 1, 1),
	   pch=c(NA_integer_, NA_integer_, NA_integer_, NA_integer_))
dev.off()
		
## Plot Frequencies on subplots
#Display to output
#cat("Plotting Subplots...\n")
		
FigureFileSubName <- paste(FigureOutPutDir, FigureFileName, "FigSub.jpg",sep="")
jpeg(FigureFileSubName, height=2000, width=800, res=200)
par("ask"=FALSE)
layout(matrix(c(1, 2, 3, 4), 4, 1, byrow=F))
plot(t/60, MaxHipFFT_LP, col="red", lwd=1, type="l", xlab="Time (min)", ylab="Frequency (Hz)", main="Fast Fourier Transform")
#lines(t_AV/60,FFT_AV,col="red",lwd=4,pch=8,type="b")
grid()
		
plot(t/60, SmFF, col="blue", lwd=1, type="l", xlab="Time (min)", ylab="Frequency (Hz)", main="Extrema Count")
#lines(t_AV/60,EC_AV,col="blue",lwd=4,pch=8,type="b")
grid()
		
plot(t/60, SmFFTab, col="green", lwd=1, type="l", xlab="Time (min)", ylab="Frequency (Hz)", main="Peak Delta")
#lines(t_AV/60,PD_AV,col="green",lwd=4,pch=8,type="b")
grid()
				   
plot(t/60, SmFF_RT_GP, col="magenta", lwd=1, type="l", xlab="Time (min)", ylab="Frequency (Hz)", main="RT+GP")
grid()
		
dev.off()
		
## Write result to csvfile OriginalFilename + Freq
#cat("Writing to file...\n")
FreqFilename <- paste(FigureOutPutDir, FigureFileName, "Freq.csv",sep="")
OutputMatrix <- cbind(as.vector(MaxHipFFT_LP), as.vector(SmFF), as.vector(SmFFTab), as.vector(SmFF_RT_GP), t)
colnames(OutputMatrix) <- c("Fast Fourier Transform", "Extrema Count", "Peak Delta", "RT+GP", "time")
write.csv(OutputMatrix, FreqFilename, row.names=FALSE,quote=FALSE)
		

if(method == "FFT"){
    MaxHipFFT_LP <- as.vector(MaxHipFFT_LP)
    datavector <- array(MaxHipFFT_LP, dim=c(length(MaxHipFFT_LP), 1))
}
if(method == "Extrema"){
    SmFF <- as.vector(SmFF)
    datavector <- array(SmFF, dim=c(length(SmFF), 1))
}
if(method == "PeakDet"){
    SmFFTab <- as.vector(SmFFTab)
    datavector <- array(SmFFTab, dim=c(length(SmFFTab), 1))
}
if(method == "RT+GP"){
    SmFF_RT_GP <- as.vector(SmFF_RT_GP)
    datavector <- array(SmFF_RT_GP, dim=c(length(SmFF_RT_GP), 1))
}
colnames(datavector) <- FigureFileName
datav <- data.frame(datavector=datavector, t=t)
return(datav)
}
	
	
#create frequencyMatrix

require("R2HTML") || stop("package R2HTML is required!")
require(signal) || stop("package signal is required!")

FigureOutPutDir <- file.path(outputPath, "FreqOutputFiles",sep="")
dir.create(FigureOutPutDir, showWarnings=FALSE)
	
filenames <- list.files(inputPath)
frequencyMatrix <- vector()
frequencyName <- vector()
tl <- 0
cat("Processing...\n")
tg <- 0
for(i in c(1:length(filenames))){
    filename <- filenames[i]
    if(grep("_", filename)==1){
        tg <- tg+1
        genotype <- strsplit(filename, "_")[[1]][1]
        frequencyName <- c(frequencyName, genotype)
    }else{
        cat("\nFile: ", filename, " has a wrong format of file name, please re-format the file name like genotype_***.txt\n")
    }
}
	
if(tg == length(filenames)){
    for(i in c(1:length(filenames))){
        filename <- filenames[i]
        cat("\nFile: ", filename, "is in process...\n")
        datav <- FrequencyAnalysis(filename, inputPath, outputPath, method, Threshold,
								   DeltaPeakDt, MinFrameBtwnMax, MinDelta, longPeriod,
								   AvWindowSize, fps, ZP_Length, WindowSize, MaxCompWin)
        t <- datav[,2]
        if(max(t) < maxTime){
            stop("The max time in file ", filename, " is ", max(t), " please reset the parameter maxTime!\n")
        }
		
        datavector <- datav[datav[,2]>=minTime & datav[,2]<=maxTime,1]
        frequencyMatrix <- cbind(frequencyMatrix,datavector)
    }
    t <- t[t>=minTime & t<=maxTime]
    colnames(frequencyMatrix) <- c(1:ncol(frequencyMatrix))
    rownames(frequencyMatrix) <- t
    outputFilename <- file.path(outputPath, "frequencyMatrix.txt")
    write.table(frequencyMatrix, outputFilename, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
	
	
	target <- HTMLInitFile(outputPath, filename="outputDescription_createFrequencyMatrix", 
						   BackGroundColor="#BBBBEE")
	HTML("<h2>The output of createFrequencyMatrix</h2>", file=target)
	HTML("<ul><li>The frequency matrix of worm thrashing over time</li>", file=target)
	de <- paste("<a href=", outputFilename, ">frequencyMatrix.txt</a> ",
				"combines the analysis results of the frequency of worm ",
				"thrashing over time for all Tracker files.<br/><br/>", sep="")
	HTML(de, file=target)
	
    outputFilename <- file.path(outputPath, "annotationfile.txt")
    frequencyName <- array(frequencyName, dim=c(1, length(frequencyName)))
    colnames(frequencyName) <- c(1:length(frequencyName))
    rownames(frequencyName) <- "genotype"
    write.table(frequencyName, outputFilename, row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
	
	HTML("<li>Genotype information</li>", file=target)
	de <- paste("<a href=", outputFilename, ">annotation.txt</a> ",
				"contains all genotype information extracted from file names of ",
				"all Tracker files.<br/><br/>", sep="")
	HTML(de, file=target)
	HTML("<li>Detailed information for each tracker file</li>", file=target)
	de <- paste("The function creates a folder for each tracker file " ,
				"which can be found in the folder <a href=", FigureOutPutDir,
				">", FigureOutPutDir, "</a>.<br/><br/>", "Each folder contains the following three files:", sep="")
	HTML(de, file=target)
    de <- paste('<ul><li>The image of scatter plot of one animal plotted as ',
				'"Frequency vs Time(min)" with all four counting methods.</li>', sep="")
	HTML(de, file=target)
	de <- paste("<li>The image of scatter plot in which the counting methods are ",
				"broken up into four different plots.</li>", sep="")
	HTML(de, file=target)
	de <- "<li>The CSV file of raw data orgniazed by different methods</li></ul></ul>"
	HTML(de, file=target)
	
	
	HTMLEndFile()
	
    cat("Processing completed!\n\n") 
    cat("Please see the frequency matrix \n")
    cat("and detailed information for each animal in the outputPath directory!\n")
	
		
	result <- list(frequencyMatrix=frequencyMatrix,annotation=frequencyName)
	
	if (interactive()) browseURL(file.path(outputPath,"outputDescription_createFrequencyMatrix.html"))
	
	return(result)
}else{
    stop("Please re-format the tracker file names showned above to genotype_***.txt and then run the function again!\n")
}
	
}

