SwimR <-
function(expfile, annfile, projectname, outputPath,color="red/green", 
data.collection.interval=0.067, window.size=150, mads=4.4478, 
quantile=0.95, interval=20, degree=0.2, paralysis.interval=20, 
paralysis.degree=0.2, rev.degree=0.5){
	


	
resultAll <- list()	

if(length(which(color %in% c("red/green", "red/blue", "yellow/blue", "white/black"))) == 0){
	
    stop("The inputted color is invalid!")
}
	
#output file names
heatmap1 <- 
    file.path(outputPath, paste(projectname, "_heatmap_withingroup_ordered_globalcentering.jpg", sep=""))
	
heatmap1_txt <- 
    file.path(outputPath, paste(projectname, "_heatmap_withingroup_ordered.txt", sep=""))
	
scatterplot <-	file.path(outputPath, paste(projectname, "_scatter.jpg", sep=""))

scatterplot.nooutliers.smoothed <- 
	file.path(outputPath, paste(projectname, "_nooutliers_smoothed_scatter.jpg", sep=""))
	
histogram.nooutliers.smoothed <- 
	file.path(outputPath, paste(projectname, "_histogram.nooutliers.smoothed.jpg", sep=""))
	
hist.data.file.prefix <- 
	file.path(outputPath, paste(projectname, "_histogram.nooutliers.smoothed.data", sep=""))

sample.t_half.file <- 
	file.path(outputPath, paste(projectname, "_sample_t_half.txt", sep=""))
	
group.data.file <- 
	file.path(outputPath, paste(projectname, "_group_data.txt", sep=""))
	
individual.data.file <- 
	file.path(outputPath, paste(projectname, "_individual_data.txt", sep=""))
	
individual.data.file1 <- 
	file.path(outputPath, paste(projectname, "_individual_data1.txt", sep=""))
	
group.means.file <- 
	file.path(outputPath, paste(projectname, "_nooutliers_smoothed_scatter_data.txt", sep=""))
	
intermediate.results <- 
	file.path(outputPath, paste(projectname, "_intermediate.results.txt", sep=""))
	
	
#Subfunctions
#function for removing outliers based on MAD (median absolute deviation)
#x is a vector
#y is the number of MADs away from median. 
#MAD is a robust measure of the variability. 1 standard deviation = 1.4826 MAD.	
remove.outliers <- function(x, k){
    good <- which(abs((x-median(x))/mad(x))<k)
    return (good)
}
	
#function for calculating moving average
#x is a vector
#k is the moving window size
moving.average <- function(x, k) { 
    n <- length(x) 
    y <- rep(0, n) 
    for (i in c(k:n)){ 
        y[i] <- mean(x[(i-k+1):i])
    }
    return(y[k:n])
} 
	
#function for calculating t_half, the index when the value first goes above mean
#x is a vector
t_half_max <- function(x) {
    max<-max(x) #should this be max?
    t_half<-min(which(x>max/2))
    return(t_half)
}
	
t_half_1time <- function(x) {
    current<-0
    t_half<-0
    first.peak<-0
    for (i in 1:length(x)){
        if (x[i]<current){
            first.peak<-current
            break
        }
        else{
            current<-x[i]
        }
    }
    for (j in (i+1):length(x)){
        if (x[j]<first.peak/2){
            t_half<-j
            break
        }
    }
    if (t_half==0){
        t_half=1
    }
    return (c(i-1, t_half))
}
	
t_half <- function(x, y, z) {
    current <- 0
    t_half <- 0
    first.peak <- 0
    for (i in c(1:length(x))){
        if (x[i] < current){
            first.peak <- current
            break
        }
        else{
            current <- x[i]
        }
    }
    for (j in c(i:(length(x)-y+1))){
        if (max(x[j:(j+y-1)]) < (first.peak*z)){
            t_half <- j
            break
        }
    }
    if (t_half==0){
        t_half=1
    }
    return (c(i-1, t_half))
}
	
paralysis <- function(x, w, y, z, v) {
    start <- 0
    count <- 0
    range <- max(x)-min(x)
    threshold <- min(x)+(z*range)
    result <- c()
    result1 <- "NA"
    result[1] <- max(x)
    result[2] <- min(x)
    result[3] <- range
		
    rev_ct <- 0
    new_rev <- 1
    first_rev <- 1
    total_rev_time <- 0
    total_rev_freq <- 0
    for (i in c(w:length(x))){
        if (start==0 & x[i]<threshold){
            start <- i
        }
        if (start>0 & x[i]<threshold){
            count<-count+1
            if (count>y){
                break
            }
        }
        if (start>0 & x[i]>=threshold){
            count <- 0
            start <- 0
        }
    }
    result[4] <- as.numeric(names(x[start]))
    t_range<-as.numeric(names(x[length(x)]))-result[4]
    result[5]<-t_range
    first_rev_time <- result[4]
    if (start!=0){
        for (j in c(start:length(x))){
            if (x[j]<v*max(x) & new_rev==0){
                new_rev <- 1
                rev_end <- j
                result1 <- paste(result1, as.numeric(names(x[j])), sep="-")
                total_rev_time <- as.numeric(names(x[rev_end]))-as.numeric(names(x[rev_start]))
            }
            if (x[j]>v*max(x) & new_rev==0){
                total_rev_freq <- total_rev_freq+x[j]-(v*max(x))
                if (j == length(x)){
                    result1<-paste(result1, "-", sep="")
                }
            }
            if (x[j]>v*max(x) & new_rev==1 & first_rev==0){
                new_rev <- 0
                rev_ct <- rev_ct+1
                rev_start <- j
                total_rev_freq <- total_rev_freq+x[j]-(v*max(x))
                result1 <- paste(result1, as.numeric(names(x[j])), sep=";")
             }
             if (x[j]>v*max(x) & new_rev==1 & first_rev==1){
                new_rev <- 0
                rev_ct <- rev_ct+1
                rev_start <- j
                total_rev_freq <- total_rev_freq+x[j]-(v*max(x))
                first_rev <- 0
                first_rev_time <- as.numeric(names(x[j]))
                result1 <- first_rev_time
             }
        }
    }
		
    result[6] <- rev_ct
    result[7] <- first_rev_time-result[4]
    if (result[7]==0){
        results[7] <- result[7]/rev_ct
    }
    result[8] <- total_rev_time
    result[9] <- total_rev_time/rev_ct
    result[10] <- total_rev_freq
    return (list(r1=result, r2=result1))
}
	
	
	
require(gplots) || stop("package gplots is required!")
require(heatmap.plus) || stop("package heatmap.plus is required!")
require("R2HTML") || stop("package R2HTML is required!")
	
	
ann.row <- 1
	
cat("Processing...\n")

sink(intermediate.results, append=TRUE)
print ("parameters used")
print (paste("ann.row: ", ann.row, sep=""))
print (paste("window.size for smoothing: ", window.size,sep=""))
print (paste("mads for outlier removal: ", mads,sep=""))
print (paste("quantile for setting up the ceiling: ", quantile,sep=""))
print (paste("sustained paralytic interval size: ", interval,sep=""))
print (paste("paralytic degree: ", degree,sep=""))
print (paste("frequency range-based paralytic interval: ", paralysis.interval, sep=""))
print (paste("frequency range-based paralytic degree: ", paralysis.degree, sep=""))
print (paste("reversion degree: ", rev.degree, sep=""))
	
#read data
	
data <- read.table(expfile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
data <- as.matrix(data)
	
ann <- read.table(annfile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
ann <- rbind(ann,ann)
rownames(ann)[2] <- ""
ann <- as.matrix(ann)
	
	
#calculate sum for each sample
sum <- apply(data, 2, sum)
print ("sum for each sample")
print (sum)
	
#identify group ids
groups <- unique(ann[ann.row,])
	
	
#remove outliers for each group
#cat("Identify outlier samples...\n")
selected.samples <- vector()
for (i in 1:length(groups)){
	sumg1 <- sum[which(ann[ann.row,]==groups[i])]
#bimo <- bimodality.test(sumg1)
#bimop <- bimo@p_value
#print (paste("The p value of bimodal test for Group: ",groups[i]," is ",bimop,sep=""))
	good.samples <- names(remove.outliers(sumg1,mads))
	selected.samples <- c(selected.samples,good.samples)
}
	
print ("selected samples after removing outliers")
print (selected.samples)
print ("removed samples")
print (setdiff(colnames(data), selected.samples))
	
#get smoothed data for good samples
#cat("Get smoothed data for good samples...\n")
smoothed.data <- apply(data[,selected.samples], 2, moving.average, window.size)
rownames(smoothed.data) <- rownames(data)[window.size:nrow(data)]
	
#calculate t_half (fall below their x% first peak thrashing value for a period (interval)) 
#for each good samples and for each group
#cat("Calculate t_half (fall below their x% first peak thrashing value for a period (interval)) 
#for each good samples and for each group...\n")
t_halfs <- apply(smoothed.data, 2, t_half, round(interval/data.collection.interval), degree)
first_peak_time <- as.numeric(rownames(smoothed.data)[t_halfs[1,]])
names(first_peak_time) <- colnames(smoothed.data)
t_half_time <- as.numeric(rownames(smoothed.data)[t_halfs[2,]])
t_half_time[t_half_time==10] <- NA
names(t_half_time) <- colnames(smoothed.data)
sample_t_half <- cbind(first_peak_time, t_half_time)
write.table(sample_t_half, sample.t_half.file, sep="\t", quote=FALSE)
resultAll$sample_t_half <- sample_t_half
	
#get paralysed animals and corresponding data
paralysed.samples <- names(which(!is.na(sample_t_half[, 2])))
if(length(paralysed.samples)>0){
    paralysed.first.peak.index <- t_halfs[1, paralysed.samples]
    paralysed.data <- smoothed.data[, paralysed.samples]
    if(length(paralysed.samples)==1){
        paralysed.data <- array(paralysed.data,dim=c(length(paralysed.data), 1))
        rownames(paralysed.data) <- rownames(smoothed.data)
    }
    paralysed.data.max <- apply(paralysed.data, 2, max)
    paralysed.data.min <- apply(paralysed.data, 2, min)
		
    paralysed.data.range <- paralysed.data.max-paralysed.data.min
	
#calculate paralysis and reversion for individual animals
	
    results <- c()
    results1 <- c()
    for (m in c(1:length(paralysed.samples))){
        result <- paralysis(paralysed.data[,m], paralysed.first.peak.index[m], 
							round(paralysis.interval/data.collection.interval), paralysis.degree, rev.degree)
        results <- cbind(results, result$r1)
        results1 <- cbind(results1, result$r2)
	}
    colnames(results) <- paralysed.samples
    rownames(results) <- c("F_max", "F_min", "F_range", "T_p_start", "T_p2end", "R_count", 
							   "T_p2r", "T_r_total", "T_r_average", "R_amp")
    colnames(results1) <- paralysed.samples
    rownames(results1) <- "R_instances"
    write.table(results, individual.data.file, sep="\t", quote=FALSE)
    write.table(results1, individual.data.file1, sep="\t", quote=FALSE)
	resultAll$individual.data <- results
	resultAll$individual.data1 <- results1
}
	
#calculate group information
#cat("Calculate group information...\n")
group_data_label <- c("freq_max_mean", "freq_max_sd", "freq_min_mean", "freq_min_sd", "freq_range_mean",
					  "freq_range_sd", "paralytic_count", "non-paralytic_count", "t_half_mean",
					  "t_half_sd", "t_p_start_mean", "t_p_start_sd", "t_p2end_mean", "t_p2end_sd",
					  "rev_count", "rev_percent", "rev_frequency_mean", "rev_frequency_sd", "t_p2r_mean",
					  "t_p2r_sd", "t_r_total_mean", "t_r_total_sd", "t_r_average_mean", "t_r_average_sd",
					  "r_amp_mean", "r_amp_sd")
group_data <- matrix(NA, length(groups), length(group_data_label))
colnames(group_data) <- group_data_label
print ("paralytic and non-paralytic samples in each group")
for (i in 1:length(groups)){
    selected.columns <- names(which(ann[ann.row,selected.samples]==groups[i]))
    selected.t_half <- sample_t_half[selected.columns,2]
    non_paralytic_samples <- names(selected.t_half[which(is.na(selected.t_half))])
    paralytic_samples <- names(selected.t_half[which(!is.na(selected.t_half))])
    group.name <- paste("GROUP:", groups[i], sep="")
    print(group.name)
    print("non-paralytic samples")
    print(non_paralytic_samples)
    print("paralytic samples")
    print(paralytic_samples)
    if(length(paralytic_samples)>0){
        selected.t_half.nona <- selected.t_half[paralytic_samples]
        group_data[i, "freq_max_mean"] <- mean(results["F_max", paralytic_samples])
        group_data[i, "freq_max_sd"] <- sd(results["F_max", paralytic_samples])
        group_data[i, "freq_min_mean"] <- mean(results["F_min", paralytic_samples])
        group_data[i, "freq_min_sd"] <- sd(results["F_min", paralytic_samples])
        group_data[i, "freq_range_mean"] <- mean(results["F_range", paralytic_samples])
        group_data[i, "freq_range_sd"] <- sd(results["F_range", paralytic_samples])
        group_data[i, "paralytic_count"] <- length(paralytic_samples)
        group_data[i, "non-paralytic_count"] <- length(non_paralytic_samples)
        group_data[i, "t_half_mean"] <- mean(as.numeric(selected.t_half.nona))
        group_data[i, "t_half_sd"] <- sd(as.numeric(selected.t_half.nona))
        group_data[i, "t_p_start_mean"] <- mean(results["T_p_start", paralytic_samples])
        group_data[i, "t_p_start_sd"] <- sd(results["T_p_start", paralytic_samples])
        group_data[i, "t_p2end_mean"] <- mean(results["T_p2end", paralytic_samples])
        group_data[i, "t_p2end_sd"] <- sd(results["T_p2end", paralytic_samples])
		
        non_revertant_samples <- names(which(results["R_count", paralytic_samples]==0))
        revertant_samples <- names(which(results["R_count", paralytic_samples]>0))
        print("non-revertant samples")
        print(non_revertant_samples)
        print("revertant samples")
        print(revertant_samples)
		
        if (length(revertant_samples)>0){
            group_data[i, "rev_count"] <- length(revertant_samples)
            group_data[i, "rev_percent"] <- group_data[i, "rev_count"]/group_data[i, "paralytic_count"]
            group_data[i, "rev_frequency_mean"] <- mean(results["R_count", revertant_samples])
            group_data[i, "rev_frequency_sd"] <- sd(results["R_count", revertant_samples])
            group_data[i, "t_p2r_mean"] <- mean(results["T_p2r", revertant_samples])
            group_data[i, "t_p2r_sd"] <- sd(results["T_p2r", revertant_samples])
            group_data[i, "t_r_total_mean"] <- mean(results["T_r_total", revertant_samples])
            group_data[i, "t_r_total_sd"] <- sd(results["T_r_total", revertant_samples])
            group_data[i, "t_r_average_mean"] <- mean(results["T_r_average", revertant_samples])
            group_data[i, "t_r_average_sd"] <- sd(results["T_r_average", revertant_samples])
            group_data[i, "r_amp_mean"] <- mean(results["R_amp", revertant_samples])
            group_data[i, "r_amp_sd"] <- sd(results["R_amp", revertant_samples])
			}
        }else{
            group_data[i, "freq_max_mean"] <- NA
            group_data[i, "freq_max_sd"] <- NA
			group_data[i, "freq_min_mean"] <- NA
            group_data[i, "freq_min_sd"] <- NA
            group_data[i, "freq_range_mean"] <- NA
            group_data[i, "freq_range_sd"] <- NA
            group_data[i, "paralytic_count"] <- length(paralytic_samples)
            group_data[i, "non-paralytic_count"] <- length(non_paralytic_samples)
            group_data[i, "t_half_mean"] <- NA
            group_data[i, "t_half_sd"] <- NA
            group_data[i, "t_p_start_mean"] <- NA
            group_data[i, "t_p_start_sd"] <- NA
            group_data[i, "t_p2end_mean"] <-	NA
            group_data[i, "t_p2end_sd"] <- NA
            group_data[i, "rev_count"] <- 0
            group_data[i, "rev_percent"] <- NA
            group_data[i, "rev_frequency_mean"] <- NA
            group_data[i, "rev_frequency_sd"] <- NA
            group_data[i, "t_p2r_mean"] <- NA
            group_data[i, "t_p2r_sd"] <- NA
            group_data[i, "t_r_total_mean"] <- NA
            group_data[i, "t_r_total_sd"] <- NA
            group_data[i, "t_r_average_mean"] <- NA
            group_data[i, "t_r_average_sd"] <- NA
            group_data[i, "r_amp_mean"] <- NA
            group_data[i, "r_amp_sd"] <- NA
        }
}
sink()

rownames(group_data) <- groups
write.table(group_data, group.data.file, sep="\t", quote=FALSE)
resultAll$group_data <- group_data
	
#scatter plot for original data
#cat("Plot...\n")
jpeg(scatterplot,width=600,height=600,res=100)
par("ask"=FALSE)
for (i in c(1:length(groups))){
    selected.columns <- names(which(ann[ann.row,]==groups[i]))
    if (i==1){
        plot(as.numeric(rownames(data)), rowMeans(data[, selected.columns]), 
			 cex=0.1, ylim=c(0,2), col=i, xlab="Time", ylab="Frequency")
    } 
    else{
        points(as.numeric(rownames(data)), rowMeans(data[, selected.columns]), cex=0.1, col=i)
    }
}
legend(400, 2, groups, cex=0.9, col=c(1:6), pch=1)
dev.off()
	
	
#scatter plot for outlier-removed, smoothed data
jpeg(scatterplot.nooutliers.smoothed,width=600,height=600,res=100)
par("ask"=FALSE)
group_means <- matrix(0, nrow(smoothed.data), length(groups)*2)
colnames(group_means) <- rep("TEMP", length(groups)*2)
for (i in c(1:length(groups))){
    selected.columns <- names(which(ann[ann.row, selected.samples]==groups[i]))
    if (i==1){
        plot(as.numeric(rownames(smoothed.data)), rowMeans(smoothed.data[, selected.columns]), 
			 cex=0.1, ylim=c(0,2), col=i, xlab="Time", ylab="Frequency")
    }
    else{
        points(as.numeric(rownames(smoothed.data)), rowMeans(smoothed.data[, selected.columns]), cex=0.1, col=i)
    }
    group_means[,2*i-1] <- rowMeans(smoothed.data[,selected.columns])
    group_means[,2*i] <- apply(smoothed.data[,selected.columns],1,sd)
}
legend(400, 2, groups, cex=0.9, col=c(1:6), pch=1)
dev.off()
	
for (i in c(1:length(groups))){
    colnames(group_means)[2*i-1] <- paste(groups[i], "mean")
    colnames(group_means)[2*i] <- paste(groups[i], "stdev")
}
rownames(group_means) <- rownames(smoothed.data)
write.table(group_means, group.means.file, sep="\t", quote=FALSE)
resultAll$group_means <- group_means
	
#histogram for outlier-removed,smoothed data
jpeg(histogram.nooutliers.smoothed,width=1000,height=600,res=100)
par("ask"=FALSE)
rowN <- ceiling(length(groups)/2)
par(mfrow=c(rowN,2))

nooutliers.smoothed <- list()
for (i in c(1:length(groups))){
    selected.columns <- names(which(ann[ann.row,selected.samples]==groups[i]))
    hist.data <- hist(smoothed.data[, selected.columns], main=groups[i], ylab="Count", xlab="Frequency")
    hist.data.file <- paste(hist.data.file.prefix,groups[i],"txt",sep=".")
    hist.table <- cbind(hist.data$mids, hist.data$counts, hist.data$counts/sum(hist.data$counts))
    colnames(hist.table) <- c("mids", "counts", "density")
    write.table(hist.table, hist.data.file, sep="\t", quote=FALSE, row.names=FALSE)
    nooutliers.smoothed[[i]] <- hist.table
}
dev.off()
names(nooutliers.smoothed) <- groups
resultAll$nooutliers.smoothed.data <- nooutliers.smoothed
	
#order data within group based on t-half
ordered.samples<-vector()
for (i in 1:length(groups)){
    selected.columns <- names(which(ann[ann.row, selected.samples]==groups[i]))
    ordered.group <- selected.columns[order(as.numeric(sample_t_half[selected.columns, "t_half_time"]))]
    ordered.samples <- c(ordered.samples, ordered.group)
}
group.ordered.data <- smoothed.data[, ordered.samples]
write.table(group.ordered.data, heatmap1_txt, quote=FALSE, sep="\t")
resultAll$group.ordered.data <- group.ordered.data
	
group.ordered.data <- group.ordered.data[nrow(group.ordered.data):1,]
	
#set a ceiling for the data to make the heatmap more colorful
ceiling <- quantile(smoothed.data, quantile)
group.ordered.data <- ifelse(abs(group.ordered.data)>ceiling, ceiling*(group.ordered.data/abs(group.ordered.data)), 
							 group.ordered.data)
	
#get color matrix from annotation file
color.matrix <- ann[, selected.samples]
for (i in c(1:nrow(color.matrix))){
    row.groups <- unique(color.matrix[i,])
    if(color == "red/green" | color == "red/blue"){
        row.colors <- rainbow(length(row.groups))
    }
    if(color == "yellow/blue"){
        row.colors <- topo.colors(length(row.groups))
    }
    if(color == "white/black"){
        row.colors <- gray(seq(0, 1-(1/length(row.groups)), 1/length(row.groups)))
    }
    for (j in 1:length(row.groups)){
        color.matrix[i, color.matrix[i,]==row.groups[j]] <- row.colors[j]
    }
}
clab <- matrix(as.character(t(color.matrix)), nrow=ncol(color.matrix), ncol=nrow(color.matrix))
colnames(clab) <- rownames(color.matrix)
	
#generating column labels
col.lab <- paste(colnames(group.ordered.data), sample_t_half[ordered.samples,2], sep="_")
	
#heatmap plotting

jpeg(heatmap1, width=600, height=1300, res=100)
par("ask"=FALSE)
if(color == "red/green"){
    colp <- colorpanel(256, low="green", high="red", mid="black")
}
if(color == "red/blue"){
    colp <- colorpanel(256, low="blue", high="red", mid="white")
}
if(color == "yellow/blue"){
    colp <- colorpanel(256, low="blue", high="yellow", mid="black")
}
if(color == "white/black"){
    colp <- colorpanel(256, low="white", high="black")
}
heatmap.plus(group.ordered.data, margins=c(12,12), Rowv=NA, Colv=NA, ColSideColors=clab,
			 labRow=NA, labCol=ann[1,], cexCol=1.5, col=colp, scale="none")
dev.off()
cat("Processing completed!\n")
cat("Please see the detailed information in the outputPath directory!\n")
	
target <- HTMLInitFile(outputPath, filename="outputDescription_SwimR", BackGroundColor="#BBBBEE")
HTML("<h2>The output of SwimR</h2>", file=target)

de <- paste("<ul><li><a href=", group.data.file, ">", projectname, "_group_data.txt</a> ",
			"contains the detailed information of each group in the input data.</li>", sep="")
HTML(de, file=target)	
	
HTML("<li>Inforamation for individual animal</li>", file=target)
	
de <- paste("<a href=", sample.t_half.file, ">", projectname, "_sample_t_half.txt</a> ",
			"contains each animal and their corresponding latency to paralyze. ",
			"For non-paralyzers, N/A will be listed.<br/>", sep="")
HTML(de, file=target)
	
	
de <- paste("<a href=", individual.data.file, ">", projectname, "_individual_data.txt</a> ",
			"is a TXT file that returns reversion information for individual animals.<br/>", sep="")
HTML(de, file=target)
	

de <- paste("<a href=", individual.data.file1, ">", projectname, "_individual_data1.txt</a> ",
			"is a TXT file in which R_instances row tells the user exactly when the animal  ",
			"reverted.<br/>", sep="")
HTML(de, file=target)


HTML("<li>Heat map</li>", file=target)
	
de <- paste("<a href=", heatmap1, ">", projectname, "_heatmap_withingroup_ordered_globalcentering.pdf</a> ",
			"is a PDF of the heat map of all samples included in the data matrix after ",
			"outlier exclusion, smoothing, ordering based on the latency to paralyse, and ",
			"centering the color based on the quantile percent that can be set by the user ",
			"by changing the color parameter.<br/><br/>", sep="")
HTML(de, file=target)
	
de <- paste("<img src=", heatmap1, ">", sep="")
HTML(de, file=target)
	
de <- paste("<a href=", heatmap1_txt, ">", projectname, "_heatmap_withingroup_ordered.txt</a> ",
			"is a TXT file of the raw data used to plot the heat map after outlier exclusion, ",
			"smoothing and ordering based on the latency to paralyze.<br/>", sep="")
HTML(de, file=target)

	
HTML("<li>Histogram</li>", file=target)
de <- paste("<a href=", histogram.nooutliers.smoothed, ">", projectname, "_histogram.nooutliers.smoothed.jpg</a> ",
			"is a JPEG file of all frequency data points broken up into increasing 0.1 Hz bins ",
			"after exclusion and smoothing and then plotted as the fraction of the total as a histogram.<br/><br/>", sep="")
HTML(de, file=target)

de <- paste("<img src=", histogram.nooutliers.smoothed, ">", sep="")
HTML(de, file=target)

for(i in c(1:length(groups))){
    hist.data.file <- paste(hist.data.file.prefix,groups[i],"txt",sep=".")
    de <- paste("<a href=", hist.data.file, ">", projectname, "_histogram.nooutliers.smoothed.data.",
			groups[i], ".txt</a> contains the raw data used to plot the histogram for genotype ",
			groups[i], ".<br/>", sep="")
    HTML(de, file=target)
}
	
HTML("<li>Scatter plot</li>", file=target)
de <- paste("<a href=", scatterplot, ">", projectname, "_scatter.jpg</a> ",
			"is a JPEG image of the average frequency plotted against time after ",
			"outlier exclusion, but no smoothing.<br/><br/>", sep="")
HTML(de, file=target)

de <- paste("<img src=", scatterplot, ">", sep="")
HTML(de, file=target)

de <- paste("<a href=", scatterplot.nooutliers.smoothed, ">", projectname, "_nooutliers_smoothed_scatter.jpg</a> ",
			"is a JPEG image of the average frequency plotted against time after ",
			"outlier exclusion and smoothing.<br/><br/>", sep="")
HTML(de, file=target)
	
de <- paste("<img src=", scatterplot.nooutliers.smoothed, ">", sep="")
HTML(de, file=target)
	
de <- paste("<a href=", group.means.file, ">", projectname, "_nooutliers_smoothed_scatter_data.txt</a> ",
			"is a TXT file of the raw data used to plot the smoothed scatter plot.  ",
			"Column one is time, two is average frequency and three is standard deviation.<br/>", sep="")
HTML(de, file=target)
	
de <- paste("<li><a href=", intermediate.results, ">", projectname, "_intermediate.results.txt</a> ",
		"describes some key features of your samples after running SwimR.</li></ul>", sep="")
HTML(de, file=target)	
	
HTMLEndFile()
	
if (interactive()) browseURL(file.path(outputPath,"outputDescription_SwimR.html"))
	
return(resultAll)
}

