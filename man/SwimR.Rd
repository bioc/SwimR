\name{SwimR}
\alias{SwimR}

\title{
Analyzing and visualizing worm swimming data
}
\description{
 This function analyzes and visualizes worm swimming data returned by 
 \code{link{createFrequencyMatrix}}. It places a particular emphasis 
 on identifying paralysis and quantifying the kinetic elements of paralysis 
 during swimming. 
}

\arguments{
  \item{expfile}{
    \code{expfile} is the path of the frequency matrix returned by the 
\code{link{createFrequencyMatrix}} function.
  }
  \item{annfile}{
    \code{annfile} is the path of annotation file returned by the 
\code{link{createFrequencyMatrix}} function.
 }
  \item{projectname}{
    \code{projectname} is the name of the project.
  }
  \item{outputPath}{
    \code{outputPath} is a directory which saves the plots and files returned 
 by the function.
  }
  \item{color}{
    The function provides four colors to plot the heat map plot:"red/green", 
 "red/blue", "yellow/blue" and "white/black". The default \code{color} is "red/green".
  }
  \item{data.collection.interval}{
    \code{data.collection.interval} is the time interval between two points 
 and the default is 0.067.
  }
  \item{window.size}{
    \code{window.size} is the size of the window for the running average 
 that is calculated to smooth the data. The default is 150.
  }
  \item{mads}{
    \code{mads} is the number of median absolute deviations that a given animal 
 must deviate from the median sum of frequencies to be called an outlier. 
 The default is 4.4478.
  }
  \item{quantile}{
    \code{quantile} is the proportion of data points that are used in calculating 
 the color scheme for the heat map and the default is 0.95.
  }
  \item{interval}{
    \code{interval} is the minimum time that a given animal must lie below a threshold 
 to be called a paralyzed worm for the first calculation and the default is 20.
  }
  \item{degree}{
    \code{degree} is the paralytic degree for the first calculation and 
 the default is 0.2.
  }
  \item{paralysis.interval}{
    \code{paralysis.interval} is the same as \code{interval} but for the second 
 calculation and the default is 20.
  }
  \item{paralysis.degree}{
    \code{paralysis.degree} is the paralytic degree for the second calculation 
 and the default is 0.2.
  }
  \item{rev.degree}{
    \code{rev.degree} is the threshold that an animal must cross to be called 
 a revertant and the default is 0.5.
  }
}

\value{
  The SwimR function returns a list object which contains the following information:
  \item{sample_t_half}{
    \code{sample_t_half} contains each animal and their corresponding latency to paralyze. 
 For non-paralyzers, N/A will be listed. 
  }
  \item{group_data}{
    The columns of \code{group_data} is defined as follow.
    "freq_max_mean": Mean maximal swimming frequency;
    "freq_max_sd": Standard deviation of Mean maximal swimming frequency;
    "freq_min_mean": Mean minimum swimming frequency;
    "freq_min_sd": Standard deviation of Mean minimum swimming frequency;
    "freq_range_mean": Mean range between maximum and minimum;
    "freq_range_sd": Standard deviation of Mean range between maximum and minimum;
    "paralytic_count": The number of paralyzed animals amongst the samples;
    "non-paralytic_count": The number of non-paralyzed animals amongst the samples;
    "t_half_mean": Mean latency to cross the paralytic threshold set by the users 
 (default is 20% of Frequency range) and stay below it for the interval specified 
 interval (default is 20 seconds);
    "t_half_sd": Standard deviation of \code{t_half_mean};
    "t_p_start_mean": The mean time point (in seconds) at which each animal crosses 
 a frequency that is min+paralytic threshold and stays below that threshold for 
 the paralytic interval;
    "t_p_start_sd": Standard deviation of \code{t_p_start_mean};
    "t_p2end_mean": The average range of time after paralysis;
    "t_p2end_sd": Standard deviation of \code{t_p2end_mean};
    "rev_count": The number of revertants amongst the samples as defined by 
 the threshold set by the user (default is animals have to recross 50% of their 
 frequency range for any length of time;
    "rev_percent": The number of revertants;
    "rev_frequency_mean": The number of reversion events;
    "t_p2r_mean": Mean time between 1st reversion and \code{t_p_start_mean};
    "t_p2r_sd": Standard deviation of \code{t_p2r_mean};
    "t_r_total_mean": Mean of total time spent in reversion for all revertants;
    "t_r_total_sd": Standard deviation of \code{t_r_total_mean};
    "t_r_average_mean": Mean length of an individual reversion event;
    "t_r_average_sd": Standard deviation of \code{t_r_average_mean};
    "r_amp_mean": Mean of total amplitude of reversion for all revertants, 
 where amplitude is defined by the area beyond the reversion threshold set 
 by user (default is 50% Freq range) during reversion, calculated by summing up 
 discrete values for each measurement (same unit as frequency);
    "r_amp_sd": Standard deviation of \code{r_amp_mean}.
  }
  \item{group.ordered.data}{
    \code{group.ordered.data} contains the data after outlier exclusion, smoothing and 
 ordering based on the latency to paralyze.
  }
  \item{individual.data}{
    \code{individual.data} contains reversion information for individual animals. 
 The definitions are identical to the \code{group_data}, but "R_count" is the number of 
 reversion events for that animal. If there is no paralyzed animal, it will not be returned.
  }
  \item{individual.data1}{
    For animals that paralyzed: The R_instances row tells the user exactly 
 when the animal reverted. For animals that did not revert, N/A will be listed. 
 If there is no paralyzed animal, it will not be returned.
  }
  \item{group_means}{
    \code{group_means} contains average frequency and standard deviation for each group.
 The row names are the time.
  }
  \item{nooutliers.smoothed.data}{
    \code{nooutliers.smoothed.data} is a list object which contains all frequency data 
 points broken up into increasing 0.1 Hz bins after exclusion and smoothing for each of genotypes.
  }
}

\details{
  The SwimR function outputs 13 files:
  1. \code{output_SwimR.html} contains a summary of all output files.

  2. \code{P_sample_t_half.txt} is a TXT file which contains the information of \code{sample_t_half}. 
 "P" of "P_sample_t_half.txt" is the \code{projectname} inputted by users.

  3. \code{P_group_data.txt} is a TXT file which contains the information of \code{group_data}.

  4. \code{P_heatmap_withingroup_ordered_globalcentering.jpg} is a JPEG file of the 
 heat map of all of the samples included in the data matrix after outlier exclusion, 
 smoothing, ordering based on the latency to paralyse, and centering the color 
 based on the quantile percent that can be set by the user in the parameters section 
 of SwimR.

  5. \code{P_heatmap_withingroup_ordered.txt} is a TXT file of the raw data used to plot the heat map 
 which is the same with \code{group.ordered.data}.

  6. After exclusion and smoothing, \code{P_histogram.nooutliers.smoothed.jpg} 
 is a JPEG file of all frequency data points broken up into increasing 0.1 Hz 
 bins and then plotted as the fraction of the total as a histogram.

  7. \code{P_histogram.nooutliers.smoothed.data.G.txt} is a TXT file of 
 the raw data used to plot the histogram, which is the same with \code{nooutliers.smoothed.data}.
  "G" in the "P_histogram.nooutliers.smoothed.data.G.txt" is the genotype in the annotation file.

  8. \code{P_individual_data.txt} is a TXT file which contains the information of \code{individual.data}.
 If there is no paralyzed animal, this file will not be outputted.

  9. \code{P_individual_data1.txt} is a TXT file which contains the information of \code{individual.data1}. 
 If there is no paralyzed animal, this file will not be outputted.

  10. \code{P_intermediate.results.txt} describes some key features of your 
 samples after running SwimR, and is a great way to get a quick look at the 
 incidence of paralysis amongst your samples. At the top of the file, it lists 
 the parameters used in the subsequent calculations. Below that, it lists the 
 summed frequency values for each of the animals included in the sample. And 
 then the p value of the bimodal test for each genotype was listed. Below that, 
 it lists each of the animals included and excluded after outlier detection. 
 After that, it lists which animals were considered paralyzed and which not. 
 For paralyzed animals, it then lists which of them were called revertants.

  11. \code{P_scatter.jpg} is a JPEG image of the average frequency plotted 
 against time after outlier exclusion, but w/o smoothing.
 
  12. \code{P_nooutliers_smoothed_scatter.jpg} is a JPEG image of the average 
 frequency plotted against time after outlier exclusion and smoothing.
 
  13. \code{P_nooutliers_smoothed_scatter_data.txt} is a TXT file of the raw data 
 used to plot the smoothed scatter.pdf, which is the same with \code{group_means}.
}

\author{
    Jing Wang, Andrew Hardaway and Bing Zhang
}

\examples{
    expfile <- system.file("extdata", "SwimExample", "SwimR_Matrix.txt", package="SwimR")
    annfile <- system.file("extdata", "SwimExample", "SwimR_anno.txt", package="SwimR")
    projectname <- "SwimR"
    outputPath <- getwd()
    result <- SwimR(expfile, annfile, projectname, outputPath, color = "red/green",
 data.collection.interval = 0.067, window.size = 150, mads = 4.4478, quantile = 0.95, 
 interval = 20, degree = 0.2, paralysis.interval = 20, paralysis.degree = 0.2, 
 rev.degree = 0.5)
}

\seealso{
    \code{\link{createFrequencyMatrix}}
}

\keyword{ methods }