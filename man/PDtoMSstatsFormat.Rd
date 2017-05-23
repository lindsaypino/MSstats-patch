\name{PDtoMSstatsFormat}
\alias{PDtoMSstatsFormat}
\title{Generate MSstats required input format for Proteome discoverer output}
\description{Convert Proteome discoverer output into the required input format for MSstats.}

\usage{
PDtoMSstatsFormat(input,
      annotation,
      useUniquePeptide=TRUE,
      summaryforMultipleRows=max,
      fewMeasurements="remove",
      removeOxidationMpeptides=FALSE,
      removeProtein_with1Peptide=FALSE)	
}

\arguments{
  \item{input}{name of Proteome discover PSM output, which is long-format. "Protein.Group.Accessions", "#Proteins", "Sequence", "Modifications", "Charge", "Intensity",  "Spectrum.File" are required.}
  \item{annotation}{name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, Run information. 'Run' will be matched with 'Spectrum.File'.}
  \item{useUniquePeptide}{TRUE(default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.}
  \item{summaryforMultipleRows}{max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
}
  \item{fewMeasurements}{'remove'(default) will remove the features that have 1 or 2 measurements across runs.}
  \item{removeOxidationMpeptides}{TRUE will remove the modified peptides including 'Oxidation (M)' in 'Modifications' column. FALSE is default.}
  \item{removeProtein_with1Peptide}{TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.}

}


\author{Meena Choi, Olga Vitek. 

Maintainer: Meena Choi (\email{mnchoi67@gmail.com})
}
