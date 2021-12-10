\name{morph3dplot}

\alias{morph3dplot}

\title{
A plot function that converts a 3D array of voxels into an interactive rgl plot.
}

\description{
This plot functionality can be used for any 3D array data that represent voxAels. This functionality is called directly by morph3d() but can be called independently by passing it 3D array data that has been processed by morph3dprep().
}

\usage{
morph3dplot(data = NULL, CELLID = TRUE, LEGEND = FALSE, ORIGTRANSP = TRUE)
}

\arguments{
  \item{data}{
This is the voxel data to be plotted and needs to be processed with morph3dprep() first.
}
  \item{CELLID}{
A Boolean flag used to control whether the unique voxel IDs are plotted or not.
}
  \item{LEGEND}{
A Boolean flag that determines whether a legend will be drawn on the ouput interactive plot or not. When this is set to TRUE, plotting is much slower.
}
  \item{ORIGTRANSP}{
A Boolean flag that is used to set the transparency of original binary input voxel data. It can be used to help produce better visualizations but is not an overly important argument for most.
}
}

\details{
Output classes are coloured (with varying degrees of transparency to make classes easier to see) as follows. Mass = "green", Skin = "black", Crumb = "brown", Circuit = "orange", Antenna = "pink", Bond = "cornflowerblue", Void-volume = "navy", and Void = "seagreen".
}

\value{
The output is an interactive 3D rgl plot of the voxels supplied as input. Morphological segmentations produced by morph3d() will always have elements drawn in consistent color for ease of interpretation.
}

\references{
Remmel, T.K. 2021. (Unpuplished manuscript)
}

\author{
Tarmo K. Remmel
}

\note{
No further notes available at this time.
}

\seealso{
See Also \code{\link{morph3d}}, \code{\link{morph3dprep}} 
}

\examples{
LEdemoOut <- morph3d(DATACUBE = LEdemo, VERBOSE = FALSE, FINALPLOT = FALSE)
morph2plot <- morph3dprep(INCUBE=LEdemoOut$Morphology, ORIG=FALSE)
# PLOTTING IS TIME INTENSIVE 
\donttest{morph3dplot(data=morph2plot, CELLID=FALSE, LEGEND=FALSE, ORIGTRANSP=TRUE)}
}

\keyword{ iplot }
\keyword{ array }
