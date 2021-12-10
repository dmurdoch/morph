\name{morph3d}
\alias{morph3d}

\title{
Performs a morphological segmentation of a binary 3D array
}

\description{
This is the extension of 2D morphological segmentation of binary pixels in a planar raster setting to the 3D domain of voxels. A 3D binary array with 1 (belonging to the feature of interest) and 0 (not belonging to the feature of interest) is the only input required for the mathematical morphological segmentation to proceed. The result is a list object containing the original input data, the categorical segmentation layer, and various summary outputs that may be necessary at later interpretation stages. The algorithm assigns voxels into mutually exclusive categories: outside, mass, skin, crumb, circuit, antenna, bond, void-volume, and void based on their relative positioning and connectivity within the feature of interest.
}

\usage{
morph3d(DATACUBE = NULL, VERBOSE = FALSE, PLOT = FALSE, FINALPLOT = TRUE, PLOTIDS = FALSE)
}

\arguments{
  \item{DATACUBE}{
This is the input 3D data array. The dimensions do not need to be equal, but the data must be binary (0,1) and be of numeric type.
}

  \item{VERBOSE}{
This is a Boolean flag that turns on | off additional feedback as the code is running. This is generally used for testing and debugging and can easily be set to FALSE with no concern over how the code runs.
}

  \item{PLOT}{
This is a Boolean flag that turns on | off 3D interactive environment plots for each mutually exclusive morphological class. If these plots are not explicitly required, set this to FALSE, as the plotting is the bottleneck that dramatically slows processing.
}

  \item{FINALPLOT}{
This is a Boolean flag that turns on | off the 3D interactive environment plot for the final voxel segmentation result. If this is set to FALSE, it can readily be plotted in the future, since the results are saved by the function to a compound list object and the plotting function is supplied with this package.
}
  \item{PLOTIDS}{
This is a Boolean flag that turns on | off the plotting of unique voxel IDS at the centroids of each voxel during plotting. These can be very useful for identifying specific voxels and associations, but can also clutter the view if they are not expressly required.
}
}

\details{
Each voxel can belong to only a single 3D morphological class. The concept expands the planar morphological segmentation avaliable in the MSPA software package and establishes corresponding terminology to distinguish the 3D case. The processessing logic is presented and a fully-functional environment allows any binary 3D data to be processed with this engine.}

\value{
The function returns a LIST object with the following items:
\item{OriginalData }{This is the original input 3D data array of 0,1 values.}
\item{Graph }{This is the network graph object in igraph format for the feature of interest.}
\item{VoxelIDs }{This is a 3D array with the same dimensions as the input data and contains the unique voxel ID value for each voxel.}
\item{ObjectID }{If there are multiple disconnected components for the input feature of interest, each component is assigned a unique object ID. This 3D array with the same dimensions as the input data contains those unique object ID values for each voxel.}
\item{Morphology }{This is a 3D array with the same dimensions as the input data and is a categorical voxel map identifying which morphological class each voxel belongs to. The codes are: 1=outside, 2=mass, 3=skin, 4=crumb, 5=circuit, 6=antenna, 7=bond, 8=void-volume, and 9=void.}
\item{Cores }{This is a 3D array with the same dimensions as the input data and contains the codes for unique cores (of mass voxels) that may exist within the feature of interest.}
\item{ExpCores }{This is the same as Cores, but expanded outward to enlarge the volume represented by each core. This is used for internal processing but may be of use or benefit to some, so it is also provided here.}
\item{Summary }{This is a tabular summary of the morphological segmentation by class. The call name, number of voxels, and the percentage of the total volume is summarized for the actual data provided in the Morphology item.}
\item{Egg }{This is a temporary output that is used for checking and debugging. It is likely to be eliminated in the future.}
\item{Bgrnd }{This 3D array identifies the voxels outside the feature of interest, but that are not part of the voids internal to the structure itself.}
\item{VOIDvolume }{This 3D array identifies the voxels that form the holes, gaps, or voids within the structure.}
}

\references{
Remmel, T.K. (Accepted December 2021). Extending morphological pattern segmentation to 3D voxels. Landscape Ecology.
}

\author{
Tarmo K. Remmel
}

\note{
Current implementation consideres voxels as neighbours if they share a face only (orthogonal neighbours). Expansion to include diagonal neighbours, where voxels share an edge or a point, are being considered. The use of rgl provides functionality to produce interactive 3D plots that greatly facilitate visualization, but dramatically slow the processing as input volumes increase in size. Avoid plotting unless explicitly required. You are able to plot from the stored list object at anytime afterwards too, it is not necessary to plot during the segmentation call.
}

\seealso{
See Also \code{\link{morph3dplot}}, \code{\link{morph3dprep}} 
}

\examples{
# SET FINALPLOT=TRUE IF THE INTERACTIVE 3D PLOT IS DESIRED
# THAT CONTAINS THE FINAL SEGMENTATION CATEGORIES
# IF A PLOT FOR EACH CATEGORY IS REQUIRED, SET PLOT=TRUE 
LEdemoOut <- morph3d(LEdemo, PLOT=FALSE, FINALPLOT=FALSE)
str(LEdemoOut)
}

\keyword{ manip }
\keyword{ array }
