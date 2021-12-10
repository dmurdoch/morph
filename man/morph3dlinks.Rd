\name{morph3dlinks}
\alias{morph3dlinks}

\title{
An internal function that computes the spatial neighbours for voxels. It is an internal function called by morph3d() and not intended to be called on its own.
}

\description{
An internal helper function to compute voxel neighbours in orthogonal 3D space and to return the information in data frame format back to morph3d(), the function that calls this function.
}
\usage{
morph3dlinks(VOLOBJ = NULL, VOXELIDS = NULL, VERBOSE = FALSE)
}

\arguments{
  \item{VOLOBJ}{
This is a 3D array that is larger in all dimensions, within which the actual data array is inset. This ensures that outer margins will be properly processed.
}
  \item{VOXELIDS}{
This is a 3D array of the unique voxel IDS, which are used to track individual voxels and their connectivity.
}
  \item{VERBOSE}{
This is a Boolean flag for turning on | off additional feedback during processing.
}
}

\details{
This internal function simply makes the total length of the morph3d() function shorter.
}

\value{
This function returns an intermediate data frame of voxel neighbours that is subsequently pared down as morph3d() continues.
}

\references{
No references at this time.
}

\author{
Tarmo K. Remmel
}

\note{
No additional notes at this time.
}

\seealso{
See Also \code{\link{morph3d}}
}

\examples{
# THIS IS AN INTERNAL FUNCITON AND USAGE EXAMPLES ARE NOT NECESSARY
}

\keyword{ manip }
