\name{morph3dprep}

\alias{morph3dprep}

\title{
Prepare 3D voxel data for interactive 3D plotting.
}

\description{
Prepares 3D array data for plotting.
}

\usage{
morph3dprep(INCUBE = NULL, ORIG = FALSE)
}

\arguments{
  \item{INCUBE}{
This is the input 3D array of voxels.
}
  \item{ORIG}{
A Boolean flag to control whether the data is the original input binary data (TRUE), or some other morphological segmentation (FALSE).
}
}

\details{
No additional details at this time.
}

\value{
The output is...
}

\references{
Remmel, T.K. 2021 (Unpublished manuscript).
}

\author{
Tarmo K. Remmel
}

\note{
No additional notes at this time.
}

\seealso{
See Also \code{\link{morph3d}}, \code{\link{morph3dplot}}
}

\examples{
out <- morph3dprep(INCUBE = LEdemo, ORIG = FALSE)
print(out)
}

\keyword{ manip }
