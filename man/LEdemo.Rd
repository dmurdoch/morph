\name{LEdemo}
\alias{LEdemo}
\docType{data}
\title{
Landscape Ecology Demo Data (3D Voxels)
}
\description{
This is a 3D numeric array with dimensions [10,10,10], representing 1000 voxels for demonstration purposes. The data was manually produced to ensure that all morphological classes are present, to allow demonstration and illustration within the Landscape Ecology paper in which the morph3d() functionality was introduced.
}
\usage{data("LEdemo")}
\format{
  The format is:
 num [1:10, 1:10, 1:10] 1 1 0 0 0 0 0 0 0 0 ...
}
\details{
This data can be used for testing, demonstrations, or for quickly getting the package up and running, before getting your own data into the proper format.
}
\source{
This data was produced by the author of the package for demonstration purposes and was used in the paper referenced below as the example data.
}

\references{
Remmel, T.K. 2021. (Unpublished manuscript)
}

\examples{
data(LEdemo)
str(LEdemo)
}

\keyword{datasets}
