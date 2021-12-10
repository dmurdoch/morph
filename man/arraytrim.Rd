\name{arraytrim}

\alias{arraytrim}

\title{
Trim excess white-space (0 planes) from a 3D binary array.
}

\description{
A helper function for reducing the size of 3D arrays if necessary. This reduction can improve memory utilization and speed processing, since voxels outside the feature of interest are eliminated.
}

\usage{
arraytrim(VOLOBJ = NULL)
}

\arguments{
  \item{VOLOBJ}{
This is a 3D binary (0,1) numeric array, which represent voxels. If there are outter margins that are all zeros, they will be trimmed off to reduce the size of the array.
}
}

\details{
The function scans inward, starting at the outer margins, for planes that are all zeros. If they are, those planes will be trimmed at the end. Only planes of zeros that connect contiguously to the outer margins will be trimmed.
}

\value{
The output is an array with the excess outter zero planes removed.
}

\references{
Remmel, T.K. 2021. (Unpuplished manuscript)
}

\author{
Tarmo K. Remmel
}

\note{
It can happen that trimming reduces the resulting dimensions to be <3, and thus, prior to running morph3d(), ensure that your data still fits the requirements.
}

\seealso{
See Also \code{\link{morph3d}}, \code{\link{morph3dprep}} 
}

\examples{
LEdemoSmall <- arraytrim(LEdemo)
}

\keyword{ array }
