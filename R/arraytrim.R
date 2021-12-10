arraytrim <-
function(VOLOBJ=NULL) {

  #--------------------------------------------------------------
  #
  # TITLE:     arraytrim()
  # FILENAME:  morph3d_Final.r
  # AUTHOR:    TARMO K REMMEL
  # DATE:      26 November 2021
  # CALLS:     NA
  # CALLED BY: NA
  # NEEDS:     NA
  # NOTES:     Used to trim 0 (zeros, white-space from a 3D array).
  # ARGS:      VOLOBJ: a binary 3D array (0,1). If there are complete
  #            planes of 0 on the outter margins of the array, these
  #            are trimmed off.
  #
  #--------------------------------------------------------------

  # INPUT 3D OBJECT (0,1), WHERE 1 IS THE FEATURE OF INTEREST AND 0 OTHERWISE

  # STORE THE DIMENSIONS OF THE INPUT 3D ARRAY
  dimensions <- dim(VOLOBJ)
  dimx <- dimensions[1]
  dimy <- dimensions[2]
  dimz <- dimensions[3]

  # PREPARE VECTORS THAT WILL INDICATE WHICH PLANES TO TRIM
  # INITIALIZE ALL TO -1
  # LATER -1 MEANS DO NOT CUT, >0 IDENTIFIES THE PLANE THAT COULD BE CUT, IF IT IS OUTSIDE AND CONNECTED TO THE EDGE
  trimx <- rep(x=-1, times=dimx)
  trimy <- rep(x=-1, times=dimy)
  trimz <- rep(x=-1, times=dimz)

  # INDENTIFY WHICH PLANES NEED TO BE CUT AND UPDATE THE APPROPRIATE VECTORS FOR X,Y,Z DIMENSIONS
  for(plane in 1:dimx) {
    if(sum(VOLOBJ[plane,,]) == 0) {
        trimx[plane] <- plane
    }
    if(sum(VOLOBJ[,plane,]) == 0) {
      trimy[plane] <- plane
    }
    if(sum(VOLOBJ[,,plane]) == 0) {
        trimz[plane] <- plane
    }
  }

  # IDENTIFY LOW-END INDICES FOR PERFORMING THE TRIMMING BASED ON THE THREE VECTORS PRODUCED ABOVE
  xindexlow <- NULL
  for(step in 1:length(trimx)) {
    if(trimx[step] != -1) {
      xindexlow <- step
    } else {
      break
    }
  }
  xindexlow <- xindexlow + 1

  yindexlow <- NULL
  for(step in 1:length(trimy)) {
    if(trimy[step] != -1) {
      yindexlow <- step
    } else {
      break
    }
  }
  yindexlow <- yindexlow + 1

  zindexlow <- NULL
  for(step in 1:length(trimz)) {
    if(trimz[step] != -1) {
      zindexlow <- step
    } else {
      break
    }
  }
  zindexlow <- zindexlow + 1

  # IDENTIFY HIGH-END INDICES FOR PERFORMING THE TRIMMING BASED ON THE THREE VECTORS PRODUCED ABOVE
  xindexhigh <- NULL
  for(step in length(trimx):1) {
    if(trimx[step] != -1) {
      xindexhigh <- step
    } else {
      break
    }
  }
  xindexhigh <- xindexhigh - 1

  yindexhigh <- NULL
  for(step in length(trimy):1) {
    if(trimy[step] != -1) {
      yindexhigh <- step
    } else {
      break
    }
  }
  yindexhigh <- yindexhigh - 1

  zindexhigh <- NULL
  for(step in length(trimz):1) {
    if(trimz[step] != -1) {
      zindexhigh <- step
    } else {
      break
    }
  }
  zindexhigh <- zindexhigh - 1

  # CORRECT FOR NULL DIMENSIONS WHEN TRIMMING IS NOT REQUIRED
  cat("\n", xindexlow, xindexhigh, yindexlow, yindexhigh, zindexlow, zindexhigh, sep=" - ")
  if(length(xindexlow)==0) {
    xindexlow <- 1
  }
  if(length(yindexlow)==0) {
    yindexlow <- 1
  }
  if(length(zindexlow)==0) {
    zindexlow <- 1
  }
  
  
  if(length(xindexhigh)==0) {
    xindexhigh <- dimx
  }
  if(length(yindexhigh)==0) {
    yindexhigh <- dimy
  }
  if(length(zindexhigh)==0) {
    zindexhigh <- dimz
  }
  
  cat("\n", xindexlow, xindexhigh, yindexlow, yindexhigh, zindexlow, zindexhigh, sep=" - ")
  
  # PERFORM AND RETURN THE SUBSET
  return(VOLOBJ[xindexlow:xindexhigh, yindexlow:yindexhigh, zindexlow:zindexhigh])

  } # END FUNCTION: arraytrim


