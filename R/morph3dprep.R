morph3dprep <-
function(INCUBE=NULL, ORIG=FALSE) {

  #--------------------------------------------------------------
  #
  # TITLE:     morph3dprep()
  # FILENAME:  morph3d_Final.r
  # AUTHOR:    TARMO K REMMEL
  # DATE:      30 August 2021
  # CALLS:     melt()
  # CALLED BY: NA
  # NEEDS:     reshape2
  # NOTES:     Given a 3D array, the dimensions are used to melt(),
  #            that is restructure, the data into a data frame with
  #            x,y,z,value columns of integers. This is required to
  #            convert the output from morph3d2() for use with
  #            morph3dplot().
  #
  #            For the original 3D array (0,1), to prepare for plotting
  #            do the following (assume 'stk' is the array):
  #
  #            orig <- morph3dprep(stk, ORIG=TRUE)
  #            morph3dplot(orig)
  #
  # ARGS:      INCUBE = A 3D array with categories representing the
  #            3D morphological elements that will be plotted.
  #
  #            ORIG = TRUE|FALSE (True when preparing to plot the original array,
  #            and FALSE when preparing to plot the 3D morphological elements.
  #
  #--------------------------------------------------------------
  
  # MELT THE CUBE INTIO x,y,z,value FORMAT
  cubexyz <- as.integer(unlist(melt(INCUBE)))

  # SET THE CORRECT DIMENSIONS
  dim(cubexyz) <- c(prod(dim(INCUBE)), 4)

  # CONVERT TO DATA FRAME
  cubexyz <- as.data.frame(cubexyz)

  # ADD APPROPRIATE HEADINGS
  names(cubexyz) <- c("x", "y", "z", "value")

  # IF THE ORIGINAL (0,1) 3D ARRAY, +1 TO VALUES TO AVOID THE ZERO (FOR PLOTTING)
  if(ORIG) {
    cubexyz$value <- cubexyz$value + 1
  }

  # RETURN THE RESTRUCTURED 3D ARRAY AS x,y,z,value DATA.FRAME
  return(cubexyz)

} # END FUNCTION: morph3dprep
