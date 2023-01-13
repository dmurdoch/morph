morph3dplot <-
function(data=NULL, CELLID=TRUE, LEGEND=FALSE, ORIGTRANSP=TRUE) {

  #--------------------------------------------------------------
  #
  # TITLE:     morph3dplot()
  # FILENAME:  morph3d_Final.r
  # AUTHOR:    TARMO K REMMEL
  # DATE:      11 November 2021
  # CALLS:     clear3d, bg3d, cube3d, shade3d, wire3d, legend3d
  # CALLED BY: NA
  # NEEDS:     rgl
  # NOTES:     Given a 3D morphological segmentation (3D array that
  #            that has been melted into an x,y,z,value data.frame,
  #            produces a 3D interactive plot for visualization.
  # ARGS:      LEGEND = TRUE|FALSE (Should a legend be drawn?)
  #            data = a data.frame with x,y,z,value columns that are
  #            integers. The values need to be as follows:
  #              OUTSIDE - TRANSPARENT - 1
  #              MASS - GREEN - 2
  #              SKIN - BLACK - 3
  #              CRUMB - BROWN - 4
  #              CONNECTOR - ORANGE - 5
  #
  #          morp    Good for plotting original with cell IDs and slight transparency
  #              morph3dplot(orig, CELLID=TRUE, ORIGTRANSP=TRUE)
  #--------------------------------------------------------------

  # SAVE GRAPHIC PARAMETERS AND RESTATE THEM ON EXIT
  opar <- par3d(no.readonly=TRUE)
  on.exit(par3d(opar))

  par3d(skipRedraw = TRUE)
  
  # GENERATE CUSTOM COLOUR PALETTE
  cols <- c("grey", "green", "black", "brown", "orange", "pink", "cornflowerblue", "navy", "seagreen", "olivedrab")

  # CLEAR THE DISPLAY AND SET THE BACKGROUND TO WHITE
  clear3d()
  bg3d('white')

  # AN AN ALPHA COLUMN TO CONTROL TRANSPARENCY PER EACH VOXEL (TIED TO COLOUR)
  data <- cbind(data,0.1)
  names(data)[5] <- "alph"
  # OUTSIDE (TRANSPARENT = 1)
  data$alph[data$value == 1] <- 0.1
  # MASS (GREEN = 2)
  data$alph[data$value == 2] <- 0.7
  # SKIN (BLACK = 3)
  data$alph[data$value == 3] <- 0.1
  # CRUMB (BROWN = 4)
  data$alph[data$value == 4] <- 0.8
  # CIRCUIT (ORANGE = 5)
  data$alph[data$value == 5] <- 0.3
  # ANTENNA (PINK = 6)
  data$alph[data$value == 6] <- 0.5
  # BOND (CORNFLOWERBLUE = 7)
  data$alph[data$value == 7] <- 0.3
  # VOID-VOLUME (NAVY = 8)
  data$alph[data$value == 8] <- 0.9
  # VOID (SEAGREEN = 9)
  data$alph[data$value == 9] <- 0.5
  
  if(ORIGTRANSP) {
    data$alph[data$value == 2] <- 0.5
  }

  # REPEAT FOR EACH VOXEL
  for(j in 1:nrow(data)) {

    # SET THE COLOUR FOR THE CURRENT VOXEL
    col <- cols[data$value[j]]

    # DRAW THE CURRENT VOXEL
    cubit <- cube3d(color=col, alpha=0.1)

    # cubit$vb is a 4 by n matrix of vertices in homogeneous coordinates. Each column is a point.
    # Cube is originally 2x2x2. We shrink it to a 1x1x1 cube in the positive octant.
    cubit$vb[cubit$vb == -1] <- 0
    # ADD X TO FIRST ROW
    cubit$vb[1,] <- cubit$vb[1,] + data$x[j]
    # ADD Y TO SECOND ROW
    cubit$vb[2,] <- cubit$vb[2,] + data$y[j]
    # ADD Z TO THIRD ROW
    cubit$vb[3,] <- cubit$vb[3,] + data$z[j]

    # DO NOT SHADE BACKGROUND VOXELS
    if(data$value[j] > 1) {
      shade3d(cubit, add=TRUE, alpha=data$alph[j])
    }

    # DRAW WIREFRAME OUTLINES OF VOXELS
    wire3d(cubit, add=TRUE, color=col)

    # ADD LEGEND IF DESIRED
    if(LEGEND) {
      legend3d("topright", legend = c("Outside", "Mass", "Skin", "Crumb", "Connection"), pch = 15, col = cols, inset=c(0.02))
    }

  } # END FOR: j

  # IF CELL IDS ARE REQUESTED
  if(CELLID) {
    if(ORIGTRANSP) {
      txtcolor <- "black"
    }
    else {
      txtcolor <- "purple"
    }
    coords <- data[,1:3]
    alph <- data[,4]
    alph[alph==1] <- 0.01
    alph[alph>0.01] <- 0.9
    text3d(coords+0.5, texts=seq(1:prod(dim(data))), cex=0.75, col=txtcolor, alpha=alph)
  } # END IF

} # END FUNCTION: morph3dplot
