morph3dlinks <-
function(VOLOBJ=NULL, VOXELIDS=NULL, VERBOSE=FALSE) {
  
  #-------------------------------------------------------------------------------
  #
  #  FILENAME:     morph3d_Final.r
  #  FUNCTION:     morph3dlinks()
  #  AUTHOR:       Tarmo K. Remmel
  #  DATE:         6 October 2021
  #  CALLS:        NA
  #  CALLED BY:    NA
  #  NEEDS:        This function requires a binary (0,1) expanded volumetric object
  #                where 0 is the background and 1 is the object of interest.
  #                The object represents 3D space as an array object for R. Convention
  #                has the X for Rows, Y for Columns, and Z for either vertical or
  #                time components.
  #  NOTES:        
  #
  #  REFS:
  #  FUNDING:      NSERC DG to Tarmo K. Remmel
  #
  #-------------------------------------------------------------------------------
 
  if(VERBOSE) {
    cat("\nStarting 3D Morphological Segmentation on object: ", substitute(VOLOBJ), ".\n\n", sep="")
  }
  
  # STORE EXPANDED ARRAY DIMENSIONS
  lrgarraydim <- dim(VOLOBJ)
  
  # CURRENTLY FIXED SIZE; NEED TO MAKE THIS AUTOMATIC
  # INITIALIZE SHIFTING ARRAYS TO ZEROS (0)
  up <- VOLOBJ * 0
  down <- VOLOBJ * 0
  left <- VOLOBJ * 0
  right <- VOLOBJ * 0
  forward <- VOLOBJ * 0
  backward <- VOLOBJ * 0

  lrgvoxelIDS <- VOLOBJ * 0
  lrgvoxelIDS[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)] <- VOXELIDS
   
  #SHIFT DOWN
  down[,,2:(lrgarraydim[3]-1)] <- lrgvoxelIDS[,,3:lrgarraydim[3]]
  down <- down[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)]
  
  #SHIFT UP
  up[,,2:(lrgarraydim[3]-1)] <- lrgvoxelIDS[,,1:(lrgarraydim[3]-2)]
  up <- up[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)]

  #SHIFT LEFT
  left[,2:(lrgarraydim[2]-1),] <- lrgvoxelIDS[,3:lrgarraydim[2],]
  left <- left[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)]
  #SHIFT RIGHT
  right[,2:(lrgarraydim[2]-1),] <- lrgvoxelIDS[,1:(lrgarraydim[2]-2),]
  right <- right[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)]

  #SHIFT FORWARD
  forward[2:(lrgarraydim[1]-1),,] <- lrgvoxelIDS[3:lrgarraydim[1],,]
  forward <- forward[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)]
  #SHIFT BACKWARD
  backward[2:(lrgarraydim[1]-1),,] <- lrgvoxelIDS[1:(lrgarraydim[1]-2),,]
  backward <- backward[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)]
  
  # EXTRACT JUST THE ORIGINAL SIZE FROM THE LARGE ARRAY
  VOLOBJ <- VOLOBJ[2:(lrgarraydim[1]-1), 2:(lrgarraydim[2]-1), 2:(lrgarraydim[3]-1)]
  
  # AGGREGATE RESULTS INTO A DATA FRAME
  new <- as.data.frame(cbind(VOLOBJ, VOXELIDS, up, down, left, right, forward, backward))
  new <- new[new$VOLOBJ==1,]

  # RETURN RESULT AS A DATA FRAME
  return(new)
  
}
