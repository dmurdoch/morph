morph3d <-
function(DATACUBE=NULL, VERBOSE=FALSE, PLOT=FALSE, FINALPLOT=TRUE, PLOTIDS=FALSE) {

  #--------------------------------------------------------------
  #
  # TITLE:     morph3d()
  # FILENAME:  morph3d_Final.R
  # AUTHOR:    TARMO K REMMEL
  # DATE:      15 November 2021 @ 1132h
  # CALLS:
  # CALLED BY:
  # NEEDS:     igraph, stringr, rgl
  # NOTES:     NOTE, INPUT DATA IS AS FOLLOWS:
  #  DATACUBE MUST BE A 3D ARRAY OF (0,1) WHERE 1 IS THE OBJECT OF INTEREST AND 0 IS NOT
  # IS.ARRAY = TRUE, IS.NUMERIC = TRUE
  #
  # ARGS:      DATACUBE = A 3D array with categories representing a
  #            a binary feature that is to be processed. This needs to
  #            be a numeric array (it does not need to be integer).
  #
  #            VERBOSE = TRUE|FALSE (True when more output should be
  #            provided on screen as the script runs, and FALSE when
  #            to run in quiet mode).
  #
  #--------------------------------------------------------------
    
  # CHECK THAT INPUT DATA IS 0,1 (BINARY) 3D NUMERIC ARRAY DATA
  
  # ** UPDATE: ADD CHECK THAT THERE ARE AT LEAST 2 VOXELS OF DATA (I.E., A SINGLE VOXEL CRUMB CAUSES A PROBLEM)
  
  if(is.numeric(DATACUBE)) {
    if(is.array(DATACUBE)) {
      if(length(dim(DATACUBE)) == 3) {
        if(all(as.integer(names(table(DATACUBE))) %in% c(0,1))) {
          message("\n\nInput data passess all initial checks for integrity.\n\n")
        } else {
          warning("\n\nERROR 004 - Input data does not contain proper 0,1 binary data.\n\n")
          return(NULL)
        }
      } else {
        warning("\n\nERROR 003 - Input data is not a 3D array.\n\n")
        return(NULL)
      }
    } else {
      warning("\n\nERROR 002 - Input data is not an array.\n\n")
      return(NULL)
    }
  } else {
    warning("\n\nERROR 001 - Input data is not numeric.\n\n")
    return(NULL)
  }
   
   
  # DEFINE A HELPER FUNCTION
  `%notin%` <- Negate(`%in%`)
  
  if(VERBOSE) {
  	cat("\nInsetting 3D data cube into a larger array to handle edge effects")
  } # END IF
  # STORE THE DIMENSIONS OF THE INPUT DATA CUBE
  dimdatacube <- dim(DATACUBE)

  # MAKE A LARGER DATA CUBE AND INSET THE ORIGINAL DATA CUBE INTO IT SUCH THAT THERE
  # ARE 0 VALUES ON ALL CUBE MARGINS TO CONTAIN THE DATA
  lrgdatacube <- array(data=0, dim=c(dimdatacube[1]+2, dimdatacube[2]+2, dimdatacube[3]+2))
  # MAKE A SECOND LARGEER DATA CUBE FOR SEPARATING SKIN FROM VOID LATER AND FILL IT WITH -1 RATHER THAN 0
  lrgdatacube2 <- lrgdatacube - 1
  # PERFORM INSETS
  lrgdatacube[2:(dimdatacube[1]+1),2:(dimdatacube[2]+1),2:(dimdatacube[3]+1)] <- DATACUBE
  #lrgdatacube
  lrgdatacube2[2:(dimdatacube[1]+1),2:(dimdatacube[2]+1),2:(dimdatacube[3]+1)] <- DATACUBE
  #print(lrgdatacube2)
  
  # PERFORM INITIALIZATIONS
  if(VERBOSE) {
  	cat("\n\nPerforming initializations of: voxelID, objectID, coreCode, and morphCode arrays")
  } # END IF
  # BUILD AN ARRAY FOR HOLDING THE INDEX OF VOXELS IDS
  voxelID <- array(data=1:prod(dim(DATACUBE)), dim=c(dimdatacube[1],dimdatacube[2],dimdatacube[3]))
  
  # INITIATE AN ARRAY FOR HOLDING UNIQUE OBJECT IDS
  objectID <- voxelID * 0

  # INITIATE AN ARRAY FOR HOLDING MASS-CORE CODES
  coreCode <- objectID * 0
  
  # INITIATE AN ARRAY FOR HOLDING EXPANDED CORE CODES
  expandedCoreCode <- coreCode
  
  # INITIATE AN ARRAY FOR HOLDING MORPHOLOGY CODES
  morphCode <- objectID

  # COMPUTE NEIGHBOURS
  if(VERBOSE) {
    cat("\nComputing neighbours")
  }
  neighbourtab <- morph3dlinks(VOLOBJ=lrgdatacube, VOXELIDS=voxelID)
  neighbourtab <- neighbourtab[,2:8]
    
  # MAINTAIN ONLY VOXEL IDS THAT EXIST IN THE FEATURE (1)
  goodIDs <- voxelID[DATACUBE==1]
  for(row in 1:dim(neighbourtab)[1]) {
    # IF THE VOXEL ID IS NOT PART OF THE FEATURE, SET IT TO ZERO
    neighbourtab[row,][!neighbourtab[row,] %in% as.vector(voxelID[DATACUBE==1])] <- 0
  } # END FOR: row
    
  # NOW RESTRUCTURE neighbourstab TO CREATE ALL OF THE PAIRS, SO FIRST 2 COLUMNS, THESE FOLLOW WITH FIRST COLUMN AND SECOND...
  newneighbourtab <- t(t(rep(neighbourtab[,1], 6)))
  newneighbourtab <- cbind(newneighbourtab, as.vector(unlist(neighbourtab[,2:7])))
    
  # REMOVE DUPLICATES DUE TO DIRECTION
  newneighbourtab <- newneighbourtab[!duplicated(t(apply(newneighbourtab, 1, sort))),]
    
  # PRODUCE MAIN GRAPH OF ALL VOXELS
  if(VERBOSE) {
    cat("\n  Generating network graph object")
  }
  maingraph <- graph_from_data_frame(newneighbourtab, directed=FALSE)

  # REMOVE THE EDGE SPLIT
  if("0" %in% names(V(maingraph))) {
    cutgraph <- delete_vertices(maingraph, "0")
  }
  else {
    # FOR NOW, IF THERE IS NO "0" VERTEX TO REMOVE, JUST COPY THE maingraph INTO THE INTO cutgraph AND CONTINUE ON
    cutgraph <- maingraph
  } # END IF-ELSE
    
  # DECOMPOSE THE cutgraph INTO A LIST OF AS MANY GRAPHS AS THERE EXIST IN THE STRUCTURE
  decompgraph <- decompose(cutgraph)
  if(VERBOSE) {
    cat("\n  There are", length(decompgraph), "discrete objects in this graph", sep=" ")
  }
  
  # PERFORM THE MORPHOLOGICAL SEGMENTATION
  if(VERBOSE) {
    cat("\n  Initiating the 3D morphological segmentation")
  } # END IF
  
  # ASSIGN UNIQUE IDS TO EACH DISCRETE CLUSTER/OBJECT: RESULTS ARE STORED IN objectID
  for(uq in 1:length(decompgraph)) {
    # VECTOR OF INTEGERS OF VOXELS BELONGING TO UNIQUE CLUSTER
    b <- as.integer(names(degree(decompgraph[[uq]])))
    numvoxels <- length(b)
    for(vox in 1:numvoxels) {
      # ADD UNIQUE VALUE uq TO EACH IDENTIFIED CLUSTER
      objectID[voxelID == b[vox]] <- uq
    } # END FOR: vox
  } # END FOR: uq

  if(PLOT) {
    # THIS IS A PLOT OF THE DISCRETE OBJECTS
    origclust <- morph3dprep(objectID, ORIG=TRUE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(origclust, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
    bgplot3d({
      plot.new()
      title(main = 'Object IDs', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT DISCRETE OBJECT CODES

  
  # ------------ MASS (CODE = 2) ------------
  
  
  # MASS VOXELS WILL HAVE 6 CONNECTIONS
  if(VERBOSE) {
    cat("\n\n  Identifying MASS voxels")
  } # END IF
  for(uq in 1:length(decompgraph)) {
    if(VERBOSE) {
      cat("\n    Processing unique object", uq, "of", length(decompgraph), sep=" ")
    }
    # VECTOR OF INTEGERS OF VOXELS BELONGING TO UNIQUE CLUSTER
    sel <- degree(decompgraph[[uq]])[degree(decompgraph[[uq]])==6]
    
    # DO THIS ONLY IF THERE ARE MASS VOXELS TO PROCESS
    if(length(sel) > 0){
      vox <- as.integer(names(sel))
      if(VERBOSE) {
        cat("\n      There are", length(sel), "MASS voxels in this unique object", sep=" ")
      }
      for(rep in 1:length(sel)) {
        # ADD CODE = 2 TO ALL MASS VOXELS
        morphCode[voxelID == vox[rep]] <- 2
      } # END FOR: rep
    } else {
      if(VERBOSE) {
        cat("\n      There are no MASS voxels in this unique object")
      }
    } # END IF
  } # END FOR: uq

  if(PLOT) {
    # THIS IS A PLOT OF THE MASS VOXELS
    plotmorph <- morphCode
    plotmorph[plotmorph!=2] <- 0
    plotmorph[plotmorph==2] <- 1
    morphs <- morph3dprep(plotmorph, ORIG=TRUE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=FALSE)
    bgplot3d({
      plot.new()
      title(main = 'MASS Voxels', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT MASS VOXELS

  
  # ------------ SKIN (CODE = 3) ------------
  
  
  # VOXELS MUST SHARE AN EDGE WITH A MASS VOXEL
  # ENSURE THAT A VOXEL IS NOT MASS ALREADY (THAT IS, DO NOT OVERWRITE IT WITH A SKIN CODE)
  if(VERBOSE) {
    cat("\n\n  Identifying SKIN voxels")
  }
  
  for (uq in 1:length(decompgraph)) {
    if(VERBOSE) {
      cat("\n    Processing unique object", uq, "of", length(decompgraph), sep=" ")
    }
  
    # VECTOR OF INTEGERS OF VOXELS BELONGING TO UNIQUE CLUSTER
    voxelIDS <- as.integer(names(degree(decompgraph[[uq]])))
    numvoxels <- length(voxelIDS)
    voxelMASS <- degree(decompgraph[[uq]])[degree(decompgraph[[uq]])==6]
    massIDS <- names(voxelMASS)
    numMASS <- length(massIDS)

    if(numMASS > 0) {
      for(rep in 1:numMASS) {
        # FIND ALL NEIGHBOURING VOXELS TO MASS VOXEL
        MASSneigh <- names(neighbors(decompgraph[[uq]], massIDS[rep]))
        # RETAIN ONLY NEIGHBOURS THAT ARE NOT ACTUALLY MASS VOXELS
        MASSneigh <- MASSneigh[!MASSneigh %in% massIDS]
        # WRITE SKIN CODE = 3) TO morphCode OBJECT FOR REMAINING VOXEL IDS ONLY WHEN THE LENGTH OF THE MASSneigh VECTOR IS > 0, (I.E., THERE ARE VOXELS TO PROCESS)
        if(length(MASSneigh) > 0) {
          for(newrep in 1:length(MASSneigh)) {
            morphCode[voxelID == MASSneigh[newrep]] <- 3
          } # END FOR: newrep
        } # END IF
      } # END FOR: rep
    } # END IF
  } # END FOR: uq
  
   
  # ------------ CRUMB (CODE = 4) ------------
  
  
  # VOXEL GROUPS ARE NOT CONNECTED
  # VOXELS GROUPS MUST NOT CONTAIN MASS
  # VOXEL GROUPS MUST NOT TOUCH SKIN ANYWHERE
  # MAKE VECTORS OF VERTEX IDS THAT ARE MASS AND ONE THAT ARE SKIN
  # IN EACH uq LOOP, ENSURE THAT NONE OF THE MASS ARE IN THE LIST OF VERTICES AND THAT NONE OF THE EDGES TOUGH THE SKIN IDS,
  # IF THESE TWO CHECKS PASS, THEN WE HAVE A CRUMB
  
  # VECTORS OF ALL MASS AND ALL SKIN VOXEL IDS
  allMASS <- voxelID[morphCode==2]
  allSKIN <- voxelID[morphCode==3]
  
  # CYCLE THROUGH EACH DISCRETE OBJECT
  if(VERBOSE) {
    cat("\n\n  Identifying CRUMB voxels")
  }
  for (uq in 1:length(decompgraph)) {
    if(VERBOSE) {
      cat("\n    Processing unique object", uq, "of", length(decompgraph), sep=" ")
    }
    if(length(V(decompgraph[[uq]]))==1) {
      if(VERBOSE) {
        cat("\n      Object has only 1 voxel, thus it is a CRUMB", sep="" )
      }
      morphCode[voxelID == as.integer(names(V(decompgraph[[uq]])))] <- 4
    } else {
      # PRODUCE VECTOR OF CURRENT VOXEL IDS FOR THE OBJECT IN QUESTION
      curVoxels <- as.integer(names(V(decompgraph[[uq]])))
      if(VERBOSE) {
        cat("\n      Object has more than 1 voxels. It has: ", length(curVoxels), sep="" )
      }
      if(TRUE %in% (allMASS %in% curVoxels)) {
        # POSSIBLE CRUMB: TEST FAILS; CHECKED IF CONNECTED TO ANY MASS
        if(VERBOSE) {
          cat("\n      Not a CRUMB since it connects to a MASS")
        }
        } else {
        if(TRUE %in% (allSKIN %in% as.integer(unlist(strsplit(attr(E(decompgraph[[uq]]), "vnames"),"|", fixed=TRUE))))) {
          # POSSIBLE CRUMB TEST FAILS; CHECKED CONNECTION TO ANY SKIN
          if(VERBOSE) {
            cat("\n      Not a CRUMB since it connects to a SKIN")
          }
        } else {
          # ASSIGN ALL VERTICES TO CRUMB
          if(VERBOSE) {
            cat("\n      There are", length(curVoxels), "CRUMB voxels in this object")
          }
          for(rep in 1:length(curVoxels)) {
            morphCode[voxelID == curVoxels[rep]] <- 4
          } # END FOR: rep
        } # END IF-ELSE
      } # END IF-ELSE
    } # END IF
  } # END FOR: uq
  
  if(PLOT) {
    # THIS IS A PLOT OF THE CRUMB VOXELS
    plotmorph <- morphCode
    plotmorph[plotmorph!=4] <- 0
    plotmorph[plotmorph==4] <- 3
    morphs <- morph3dprep(plotmorph, ORIG=TRUE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=FALSE)
    bgplot3d({
      plot.new()
      title(main = 'CRUMB Voxels', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT THE MORPHOLOGY WITH CRUMBS

  
  # ------------ CONNECTOR (CODE = 5) ------------
  
  
  # INITIATE THIS PHASE BY ASSIGNING ALL CONNECTOR VOXELS CODE = 5
  # THERE ARE 3 TYPES (ANTENNA, CIRCUIT, BOND) THAT WILL SUBSEQUENTLY BE PARSED
  # ANTENNA - BRANCH (EXTENDS FROM SKIN OF ONE MASS BUT TOUCHES NO OTHER VOXELS) ONLY ONE SKIN IS TOUCHED
  # BOND - BRIDGE (TOUCHES SKIN BETWEEN TWO DISTINCT MASS)
  # CIRCUIT - LOOP (TOUCHES SKIN AT LEAST TWICE FOR A SINGLE MASS)

  if(VERBOSE) {
    cat("\n\n  Identifying CONNECTOR voxels")
  }
  
  # ASSIGN CONNECTOR CODE 5 TO ALL POSSIBLE CONNECTOR TYPES
  # LATER, VOXELS THAT RETAIN CODE 5 WILL BECOME CALLED CIRCUIT
  morphCode[DATACUBE > 0 & morphCode == 0] <- 5
  if(VERBOSE) {
    cat("\n    There are", length(voxelID[morphCode==5]), "generic CONNECTOR voxels", sep=" ")
    cat("\n\n  Identifying OUTSIDE voxels")
  }
  # ADD CODE = 1 TO ALL REMAINING 0 VOXELS TO MAKE THE OUTSIDE
  morphCode[morphCode == 0] <- 1
  if(VERBOSE) {
    cat("\n    There are", length(voxelID[morphCode==1]), "OUTSIDE voxels", sep=" ")
  }
  
  
  # ------------ ANTENNA (CODE = 6) ------------
  
  if(VERBOSE) {
    cat("\n\n  Identifying ANTENNA voxels")
  }
  
  # LIST OF MASS VOXELS, SKIN, AND CONNECTOR VOXELS
  # MUST ORIGINATE FROM A SKIN CONNECTION AND CONNECT AT ONLY ONE POINT
  # ANTENNA-LIKE FEATURES THAT DO NOT ORIGINATE AT SKIN ARE CIRCUITS
  MASSVoxels <- voxelID[morphCode==2]
  SKINVoxels <- voxelID[morphCode==3]
  CONNECTORVoxels <- voxelID[morphCode==5]

  # CYCLE THROUGH EACH OBJECT (UNIQUE 3D OBJECTS)
  for(obj in 1:length(decompgraph)) {
    if(VERBOSE) {
    cat("\n    Processing object ", obj, " of ", length(decompgraph), " in gph.", sep="")
    }
    # TEST WHETHER MASS VOXELS EXIST IN THE SUBGRAPH
    MassVoxelsInObject <- MASSVoxels[MASSVoxels %in% names(V(decompgraph[[obj]]))]
    SkinVoxelsInObject <- SKINVoxels[SKINVoxels %in% names(V(decompgraph[[obj]]))]
    ConnectorVoxelsInObject <- CONNECTORVoxels[CONNECTORVoxels %in% names(V(decompgraph[[obj]]))]
    
    if(VERBOSE) {
      cat("\n      MassVoxelsInObject:", MassVoxelsInObject, sep=" ")
      cat("\n      SkinVoxelsInObject:", SkinVoxelsInObject, sep=" ")
      cat("\n      ConnectorVoxelsInObject:", ConnectorVoxelsInObject, sep=" ")
    }
    
    # IF THERE IS NO MASS, THERE CANNOT BE A CONNECTOR, SO CHECK FIRST
    if(length(MassVoxelsInObject) > 0 ) {
      # MASS EXISTS IN THIS OBJECT, SO DELETE MASS AND SKIN VOXELS
      if(VERBOSE) {
        cat("\n      Object ", obj, ": There are MASS and SKIN voxels to delete", sep=" ")
      }
      gph2 <- delete_vertices(decompgraph[[obj]], as.character(MassVoxelsInObject))
      gph3 <- delete_vertices(gph2, as.character(SkinVoxelsInObject))
      if(VERBOSE) {
        cat("\n      Decomposing...")
      }
      gph4 <- decompose(gph3)
      # STORE THE NUMBER OF OBJECT SUBPARTS TO PROCESS
      numSubparts <- length(gph4)
      if(VERBOSE) {
        cat("\n      The decomposed graph now has", numSubparts, "sub parts", sep=" ")
      }
      
      # LOOP THROUGH EACH UNIQUE SUB PART FOR A GIVEN OBJECT
      for(connector in 1:numSubparts) {
        if(VERBOSE) {
          cat("\n\n        On sub part", connector, "of", numSubparts, sep=" ")
        }
        # VOXEL IDS AS PART OF THIS CONNECTOR SUBGRAPH
        MassVoxelsInSubGraph <- MASSVoxels[MASSVoxels %in% names(V(gph4[[connector]]))]
        ConnectorVoxelsInSubGraph <- CONNECTORVoxels[CONNECTORVoxels %in% names(V(gph4[[connector]]))]
        SkinVoxelsInSubGraph <- SKINVoxels[SKINVoxels %in% names(V(gph4[[connector]]))]
        if(VERBOSE) {
          cat("\n          Sub-part", connector, ": MassVoxelsInSubGraph:", MassVoxelsInSubGraph, sep=" ")
          cat("\n          Sub-part", connector, ": SkinVoxelsInSubGraph:", SkinVoxelsInSubGraph, sep=" ")
          cat("\n          Sub-part", connector, ": ConnectorVoxelsInSubGraph:", ConnectorVoxelsInSubGraph, sep=" ")
          cat("\n          Vertices:", names(V(gph4[[connector]])), sep=" ")
        }
        
        if(length(names(V(gph4[[connector]]))) == 1) {
          if(VERBOSE) {
            cat("\n          Single voxel antenna\n")
          }
          voxID <- as.integer(names(V(gph4[[connector]])))
          locco <- which(voxelID==voxID, arr.ind=TRUE)
          if(VERBOSE) {
            cat("\n          Single voxel ANTENNA at ID:", voxID, "and locatoin:", locco, sep=" ")
          }
          morphCode[locco] <- 6
          
          
          # SINGLE ANTENNA VOXEL CORRECTOR CHECK (WITHIN SKIN NEIGHBOURS)
          # SHOULD CONFIRM THAT IT DOES NOT HAVE MULTIPLE SKIN NEIGHBOURS. IF MORE THAN 1 MASS NEIGHBOUR, THAN SIMPLY A CIRCUIT
          if(VERBOSE) {
            cat("\n          Single voxel antenna - checking how many MASS contacts", voxID, "has...", sep=" ")
          }
          loc <- which(voxelID==as.integer( names(V(gph4[[connector]])) ), arr.ind=TRUE)
          if(VERBOSE) {
            cat("\n          Single voxel antenna - checking", loc, sep=" ")
          }
          # CHECK ALL SIX NEIGHBOURS AND COUNT (DO NOT EXTEND BEYOND THE ARRAY LIMITS)
          counter <- 0
          xlow <- 0
          xhigh <- 0
          ylow <- 0
          yhigh <- 0
          zlow <- 0
          zhigh <- 0
          
          # ENSURE THAT I DO NOT LOOK AT A REFERENCE BEYOND THE ARRAY EXTENTS
          if(loc[1] > 1) {
            xlow <- morphCode[(loc[1]-1), loc[2], loc[3]]
            if(xlow == 2) {
              xlow <- 1
            } else {
              xlow <- 0
            }
          }
          if(loc[1] < dimdatacube[1]) {
            xhigh <- morphCode[(loc[1]+1), loc[2], loc[3]]
            if(xhigh == 2) {
              xhigh <- 1
            } else {
              xhigh <- 0
            }
          }
          if(loc[2] > 1) {
            ylow <- morphCode[loc[1], (loc[2]-1), loc[3]]
            if(ylow == 2) {
              ylow <- 1
            } else {
              ylow <- 0
            }
          }
          if(loc[2] < dimdatacube[2]) {
            yhigh <- morphCode[loc[1], (loc[2]+1), loc[3]]
            if(yhigh == 2) {
              yhigh <- 1
            } else {
              yhigh <- 0
            }
          }
          if(loc[3] > 1) {
            zlow <- morphCode[loc[1], loc[2], (loc[3]-1)]
            if(zlow == 2) {
              zlow <- 1
            } else {
              zlow <- 0
            }
          }
          if(loc[3] < dimdatacube[3]) {
            zhigh <- morphCode[loc[1], loc[2], (loc[3]+1)]
            if(zhigh == 2) {
              zhigh <- 1
            } else {
              zhigh <- 0
            }
          }
          counter <- xlow + xhigh + ylow + yhigh + zlow + zhigh
          if(counter==1) {
            # ONLY ONE SKIN CONTACT, SO THIS IS AN ANTENNA
            morphCode[voxelID==ConnectorVoxelsInSubGraph] <- 6
          } else {
            # MORE THAN ONE SKIN CONTACT, SO THIS IS A CIRCUIT
            morphCode[voxelID==ConnectorVoxelsInSubGraph] <- 5
          }
          
          
          # ADD CHECK FOR SINGLE CONNECTION TO SKIN, IF YES, THEN CONVERT TO ANTENNA
          checkcounter <- 0
          xlow <- 0
          xhigh <- 0
          ylow <- 0
          yhigh <- 0
          zlow <- 0
          zhigh <- 0
          
          # ENSURE THAT I DO NOT LOOK AT A REFERENCE BEYOND THE ARRAY EXTENTS
          if(loc[1] > 1) {
            xlow <- morphCode[(loc[1]-1), loc[2], loc[3]]
            if(xlow == 3) {
              xlow <- 1
            } else {
              xlow <- 0
            }
          }
          if(loc[1] < dimdatacube[1]) {
            xhigh <- morphCode[(loc[1]+1), loc[2], loc[3]]
            if(xhigh == 3) {
              xhigh <- 1
            } else {
              xhigh <- 0
            }
          }
          if(loc[2] > 1) {
            ylow <- morphCode[loc[1], (loc[2]-1), loc[3]]
            if(ylow == 3) {
              ylow <- 1
            } else {
              ylow <- 0
            }
          }
          if(loc[2] < dimdatacube[2]) {
            yhigh <- morphCode[loc[1], (loc[2]+1), loc[3]]
            if(yhigh == 3) {
              yhigh <- 1
            } else {
              yhigh <- 0
            }
          }
          if(loc[3] > 1) {
            zlow <- morphCode[loc[1], loc[2], (loc[3]-1)]
            if(zlow == 3) {
              zlow <- 1
            } else {
              zlow <- 0
            }
          }
          if(loc[3] < dimdatacube[3]) {
            zhigh <- morphCode[loc[1], loc[2], (loc[3]+1)]
            if(zhigh == 3) {
              zhigh <- 1
            } else {
              zhigh <- 0
            }
          }
          checkcounter <- xlow + xhigh + ylow + yhigh + zlow + zhigh
                    
          if(checkcounter==1) {
            # ONLY ONE SKIN CONTACT, SO THIS IS AN ANTENNA
            morphCode[locco] <- 6
          } # END IF
  
        } else {
          # MORE THAN ONE VOXEL - CHECK IF SINGLE CONNECTION TO SKIN. IF YES, THEN ANTENNA
          if(VERBOSE) {
            cat("\n          Possible antenna - need to check single point of connection with SKIN")
          }
          
          # ALL COMBINATIONS OF SKIN AND CONNECTOR VOXEL IDS
          combos <- expand.grid(SkinVoxelsInObject, ConnectorVoxelsInSubGraph)
          combos <- as.data.frame(cbind(combos,0))
                      
          if(dim(combos)[1] > 0) {
            for(row in 1:dim(combos)[1]) {
              combos[row,3] <- are_adjacent(decompgraph[[obj]], as.character(combos[row,1]), as.character(combos[row,2]))
            } # END FOR: row
                            
            if(dim(combos[combos[,3]==1,])[1] == 1) {
              if(VERBOSE) {
                cat("\n          The subgraph: ", connector, " is an ANTENNA.", sep="")
              }
              # ADD CODE 6 TO ANTENNA
              for(j in 1:length(ConnectorVoxelsInSubGraph)) {
                morphCode[voxelID==ConnectorVoxelsInSubGraph[j]] <- 6
              } # END FOR
              
            } else {
              if(VERBOSE) {
                cat("\n          The subgraph: ", connector, " is not an ANTENNA.", sep="")
              }
            } # END IF-ELSE
          } # END IF
        } # END IF-ELSE
       
      } # END FOR: connector
      if(VERBOSE) {
        cat("\n      Done with object", obj, sep=" ")
      }
    } else {
      # MASS VOXELS DO NOT EXIST IN THIS OBJECT, IT CANNOT HAVE CONNECTORS
      if(VERBOSE) {
        cat("\n      There are no MASS and SKIN voxels; therefore there are no possible connectors to process")
        cat("\n      Done with object", obj, sep=" ")
      }
    } # END IF-ELSE
    if(VERBOSE) {
      cat("\n\n")
    }
  } # END FOR: obj


  # ------------ BOND (CODE = 7) ------------
  
  if(VERBOSE) {
    cat("\n\n  Identifying BOND voxels")
    cat("\n    Identifying unique cores of MASS voxels")
  }
  # CYCLE THROUGH EACH DISCRETE OBJECT
  for(obj in 1:length(decompgraph)) {
    MASSVoxels <- voxelID[morphCode==2]
    MassVoxelsInObj <- MASSVoxels[MASSVoxels %in% names(V(decompgraph[[obj]]))]
        
    if(length(MassVoxelsInObj)>0) {
      if(VERBOSE) {
        cat("\n    Need to provide unique IDs to each core of MASS")
      }
      # MAKE GRAPH OF ONLY CORES IN THIS OBJECT
      objcoregraph <- decompgraph[[obj]]
      
      # PREPARE A LIST OF VOXELS TO DELETE (NON-MASS) FROM A NEW NETWORK GRAPH
      notMASS <- as.integer(names(V(objcoregraph)))[!as.integer(names(V(objcoregraph))) %in% MassVoxelsInObj]
      tmpgraph <- delete_vertices(objcoregraph, as.character(notMASS))
      # DECOMPOSE GRAPH
      tmpgraph2 <- decompose(tmpgraph)
      numcores <- length(tmpgraph2)
      if(VERBOSE) {
        cat("\n    There are", numcores, "core MASSES in this object", sep=" ")
      }
      
      # FOR EACH SUB-GRAPH, CREATE UNIQUE CORE IDS
      if(numcores>0) {
        for(core in 1:numcores) {
          MassInCore <- as.integer(names(V(tmpgraph2[[core]])))
            coreID <- max(coreCode) + 1
            for(vox in 1:length(MassInCore)) {
              coreCode[voxelID==MassInCore[vox]] <- coreID
            } # END FOR: vox
                        
            # NOW MAKE AN EXPANDED VERSION THAT INCLUEDS MASS, SKIN, AND FIRST LAYER OF CONNECTED CONNECTORS
            expandedCoreCode[voxelID==MassInCore[vox]] <- coreID
            inclSKIN <- adjacent_vertices(decompgraph[[obj]], as.character(MassInCore))
            MASSSKIN <- sapply(str_split(names(unlist(inclSKIN)), "\\.",  n = 2), `[`, 2)
            inclfirstconnector <- adjacent_vertices(decompgraph[[obj]], as.character(MASSSKIN))
            MASSSKINCON <-sapply(str_split(names(unlist(inclfirstconnector)), "\\.",  n = 2), `[`, 2)
            combined <- c(MASSSKIN, MASSSKINCON)
            combined <- unique(combined)

            # WRITE TO OBJECT
            for(vox in 1:length(combined)) {
              expandedCoreCode[voxelID==combined[vox]] <- coreID
            } # END FOR: vox
            
        } # END FOR: numcores
      } # END IF
    } else {
      if(VERBOSE) {
        cat("\n    No MASS in this object; therefore no cores")
      }
    } # END IF-ELSE
  }
  if(VERBOSE) {
    cat("\n\n    Cores, if present have been identified and coded\n")
  }
  
  if(PLOT) {
    # THIS IS A PLOT OF THE EXPANDED CORE VOXELS
    plotmorph <- expandedCoreCode
    morphs <- morph3dprep(plotmorph, ORIG=TRUE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
    bgplot3d({
      plot.new()
      title(main = 'Expanded CORE Voxels', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT THE MORPHOLOGY WITH ANTENNAE

  
  
  # SECONDARY TEST FOR ANTENNA THAT PORTURDE FROM CIRCUITS
  # IF THERE IS A STRING OF CONNECTORS CONNECTED TO A CONNECTOR THAT TOUCHES SKIN BUT NOT SKIN ITSELF, AND THAT CONNECTS TO ONLY ONE CORE, MAKE IT ANTENNA
  # CYCLE THROUGH EACH OBJECT
 
  for(obj in 1:length(decompgraph)) { # DO NOT PROCESS IF obj HAS NO MASS; NEED TO IMPLEMENT THAT CHECK
    if(VERBOSE) {
      cat("\n\n  ANTENNA FIX object:", obj, sep=" ")
    }
    
     # ONLY PROCESS OBJECTS THAT HAVE MASS WITHIN THEM
     # WITHOUT MASS, THEY CANNTO HAVE ANTENNA OR BOND CONNECTIONS
     # IF FLAG IS TRUE, THERE IS NO MASS, SO SKIP
     # IF FLAG IS FALSE, THERE IS MASS, SO PROCESS
     flag <- all(names(table(morphCode==2 & objectID==obj))==FALSE)
     if(flag) {
       if(VERBOSE) {
         cat("\n  Object has no MASS, skipping")
       }
     } else {
       if(VERBOSE) {
         cat("\n  Object has MASS, checking")
       }
       # CONNECTORS THAT ARE NOT PART OF THE EXPANDED CORE
       tmpgph <- decompgraph[[obj]]
       # DELETE THE EXPANDED CORE VOXELS (GET THE VOXEL IDS WHERE WE HAVE CONNECTORS OUTSIDE OF EXPANDED CORES)
       possible <- voxelID[morphCode==5 & expandedCoreCode == 0 & objectID == obj]
       allinobj <- names(V(decompgraph[[obj]]))
       todelete <- allinobj[allinobj %notin% possible]
           
       # DELETE SAID VOXELS FROM THE GRAPH
       tmpgph <- delete_vertices(decompgraph[[obj]], allinobj[allinobj %notin% possible])
      
       tmpgph2 <- decompose(tmpgph)
       if(VERBOSE) {
         cat("\n\n  There are", length(tmpgph2), "possible ANTENNA to test for connectivity", sep=" ")
       }
       # IF DEGREE IS ZERO, LEAVE AS CIRCUIT
       # IF DEGREE IS 1, IT IS THE END OF A STRING THAT NEEDS TO BE TESTED FOR CONNECTIVITY AS AN ANTENNA (COULD POSSIBLY CHECK TWO CONNECTIONS TO TWO CORES FOR BOND HERE TOO)
       # IF DEGREE IS >2, LEAVE IT FROM THE CHECKS, IT WILL TAKE WHATEVER CODE IS GIVEN TO THAT SUBGRAPH

       # CYCLE THROUGH EACH SUB CHAIN
       for(chain in 1:length(tmpgph2)) {
         if( length(V(tmpgph2[[chain]]) > 1)  ) {
           # IF THE CHAIN BEING ASSESSED AS AN ANTENNA HAS MORE THAN ONE VOXEL, ASSESS THE DEGREE OF EACH VOXEL
           # TO IDENTIFY THE ENDPOINTS NEEDED TO BE TESTED FOR CONNECTIVITY TO EXPANDED MASS CORES
           endpoints <- names(degree(tmpgph2[[chain]]))[degree(tmpgph2[[chain]])==1]
           
           # CHECK HERE IF ENDPOINTS CONNECT TO 1, 2, OR MORE MASSES. IF ONE, THEN MAKE THIS CHAIN AN ANTENNA
           # IF 2 OR MORE, MAKE THIS CHAIN A BOND
          
           # CAN TEST adjacent_vertices AGAINST maingraph, THEN DROP THE ADJACENCIES THAT ARE ALREADY ON THE CHAIN
           # THEN TEST WHICH expandedCore THE REMAINING RESULTS BELONG TO AND FILL THE SUMMARY TABLE THAT I DREW
           plugpoints <- adjacent_vertices(maingraph, endpoints)
           numrows <- length(plugpoints)
          
          
          if(numrows > 0) {
          
            # FOR EXAMPLE, THIS TESTS IF A VALUE IS IN AN EXPANDED CORE. IF IT IS, THEN CHECK WHICH ONE
            tab <- cbind(endpoints)
            tab <- as.data.frame(matrix(data=0, nrow=numrows, ncol=7))
            names(tab) <- c("PlugPoint", "Neigh1", "Neigh2", "Neigh3", "Neigh4", "Neigh5", "Neigh6")
            tab$PlugPoint <- as.integer(endpoints)
          
            # MAKING tab2 TO CONTAIN THE EXPANDED CORES CONNECTED TO
            tab2 <- tab
            tab2[,2:7] <- 0
        
            if(numrows > 0) {
              for(tabrow in 1:numrows) {
                N <- unlist(plugpoints[tabrow])
                connections <- as.integer(names(V(maingraph)[N]))
                numconnections <- length(connections)
              
                # BUILD AND APPEND TO tab, TO DEVELOP CONNECTIONS LIST
                tab[tabrow,2:(numconnections+1)] <- connections
                
                # NOW ASSESS WHICH OF THESE CONNECTIONS ARE TO EXPANDED CORES
                # tab2 IS A DUPLICATE OF tab BUT WITH THE CODES OF CORES CONNECTED TO (IF ANY)
                # ENTRIES OF 0 MEAN NO CONNECTION
                for(cell in 2:7) {
                  if(tab[tabrow,cell] != 0) {
                    # IS THE CONNECTION TO AN EXPANDED CORE?
                    findcoreflag <- tab[tabrow,cell] %in% voxelID[expandedCoreCode > 0]
                    if(findcoreflag) {
                      # IT DOES CONNECT TO AN EXPANDED CORE
                      tab2[tabrow, cell] <- expandedCoreCode[voxelID==tab[tabrow,cell]]
                    } else {
                      # NOT CONNECTED TO AN EXPANDED CORE
                      tab2[tabrow,cell] <- 0
                    } # END IF
                    
                  } # END IF
                } # END FOR: cell
              } # END FOR
              
              numcoreslinked <- length(unique(tab2[,2:7][tab2[,2:7]>0]))
               
              if(numcoreslinked == 1) {
                # IF 1 CORE, MAKE THIS CHAIN AN ANNTENNA
                # ALSO MUST ONLY HAVE CONNECTIONS BY ONE VOXEL, OTHERWISE IT IS A CIRCUIT ** NEED TO TEST THIS
                
                # SUM THE ROWS OF tab2 EXCLUDING THE FIRST COLUMN, EACH ROW EXCEPT ONE MUST BE ZEROS FOR THIS TO BE AN ANTENNA
                rowsums <- apply(tab2[,2:7], 1, sum)
                rowsums[rowsums>0] <- 1
                totsum <- sum(rowsums)
                if(totsum==1) {
                  writeobjects <- as.integer(names(V(tmpgph2[[chain]])))
                  # WRITE TO OBJECT
                  for(rep in 1:length(writeobjects)) {
                    morphCode[voxelID == writeobjects[rep]] <- 6
                  } # END FOR: rep
                } # END IF
              } # END IF
                            
              if(numcoreslinked > 1) {
                # IF MORE THAN 1 CORE, MAKE THIS CHAIN A BOND
                writeobjects <- as.integer(names(V(tmpgph2[[chain]])))
                # WRITE TO OBJECT
                for(rep in 1:length(writeobjects)) {
                  morphCode[voxelID == writeobjects[rep]] <- 7
                } # END FOR: rep
                
        
                # PROVIDE CORRECTOR HERE FOR SINGLE CIRCUIT VOXELS THAT CONNECT BOND TO SKIN, CONVERT THAT CIRCUIT TO BOND
                # REQUIRES IDENTIFYING THE BOND VOXELS WITH A CONNECTION TO CIRCUIT AND THEN CONVERTING THAT CIRCUIT TO BOND
                # LOGIC IS SIMILAR TO THE SINGLE ANTENNA VOXEL CORRECTOR
                for(rep in 1:length(writeobjects)) {
                  # REPEAT FOR EACH BOND VOXEL
                                   
                  loc <- which(voxelID==writeobjects[rep], arr.ind=TRUE)
                  # CHECK ALL SIX NEIGHBOURS AND COUNT (DO NOT EXTEND BEYOND THE ARRAY LIMITS)
                  ncounter <- 0
                  xlow <- 0
                  xhigh <- 0
                  ylow <- 0
                  yhigh <- 0
                  zlow <- 0
                  zhigh <- 0
                 
                  # ENSURE THAT I DO NOT LOOK AT A REFERENCE BEYOND THE ARRAY EXTENTS
                  if(loc[1] > 1) {
                    xlow <- morphCode[(loc[1]-1), loc[2], loc[3]]
                    if(xlow == 5) {
                      xlow <- 1
                    } else {
                      xlow <- 0
                    }
                  }
                  if(loc[1] < dimdatacube[1]) {
                    xhigh <- morphCode[(loc[1]+1), loc[2], loc[3]]
                    if(xhigh == 5) {
                      xhigh <- 1
                    } else {
                      xhigh <- 0
                    }
                  }
                  if(loc[2] > 1) {
                    ylow <- morphCode[loc[1], (loc[2]-1), loc[3]]
                    if(ylow == 5) {
                      ylow <- 1
                    } else {
                      ylow <- 0
                    }
                  }
                  if(loc[2] < dimdatacube[2]) {
                    yhigh <- morphCode[loc[1], (loc[2]+1), loc[3]]
                    if(yhigh == 5) {
                      yhigh <- 1
                    } else {
                      yhigh <- 0
                    }
                  }
                  if(loc[3] > 1) {
                    zlow <- morphCode[loc[1], loc[2], (loc[3]-1)]
                    if(zlow == 5) {
                      zlow <- 1
                    } else {
                      zlow <- 0
                    }
                  }
                  if(loc[3] < dimdatacube[3]) {
                    zhigh <- morphCode[loc[1], loc[2], (loc[3]+1)]
                    if(zhigh == 5) {
                      zhigh <- 1
                    } else {
                      zhigh <- 0
                    }
                  }
                  ncounter <- xlow + xhigh + ylow + yhigh + zlow + zhigh
                  
                  # IN ncounter IS 1, THEN THE APPROPRIATE CIRCUIT VOXEL NEEDS TO BE CONVETED TO BOND
                  # TAHT VOXEL IS INDICATED BY THE xlow, xhigh, ylow, yhigh, zlow, or zhigh VARIABLE BEING 1 AND NOT 0
                
                  if(ncounter==1) {
                    # AN UPDATE IS REQUIRED IN THE morphCode DATA LAYER
                    if(VERBOSE) {
                      cat("\n  An update of CIRCUIT to BOND is required")
                    }
                    if(xlow == 1) {
                      morphCode[(loc[1]-1), loc[2], loc[3]] <- 7
                    }
                    if(xhigh == 1) {
                      morphCode[(loc[1]+1), loc[2], loc[3]] <- 7
                    }
                    if(ylow == 1) {
                      morphCode[loc[1], (loc[2]-1), loc[3]] <- 7
                    }
                    if(yhigh == 1) {
                      morphCode[loc[1], (loc[2]+1), loc[3]] <- 7
                    }
                    if(zlow == 1) {
                      morphCode[loc[1], loc[2], (loc[3]-1)] <- 7
                    }
                    if(zhigh == 1) {
                      morphCode[loc[1], loc[2], (loc[3]+1)] <- 7
                    }
                  } # END IF: ncounter
                   
                } # END FOR: rep
              } # END IF
            } # END IF
          } # END IF
         } # END IF
       } # END FOR: chain
    
     } # END IF-ELSE
    } # END FOR: obj


  if(PLOT) {
    # THIS IS A PLOT OF THE ANTENNA VOXELS
    plotmorph <- morphCode
    plotmorph[plotmorph!=6] <- 0
    plotmorph[plotmorph==6] <- 5
    morphs <- morph3dprep(plotmorph, ORIG=TRUE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
    bgplot3d({
      plot.new()
      title(main = 'ANTENNA Voxels', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT THE MORPHOLOGY WITH ANTENNAE
  
  if(PLOT) {
    # THIS IS A PLOT OF THE BOND VOXELS
    plotmorph <- morphCode
    plotmorph[plotmorph!=7] <- 0
    plotmorph[plotmorph==7] <- 6
    morphs <- morph3dprep(plotmorph, ORIG=TRUE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
    bgplot3d({
      plot.new()
      title(main = 'BOND Voxels', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT THE MORPHOLOGY WITH ANTENNAE

  
  if(PLOT) {
    # THIS IS A PLOT OF THE CIRCUIT CONNECTOR VOXELS
    plotmorph <- morphCode
    plotmorph[plotmorph!=5] <- 0
    plotmorph[plotmorph==5] <- 4
    morphs <- morph3dprep(plotmorph, ORIG=TRUE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
    bgplot3d({
      plot.new()
      title(main = 'CIRCUIT Voxels', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT THE MORPHOLOGY WITH CIRCUITS
  
  
  
  # ------------ VOID-VOLUME (CODE = 8) ------------
 
 
  if(VERBOSE) {
    cat("\n\n  Identifying VOID-VOLUME voxels")
  }
  # NOTE THAT THE lrgdatacube2 IS A LARGER 3D ARRAY FILLED WITH -1
  
  # FLOOD FILL WITH -1 WHERE OUTSIDE THE SHAPE
  voidvolume <- lrgdatacube2
  for(row in 2:(dimdatacube[1]+1)) {
      for(col in 2:(dimdatacube[2]+1)) {
         for(z in 2:(dimdatacube[3]+1)) {
             if(voidvolume[row,col,z]==0 & ( voidvolume[(row+1),col,z] == -1 | voidvolume[(row-1),col,z] == -1 | voidvolume[row,(col+1),z] == -1 | voidvolume[row,(col-1),z] == -1 | voidvolume[row,col,(z+1)] == -1 | voidvolume[row,col,(z-1)] == -1 ) ) {
                 voidvolume[row,col,z] <- -1
             }
         }
      }
  }
  # SUBSET IT BACK TO PROPER SIZE
  voidvolume <- voidvolume[2:(dimdatacube[1]+1),2:(dimdatacube[2]+1),2:(dimdatacube[3]+1)]
  VOIDvolumeVoxels <- voxelID[voidvolume==0]
  # WRITE TO MORPHOLOGY OBJECT
  for(rep in 1:length(VOIDvolumeVoxels)) {
    morphCode[voxelID == VOIDvolumeVoxels[rep]] <- 8
  } # END FOR: rep
  # NOW THE TRUE OUTSIDE IS -1 AND VOID-VOLUME IS 0 AND THE OBJECT OF INTEREST IS 1
  # NEED TO CHECK THE ADJACENCY OF SKIN OR CONNECTORS TO 0 OR -1 TO DECIDE BETWEEN VOID AND SKIN RESPECTIVELY
  # THIS WILL BE CODE=9 FOR VOID IN THE NEXT SECTION
   
  if(PLOT) {
    plotmorph <- morphCode
    plotmorph[plotmorph!=8] <- 0
    plotmorph[plotmorph==8] <- 7
    morphs <- morph3dprep(plotmorph, ORIG=TRUE)
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
    bgplot3d({
      plot.new()
      title(main = 'VOID-VOLUME Voxels', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT
 
 
  # ------------ VOID (CODE = 9) ------------

  
  if(VERBOSE) {
    cat("\n\n  Identifying VOID-VOLUME voxels")
  }
  # INSERT INTO LARGER GRID TO DO THIS PROPERLY SO THAT WE DO NOT EXCEED EXTENTS
  lrgmorphCode <- array(data=1, dim=c(dimdatacube[1]+2, dimdatacube[2]+2, dimdatacube[3]+2))
  lrgmorphCode[2:(dimdatacube[1]+1),2:(dimdatacube[2]+1),2:(dimdatacube[3]+1)] <- morphCode
  
  # IF A VOID-VOLUME VOXEL NEIGHBOURS A VOXEL CODED AS SKIN, RECODE IT TO 9 FOR VOID
  # THIS CONVERTS SKIN TO VOID WHERE THE "EDGES" ARE INTERNAL
  for(row in 2:(dimdatacube[1]+1)) {
    for(col in 2:(dimdatacube[2]+1)) {
      for(z in 2:(dimdatacube[3]+1)) {
        if(lrgmorphCode[row,col,z]==8 & lrgmorphCode[(row+1),col,z] == 3 ) {
          lrgmorphCode[(row+1),col,z] <- 9
        }
        if(lrgmorphCode[row,col,z]==8 & lrgmorphCode[(row-1),col,z] == 3 ) {
          lrgmorphCode[(row-1),col,z] <- 9
        }
        if(lrgmorphCode[row,col,z]==8 & lrgmorphCode[row,(col+1),z] == 3 ) {
          lrgmorphCode[row,(col+1),z] <- 9
        }
        if(lrgmorphCode[row,col,z]==8 & lrgmorphCode[row,(col-1),z] == 3 ) {
          lrgmorphCode[row,(col-1),z] <- 9
        }
        if(lrgmorphCode[row,col,z]==8 & lrgmorphCode[row,col,(z+1)] == 3 ) {
          lrgmorphCode[row,col,(z+1)] <- 9
        }
        if(lrgmorphCode[row,col,z]==8 & lrgmorphCode[row,col,(z-1)] == 3 ) {
          lrgmorphCode[row,col,(z-1)] <- 9
        }
      }
    }
  }
 morphCode <- lrgmorphCode[2:(dimdatacube[1]+1),2:(dimdatacube[2]+1),2:(dimdatacube[3]+1)]
 
 
 if(PLOT) {
   plotmorph <- morphCode
   plotmorph[plotmorph!=9] <- 0
   plotmorph[plotmorph==9] <- 8
   morphs <- morph3dprep(plotmorph, ORIG=TRUE)
   open3d()
   morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
   bgplot3d({
     plot.new()
     title(main = 'VOID Voxels', line = 3)
   }) # END BGPLOT3D
 } # END IF: PLOT
 
 
 if(PLOT) {
   plotmorph <- morphCode
   plotmorph[plotmorph!=3] <- 0
   plotmorph[plotmorph==3] <- 2
   morphs <- morph3dprep(plotmorph, ORIG=TRUE)
   open3d()
   morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
   bgplot3d({
     plot.new()
     title(main = 'SKIN Voxels', line = 3)
   }) # END BGPLOT3D
 } # END IF: PLOT
  
 
 # ------------ FINISHING UP ------------
 
 
  # FINAL CODES ARE AS FOLLOWS:
  # 1 = OUTSIDE
  # 2 = MASS
  # 3 = SKIN
  # 4 = CRUMB
  # 5 = CIRCUIT
  # 6 = ANTENNA
  # 7 = BOND
  # 8 = VOID-VOLUME
  # 9 = VOID
  
  
  # PRODUCE SUMMARIES
  descriptions <- c("OUTSIDE","MASS","SKIN","CRUMB","CIRCUIT","ANTENNA","BOND","VOID-VOLUME", "VOID")
  summaries <- as.data.frame(matrix(NA, nrow=9, ncol=4))
  for(code in 1:9) {
    summaries[code,] <- c(code, descriptions[code], sum(morphCode==code), sum(morphCode==code)/max(voxelID)*100)
    names(summaries) <- c("Code", "Description", "NVoxels", "Percentage")
  }
  
  if(VERBOSE) {
    cat("\n\n")
    print(summaries)
    cat("\n\n")
  }
  
  # IF REQUESTED, CREATE A FINAL PLOT WITH ALL OF THE MORPHOLOGICAL ELEMENTS SHOWN
  if(FINALPLOT) {
    # THIS IS A PLOT OF THE FINAL 3D MORPHOLOGY CODES
    morphs <- morph3dprep(morphCode, ORIG=FALSE)
    # NEED TO FIX THE CELLID TO PLOT THE UNIQUE CLUSTER IDS
    open3d()
    morph3dplot(morphs, CELLID=PLOTIDS, LEGEND=FALSE, ORIGTRANSP=TRUE)
    bgplot3d({
      plot.new()
      title(main = '3D Morphology', line = 3)
    }) # END BGPLOT3D
  } # END IF: PLOT THE FINAL MORPHOLOGY
    
  # BUILD A LIST OBJECT TO RETURN THE INPUT, VOXEL IDS, MORPHOLOGY, AND OBJECT IDENTIFIERS, SUMMARY STATS
  outobj <- list(OriginalData=DATACUBE, Graph=decompgraph, VoxelIDs=voxelID, ObjectID=objectID, Morphology=morphCode, Cores=coreCode, ExpCores=expandedCoreCode, Summary=summaries, Egg=maingraph, Bgrnd=lrgdatacube2, VOIDvolume=voidvolume)

  # RETURN THE SEGMENTATION AND SUMMARIES
  return(outobj)
  
}
