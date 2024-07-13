R_grpnet_penalty <-
  function(znorm, penid, penone, pentwo, gamma, shrink){
    # grpnet_penalty.f90 translation to R
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2024-06-28
    
    shrink <- 1.0
    if(penid == 1L){
      shrink <- max(0, 1 - penone / znorm)
      shrink <- shrink / (1 + pentwo)
    } else if(penid == 2L){
      if(znorm <= gamma * penone * (1 + pentwo)){
        shrink <- max(0, 1 - penone / znorm)
        shrink <- shrink / (1 + pentwo - 1 / gamma)
      } else {
        shrink <- 1 / (1 + pentwo)
      }
    } else if(penid == 3L) {
      if(znorm <= (penone * (1 + pentwo) + penone)){
        shrink <- max(0, 1 - penone / znorm)
        shrink <- shrink / (1 + pentwo)
      } else if(znorm <= gamma * penone * (1 + pentwo)){
        shrink <- max(0, 1 - (penone / znorm) * (gamma / (gamma - 1)) )
        shrink <- shrink / (1 + pentwo - 1 / (gamma - 1) )
      } else {
        shrink <- 1 / (1 + pentwo)
      }
    }
    return(shrink)
    
  }