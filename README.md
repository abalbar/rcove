# rcove
 
The goal of **rcove** is to produce maps of dispersal area based on time series of ocean currents and simple life history characteristics in R. The current functionality supports in-water dispersal. 

# Installation

You can install **rcove** as a R package using two steps. 

**Step 1** Install the *R* package *devtools*

`if (!require("devtools")) install.packages("devtools") # to install`

**Step 2** Install *rcove*

`#install the package from *Github*`

`devtools::install_github("abalbar/rcove")`

`library(rcove) # load the library`

# Diagnostics

## rcove.R

Create map of dispersal area.

  -theta: Angle about which to rotate current velocities. Numeric
  
  -velocity: Velocity time series. Vector consisting of two columns u,v.
  
  -release_pts: Start locations for propagules. sf object
  
  -PD: Planktonic propagule duration of species of interest. Numeric
  
  -CPD: Competent propagule duration of species of interest. Numeric
  
  -land: Land to be subtracted from final concatenated polygons. sf polygon
  
  -proj: Projected crs in units m.
    -example: "+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs"
  

