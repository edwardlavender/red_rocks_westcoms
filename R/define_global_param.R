################################
################################
#### define_global_param.R

#### This code: 
# 1) Wipes the workspace, loads essential packages
# ... and defines parameters used by multiple scripts. 
# ... This is designed to be called at the start of every script. 

#### Steps preceding this code: 
# 1) WeStCOMS predictions have been generated.
# 2) WestCOMS predictions have been processed 
# ... in line with the requirements of the fvcom.tbx package:
# ... ... one folder for each variable 
# ... ... files labelled by date name 
# ... ... processed fields, such as current_speed, have been defined via fvcom.tbx. 


################################
################################
#### Global set up 

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(magrittr)
library(fvcom.tbx)


################################
################################
#### Define global parameters 

#### Define the root directory containing FVCOM predictions
# These have been processed in line with fvcom.tbx requirements. 
wc_root <- "/Volumes/Lacie_Share/Dima/FVCOM_variable_outputs/"
list.files(wc_root)


#### End of code. 
################################
################################