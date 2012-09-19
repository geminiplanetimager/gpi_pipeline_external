#  Set up environment variables for GPI IFS Software
#	bash version
#
#  See the GPI Pipeline Installation & Configuration Manual for instructions.
#
# HISTORY:
#       2010-01-25  Created. M. Perrin
#       2010-02-01  Added IDL_DIR, enclosing quotes - M. Perrin
#       2012-02-02  Split out the data from the software
#       2012-08-08  Simplification & reorganization
#       2012-08-10  Rewrite to conform to latest DRP setup - ds

# *** NOTE** This is the sample file for what the GPI environment
# variables should be. You should copy this file to some different
# filename, e.g.  export_GPI_mycomputername.csh, in your home dir 
# and edit it there. 
# Do NOT edit this file itself to customize it for your computer, 
# since this is under version control and thus your edits would just
# conflict with other users changing this same file too.
#
# Source this file from your .cshrc or .tcshrc (or profile)


#----- Required paths (editing required) --------------------
# GPI_DATA_ROOT is a helper path only.  If desired you can set 
# all paths independently.  
#  NOTE: If you are pulling all data from the vospace, do NOT 
#  put your queue and Reduced directories on vospace.
#  Make these local.
export GPI_DATA_ROOT="/path/to/your/GPI/data/"           	        # base dir for all data
export GPI_DRP_QUEUE_DIR="${GPI_DATA_ROOT}/queue/"      # where is the DRP Queue directory?
export GPI_RAW_DATA_DIR="${GPI_DATA_ROOT}/Detector/"    # where is raw data?
export GPI_REDUCED_DATA_DIR="${GPI_DATA_ROOT}/Reduced/" # where should we put reduced data?

#---- Optional paths (not genererally needed) -------
# these variables are optional - you may omit them if your 
# drp setup is standard
#export GPI_DRP_DIR="/path/to/your/GPI/code/pipeline"
#export GPI_DRP_CONFIG_DIR="${GPI_DRP_DIR}/config/"          # default config settings 
#export GPI_DRP_TEMPLATES_DIR="${GPI_DRP_DIR}/drf_templates"	# pipeline DRF template location

#export GPI_DRP_CALIBRATIONS_DIR="${GPI_REDUCED_DATA_DIR}/calibrations/"	# pipeline calibration location
#export GPI_DRP_LOG_DIR="${GPI_REDUCED_DATA_DIR}/logs/"	                    # default log dir

#---- DST install only (optional) -----
#export GPI_DST_DIR="${GPI_CODE_ROOT}/dst/"	 

#---------- make sure the startup scripts are in your $PATH   -----a
#           and the IDL code is in your $IDL_PATH
export PATH="${PATH}:/path/to/your/GPI/code/pipeline/scripts"
export IDL_PATH="${IDL_PATH}:+/path/to/your/GPI/code/"

