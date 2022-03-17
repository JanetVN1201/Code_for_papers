# This is a script for testing fMRI cortical surface models with INLA
# Things that should be changed ----
data_dir <- "" # The unzipped folder
PARDISO_license <- "" # PARDISO license
workbench_dir <- "" # Connectome workbench installation
hem <- "left" # This should be "left" or "right"
# Things that should stay the same ----

# >> Read in data files and format them correctly -----
file_names <- list.files(data_dir, full.names = T)
cifti_files <- grep("dtseries.nii", file_names, value = T)
left_surface <- grep("L.midthickness.32k_fs_LR.surf.gii",
                     file_names, value = T)
right_surface <- grep("R.midthickness.32k_fs_LR.surf.gii",
                      file_names, value = T)
task_names <- c("cue","left_foot","left_hand",
                "right_foot","right_hand","tongue")
onset_list <- sapply(c("LR","RL"), function(session_name) {
    session_files <- grep(paste0(session_name,".txt"), file_names, value = T)
    session_files <- grep("Movement", session_files, value = T, invert = T)
    session_onsets <- lapply(session_files,read.table,header = F)
    names(session_onsets) <- task_names
    contralateral <- grep(hem, task_names)
    session_onsets <- session_onsets[-contralateral]
    return(session_onsets)
}, simplify = F)

motion_list <- sapply(grep("Movement", file_names, value = T),
                      function(x) as.matrix(read.table(x,header = F)),
                      simplify = F)
names(motion_list) <- c("LR","RL")

# >> Load the libraries and run the fMRI analysis ----
library(INLA)
inla.setOption(pardiso.license = PARDISO_license)
inla.setOption(inla.mode = "experimental")
#inla.setOption(
library(ciftiTools) # on CRAN
ciftiTools.setOption('wb_path',workbench_dir)
# via github: 
#devtools::install_github("mandymejia/BayesfMRI",ref = "master")
library(BayesfMRI) 

start_time <- proc.time()[3]
bayesianGLM_result3 <-
  BayesGLM_cifti(
    cifti_fname = cifti_files,
    surfL_fname = left_surface,
    surfR_fname = right_surface,
    brainstructures = hem,
    design = NULL,
    onsets = onset_list,
    TR = 0.72,
    nuisance = motion_list,
    nuisance_include = c("drift", "dHRF"),
    scale_BOLD = TRUE,
    scale_design = TRUE,
    GLM_method = "Bayesian",
    ar_order = 6, # Prewhitening (this takes some time, and can be disabled by setting this to 0 for speed comparisons)
    ar_smooth = 6, # Prewhitening
    session_names = c("LR", "RL"), # This will perform the multisession Bayesian GLM
    ### PARAMETER TO PLAY WITH -> determines matrix size, number of nodes in mesh?!
    resamp_res = 5000, # Resampling value (we use 5000 in the Bayesian GLM validation paper)
    verbose = TRUE,
    outfile = NULL,
    return_INLA_result = TRUE,
    avg_sessions = TRUE,
    trim_INLA = FALSE, # This will trim the INLA object so that it only returns objects necessary for group analysis.
    keep = TRUE,
    twostage = TRUE
    )
bayesianGLM_result3$total_time <- proc.time()[3] - start_time

saveRDS(bayesianGLM_result3, file = file.path(data_dir,'bayesianGLM_result3.rds'))

#####################
result1 <- readRDS("~/Documents/bayesianGLM_result3.rds")
plot(result1)
