conflicted::conflicts_prefer(dplyr::filter)

##########################################################################
##                   Natal and adult dispersal evolution                ##
##                         Last version: 20/11/2025                     ##
##                                                                      ##
##    R version : 4.4.1 (2024-06-14)                                    ##
##    Koch and al.                                                      ##
##########################################################################

## Loading required packages ---------------------------------------------
library(here)
library(tidyverse)
library(dplyr)
library(MCMCglmm)
library(pedigree)
library(brms)
library(kinship2)
library(Matrix)
library(QGglmm)
library(performance)
library(BaSTA)
library(lmerTest)
library(ggeffects)
library(MASS)
library(emmeans)
library(Rmisc)
library(ggplot2)
library(ggridges)
library(MuMIn)
library(aod)

## Loading required datasets ----------------------------------------------
pedigree           <-      readRDS(here("01_Data/tbl_pedigree.rds"))
tbl_capt           <-      readRDS(here("01_Data/tbl_capture.rds"))
tbl_history        <-      readRDS(here("01_Data/tbl_history.rds"))
tbl_RE             <-      readRDS(here("01_Data/tbl_RE.rds"))
tbl_hapl           <-      readRDS(here("01_Data/tbl_RE.rds"))
tbl_Propention     <-      readRDS(here("01_Data/tbl_Propention.rds"))


tbl_dispJ <- 
  tbl_capt |> 
  filter(Stage == "SA")

tbl_dispA <-
  tbl_capt |> 
  filter(Stage == "A")
## I. Heritability ---------------------------------------------------------
# Estimates by animal model, using brms

## 1. Natal dispersal status -----------------------------------------------
# ----------------------------------------------------- 'Relatedness matrix'
Amat <- Matrix((2 * kinship(id     = pedigree[["animal"]],
                            momid  = pedigree[["ID_Mom"]],
                            dadid  = pedigree[["ID_Dad"]])
                # Matrix with the kinship coefficients between individuals
                ),
                sparse = TRUE)

# We only keep indiviuals used in animal models
Id_StatN <- unique((drop_na(tbl_dispJ, SEX, ID_Mum, PARTAN))[["animal"]])
Amat_StatN <- Amat[Id_StatN, Id_StatN]

# ------------------------------------------------------------------ 'Prior'
# We assumed fixed effects are normaly distributed
prior_StatN <- 
  set_prior("normal(0, sqrt(2))", class = "b")  

# ------------------------------------------------------------ 'Animal model'
# Fitted model
Mod_StatN <- 
  brm(
  formula = Propention ~ 
                1 + 
                SEX + 
                (1 | gr(animal, cov = Amat)) + 
                (1 | ID_Mum) + 
                (1 | PARTAN),
  data = drop_na(tbl_dispJ, 
                 SEX,
                 PARTAN, 
                 ID_Mum),
  data2 = list(Amat = Amat_StatN),                      
                            # Relatedness matrix
  family = bernoulli(link = "probit"), 
                            # Link function
  prior = prior_StatN,
  chains = 4,               # Number of differents chaines
  cores = 4,                # Number of cores (change according to used computer)
  iter = 3e5,               # Number of iterations
  warmup = 10000,           # Number of first iterations to exluded
  thin = 10                 # Interval between two saved iteration
)

# ---------------------------------------------- 'Model quality verification'
# Does adding a random effect for animal identity improve the model or not ? 
Mod_StatN_ctr <- 
  brm(
  formula = Propention ~ 
    1 + 
    SEX + 
    (1 | ID_Mum) + 
    (1 | PARTAN),
  data = drop_na(tbl_dispJ, 
                 SEX, 
                 PARTAN, 
                 ID_Mum),
  data2 = list(Amat = Amat_StatN),
  family = bernoulli(link = "probit"),
  prior = prior_StatN,
  chains = 4,
  cores = 4, 
  iter = 3e5,
  warmup = 10000,
  thin = 10
)

# LOOIC et ELPD
looic(Mod_StatN_ctr, verbose = TRUE)
looic(Mod_StatN, verbose = TRUE)
      
loo_compare(loo (Mod_StatN_ctr), loo(Mod_StatN))

# -------------------------------------------------- 'Variance decomposition'
tbl_var_StatN <- 
  purrr::map(VarCorr(Mod_StatN, 
                     summary = FALSE), 
             "sd") |> 
                          # Posterior distribution of random variables ("sd")
  map(\(x) { x^2 }) |> 
                          # Convert standards errors into variance
  tibble::as_tibble()

fixed_var_StatN <- 
  posterior_samples(Mod_StatN, 
                    pars = "^b_") |> 
                          # Posterior distribution of fixed variables ("sd")
  map(\(x) { x^2 })
                          # Convert standards errors into variance

tbl_var_StatN <- cbind(tbl_var_StatN, fixed_var_StatN)

Va_StatN <- tbl_var_StatN$animal          # Aditive genetic variance
CVa_StatN <- sqrt(Va_StatN)               # Evolvability
Vm_StatN <- tbl_var_StatN$ID_Mum          # Maternal variance
Ve_StatN <- tbl_var_StatN$PARTAN          # Cohort variance
Vr_StatN <- 1                             # Residual variance
V_int_StatN <- tbl_varStatN$b_Intercept
V_sex_StatN <- tbl_var_StatN$b_SEXM       # Sex effect variance

# ------------------------------------------------------------ 'Heritability'
# Latant scale
herit_StatN <- Va_StatN / (Va_StatN + Vm_StatN + Ve_StatN + Vr_StatN)
# saveRDS(herit_StatN,"StatN_h2_lat.rds")
# herit_StatN   <-      readRDS("03_Distribution/StatN_h2_lat.rds")

# Data scale
# Female
Post_F_StatN <- 
  pmap_dfr(list(mu = V_int_StatN, 
                                    # Prediction for females on the latent scale
                var.a = Va_StatN, 
                var.p = (Va_StatN + 
                           Ve_StatN + 
                           Vr_StatN + 
                           Vm_StatN - 
                           1)),     # - 1 because family = treeshold 
           QGparams, 
           model = "binom1.probit", 
           verbose = FALSE)

herit_Data_StatN <- Post_F_StatN[["h2.obs"]]
# saveRDS(herit_Data_StatN,"StatN_h2_don")
# herit_Data_StatN   <-      readRDS("03_Distribution/StatN_h2_don.rds")

# Male
Post_M_StatN <- 
  pmap_dfr(list(mu = V_int_StatN + V_sex_StatN, 
                                    # Prediction for males on the latent scale
                var.a = Va_StatN, 
                var.p = (Va_StatN + 
                           Ve_StatN + 
                           Vr_StatN + 
                           Vm_StatN - 
                           1)),     # - 1 because family = treeshold 
           QGparams, 
           model = "binom1.probit", 
           verbose = FALSE)

herit_Data_StatN <- Post_M_StatN[["h2.obs"]]



## 2. Adult dispersal status -----------------------------------------------
# ----------------------------------------------------- 'Relatedness matrix'
Amat <- Matrix((2 * kinship(id     = pedigree[["animal"]],
                            momid  = pedigree[["ID_Mom"]],
                            dadid  = pedigree[["ID_Dad"]])
                # Matrix with the kinship coefficients between individuals
                ),
                sparse = TRUE)

# We only keep indiviuals used in animal models
Id_StatA <- unique((drop_na(tbl_dispA, SEX, ID_Mum, PARTAN))[["animal"]])
Amat_StatA <- Amat[Id_StatA, Id_StatA]

# ------------------------------------------------------------------ 'Prior'
# We assumed fixed effects are normaly distributed
prior_StatA <- 
  set_prior("normal(0, sqrt(2))", class = "b")  

# ------------------------------------------------------------ 'Animal model'
# Fitted model
Mod_StatA <- 
  brm(
    formula = Propention ~ 
      1 + 
      Senescence_Mum + 
      (1 | gr(animal, cov = Amat)) + 
      (1 | ID_Mum) + 
      (1 | PARTAN),
    data = drop_na(tbl_dispA, 
                   Senescence_Mum,
                   PARTAN, 
                   ID_Mum),
    data2 = list(Amat = Amat_StatA),                      
    # Relatedness matrix
    family = bernoulli(link = "probit"), 
    # Link function
    control = list(adapt_delta = 0.9),
    prior = prior_StatA,
    chains = 4,               
    cores = 4,               
    iter = 3e5,             
    warmup = 10000,          
    thin = 10                 
  )

# ---------------------------------------------- 'Model quality verification'
# Does adding a random effect for animal identity improve the model or not ? 
Mod_StatA_ctr <- 
  brm(
    formula = Propention ~ 
      1 + 
      Senescence_Mum + 
      (1 | ID_Mum) + 
      (1 | PARTAN),
    data = drop_na(tbl_dispA, 
                   Senescence_Mum, 
                   PARTAN, 
                   ID_Mum),
    data2 = list(Amat = Amat_StatA),
    family = bernoulli(link = "probit"),
    control = list(adapt_delta = 0.9),
    prior = prior_StatA,
    chains = 4,
    cores = 4, 
    iter = 3e5,
    warmup = 10000,
    thin = 10
  )

# LOOIC et ELPD
looic(Mod_StatA_ctr, verbose = TRUE)
looic(Mod_StatA, verbose = TRUE)

loo_compare(loo (Mod_StatA_ctr), loo(Mod_StatA))

# -------------------------------------------------- 'Variance decomposition'
tbl_var_StatA <- 
  purrr::map(VarCorr(Mod_StatA, 
                     summary = FALSE), 
             "sd") |> 
  # Posterior distribution of random variables ("sd")
  map(\(x) { x^2 }) |> 
  # Convert standards errors into variance
  tibble::as_tibble()

fixed_var_StatA <- 
  posterior_samples(Mod_StatA, 
                    pars = "^b_") |> 
  # Posterior distribution of fixed variables ("sd")
  map(\(x) { x^2 })
# Convert standards errors into variance

tbl_var_StatA <- cbind(tbl_var_StatA, fixed_var_StatA)

Va_StatA      <-   tbl_var_StatA$animal          # Aditive genetic variance
CVa_StatA     <-   sqrt(Va_StatA)                # Evolvability
Vm_StatA      <-   tbl_var_StatA$ID_Mum          # Maternal variance
Ve_StatA      <-   tbl_var_StatA$PARTAN          # Cohort variance
Vr_StatA      <-   1                             # Residual variance
V_int_StatA   <-   tbl_var_StatA$b_Intercept
V_senM_StatA <-    tbl_var_StatA$b_Senescence_MumSenescent_m  # Senescent effect variance

# ------------------------------------------------------------ 'Heritability'
# Latant scale
herit_StatA <- Va_StatA / (Va_StatA + Vm_StatA + Ve_StatA + Vr_StatA)
# saveRDS(herit_StatA,"StatA_h2_lat.rds")
# herit_StatA   <-      readRDS("03_Distribution/StatA_h2_lat.rds")

# Data scale
# Young mother
Post_YM_StatA <- 
  pmap_dfr(list(mu = V_int_StatA, 
                                    # Prediction for young mother on the latent scale
                var.a = Va_StatA, 
                var.p = (Va_StatA + 
                           Ve_StatA + 
                           Vr_StatA + 
                           Vm_StatA - 
                           1)),     # - 1 because family = treeshold 
           QGparams, 
           model = "binom1.probit", 
           verbose = FALSE)

herit_Data_StatA <- Post_YM_StatA[["h2.obs"]]
# saveRDS(herit_Data_StatA,"StatNA_h2_don")
# herit_Data_StatA   <-      readRDS("03_Distribution/StatA_h2_don.rds")

# Senescent mother
Post_SM_StatA <- 
  pmap_dfr(list(mu = V_int_StatA + V_senM_StatA, 
                                    # Prediction for senescent mother on the latent scale
                var.a = Va_StatA, 
                var.p = (Va_StatA + 
                           Ve_StatA + 
                           Vr_StatA + 
                           Vm_StatA - 
                           1)),     # - 1 because family = treeshold 
           QGparams, 
           model = "binom1.probit", 
           verbose = FALSE)

herit_Data_StatA <- Post_SM_StatA[["h2.obs"]]
## 3. Natal dispersal distance ---------------------------------------------
# ------------------------------------------------------------------- 'Data'
tbl_DistN <- tbl_dispJ |> filter(Propention == "T")
# Data are not normaly distributed

# ----------------------------------------------------- 'Relatedness matrix'
Amat <- Matrix((2 * kinship(id     = pedigree[["animal"]],
                            momid  = pedigree[["ID_Mom"]],
                            dadid  = pedigree[["ID_Dad"]])
                # Matrix with the kinship coefficients between individuals
                ),
                sparse = TRUE)

# We only keep indiviuals used in animal models
Id_DistN <- unique((drop_na(tbl_DistN, SEX, ID_Mum, PARTAN))[["animal"]])
Amat_DistN <- Amat[Id_DistN, Id_DistN]

# ------------------------------------------------------------ 'Animal model'
# Fitted model
Mod_DistN <- 
  brm(
    formula = log(Distance) ~ 
      1 + 
      (1 | gr(animal, cov = Amat)) + 
      (1 | ID_Mum) + 
      (1 | PARTAN),
    data = drop_na(tbl_DistN, 
                   PARTAN, 
                   ID_Mum),
    data2 = list(Amat = Amat_DistN),                      
                                         # Relatedness matrix
    family = gaussian, 
                                         # Link function
    control = list(adapt_delta = 0.99),  # Avoid divergent transition after warmup
    chains = 4,               
    cores = 4,               
    iter = 2e5,             
    warmup = 20000,          
    thin = 20                 
  )

# ---------------------------------------------- 'Model quality verification'
# Does adding a random effect for animal identity improve the model or not ? 
Mod_DistN_ctr <- 
  brm(
    formula = Propention ~ 
      1 + 
      (1 | ID_Mum) + 
      (1 | PARTAN),
    data = drop_na(tbl_DistN,
                   PARTAN, 
                   ID_Mum),
    data2 = list(Amat = Amat_DistN),
    family = gaussian,
    control = list(adapt_delta = 0.99),
    chains = 4,
    cores = 4, 
    iter = 2e5,
    warmup = 20000,
    thin = 20
  )

# LOOIC et ELPD
looic(Mod_DistN_ctr, verbose = TRUE)
looic(Mod_DistN, verbose = TRUE)

loo_compare(loo (Mod_DistN_ctr), loo(Mod_DistN))

# -------------------------------------------------- 'Variance decomposition'
tbl_var_DistN <- 
  purrr::map(VarCorr(Mod_DistN, 
                     summary = FALSE), 
             "sd") |> 
                              # Posterior distribution of random variables ("sd")
  map(\(x) { x^2 }) |> 
                              # Convert standards errors into variance
  tibble::as_tibble()

fixed_var_DistN <- 
  posterior_samples(Mod_DistN, 
                    pars = "^b_") |> 
                              # Posterior distribution of fixed variables ("sd")
  map(\(x) { x^2 })
                              # Convert standards errors into variance

tbl_var_DistN <- cbind(tbl_var_DistN, fixed_var_DistN)

Va_DistN      <-   tbl_var_DistN$animal          # Aditive genetic variance
CVa_DistN     <-   sqrt(Va_DistN)                # Evolvability
Vm_DistN      <-   tbl_var_DistN$ID_Mum          # Maternal variance
Ve_DistN      <-   tbl_var_DistN$PARTAN          # Cohort variance
Vr_DistN      <-   tbl_var_DistN$residual_       # Residual variance
V_int_DistN   <-   tbl_var_DistN$b_Intercept

# ------------------------------------------------------------ 'Heritability'
# Latant scale
herit_DistN <- Va_DistN / (Va_DistN + Vm_DistN + Ve_DistN + Vr_DistN)
# saveRDS(herit_DistN,"DistN_h2_lat.rds")
# herit_DistN   <-      readRDS("03_Distribution/DistN_h2_lat.rds")

# Data scale
# Inverse link function
custom <- list(inv.link   = function(x) {exp(x)},
               # Inverse of the 'log' link function, i.e. exponential
               var.func   = function(x) {0},
               # No distribution variance
               d.inv.link = function(x) {exp(x)})
# Derived of the 'log' link function, i.e. exponential

# Posterior distribution
tbl_var_DistN <- data.frame(mu = as.vector(V_int_DistN),
                            va = as.vector(Va_DistN),
                            vp = as.vector(Va_DistN + Vm_DistN + Ve_DistN + Vr_DistN))

PostData_DistN <- do.call("rbind", 
                          apply(tbl_var_DistN, 1, function(row){
                            QGparams(mu = row[["mu"]],
                                     var.a = row[["va"]],
                                     var.p = row[["vp"]],
                                     custom.model = custom,
                                     verbose = FALSE)
                            }))

Herit_Data_DistN <- PostData_DistN[["h2.obs"]]
# saveRDS(Herit_Data_DistN,"DistN_h2_don.rds")

## 4. Natal dispersal distance ---------------------------------------------
# ------------------------------------------------------------------- 'Data'
tbl_DistA <- tbl_dispA |> filter(Propention == "T")
# Data are not normaly distributed

# ----------------------------------------------------- 'Relatedness matrix'
Amat <- Matrix((2 * kinship(id     = pedigree[["animal"]],
                            momid  = pedigree[["ID_Mom"]],
                            dadid  = pedigree[["ID_Dad"]])
                # Matrix with the kinship coefficients between individuals
                ),
                sparse = TRUE)

# We only keep indiviuals used in animal models
Id_DistA <- unique((drop_na(tbl_DistA, SEX, ID_Mum, PARTAN))[["animal"]])
Amat_DistA <- Amat[Id_DistA, Id_DistA]

# ------------------------------------------------------------ 'Animal model'
# Fitted model
Mod_DistA <- 
  brm(
    formula = log(Distance) ~ 
      1 + 
      (1 | gr(animal, cov = Amat)) + 
      (1 | ID_Mum) + 
      (1 | PARTAN),
    data = drop_na(tbl_DistA,
                   PARTAN, 
                   ID_Mum),
    data2 = list(Amat = Amat_DistA),                      
                                         # Relatedness matrix
    family = gaussian, 
                                         # Link function
    control = list(adapt_delta = 0.99),  # Avoid divergent transition after warmup
    chains = 4,               
    cores = 4,               
    iter = 2e5,             
    warmup = 20000,          
    thin = 20                 
  )

# ---------------------------------------------- 'Model quality verification'
# Does adding a random effect for animal identity improve the model or not ? 
Mod_DistA_ctr <- 
  brm(
    formula = Propention ~ 
      1 + 
      (1 | ID_Mum) + 
      (1 | PARTAN),
    data = drop_na(tbl_DistA,
                   PARTAN, 
                   ID_Mum),
    data2 = list(Amat = Amat_DistA),
    family = gaussian,
    control = list(adapt_delta = 0.99),
    chains = 4,
    cores = 4, 
    iter = 2e5,
    warmup = 20000,
    thin = 20
  )

# LOOIC et ELPD
looic(Mod_DistA_ctr, verbose = TRUE)
looic(Mod_DistA, verbose = TRUE)

loo_compare(loo (Mod_DistA_ctr), loo(Mod_DistA))

# -------------------------------------------------- 'Variance decomposition'
tbl_var_DistA <- 
  purrr::map(VarCorr(Mod_DistA, 
                     summary = FALSE), 
             "sd") |> 
  # Posterior distribution of random variables ("sd")
  map(\(x) { x^2 }) |> 
  # Convert standards errors into variance
  tibble::as_tibble()

fixed_var_DistA <- 
  posterior_samples(Mod_DistA, 
                    pars = "^b_") |> 
  # Posterior distribution of fixed variables ("sd")
  map(\(x) { x^2 })
# Convert standards errors into variance

tbl_var_DistA <- cbind(tbl_var_DistA, fixed_var_DistA)

Va_DistA      <-   tbl_var_DistA$animal          # Aditive genetic variance
CVa_DistA     <-   sqrt(Va_DistA)                # Evolvability
Vm_DistA      <-   tbl_var_DistA$ID_Mum          # Maternal variance
Ve_DistA      <-   tbl_var_DistA$PARTAN          # Cohort variance
Vr_DistA      <-   tbl_var_DistA$residual_       # Residual variance
V_int_DistA   <-   tbl_var_DistA$b_Intercept

# ------------------------------------------------------------ 'Heritability'
# Latant scale
herit_DistA <- Va_DistA / (Va_DistA + Vm_DistA + Ve_DistA + Vr_Dist)A
# saveRDS(herit_DistA,"DistA_h2_lat.rds")
# herit_DistA   <-      readRDS("03_Distribution/DistA_h2_lat.rds")

# Data scale
# Inverse link function
custom <- list(inv.link   = function(x) {exp(x)},
               # Inverse of the 'log' link function, i.e. exponential
               var.func   = function(x) {0},
               # No distribution variance
               d.inv.link = function(x) {exp(x)})
# Derived of the 'log' link function, i.e. exponential

# Posterior distribution
tbl_var_DistA <- data.frame(mu = as.vector(V_int_DistA),
                            va = as.vector(Va_DistA),
                            vp = as.vector(Va_DistA + Vm_DistA + Ve_DistA + Vr_DistA))

PostData_DistN <- do.call("rbind", 
                          apply(tbl_var_DistA, 1, function(row){
                            QGparams(mu = row[["mu"]],
                                     var.a = row[["va"]],
                                     var.p = row[["vp"]],
                                     custom.model = custom,
                                     verbose = FALSE)
                          }))

Herit_Data_DistA <- PostData_DistA[["h2.obs"]]
# saveRDS(Herit_Data_DistA,"DistA_h2_don.rds")









## II. Survival ------------------------------------------------------------
# Capture probabilities allowed to change for each year
pCapt <- c(1999, 
           2000, 
           2001, 
           2002, 
           2003, 
           2004, 
           2005, 
           2006, 
           2007, 
           2008, 
           2009, 
           2010, 
           2011, 
           2012, 
           2013, 
           2014, 
           2015, 
           2016, 
           2017, 
           2018, 
           2019, 
           2020, 
           2021, 
           2022)

# Fitted model:
# Selected covariables
LCov <- c("F_J_D", "F_S_D", "F_J_P", "F_S_P", "M_J_D", "M_S_D", "M_J_P", "M_S_P")

tbl_history <- 
  tbl_history |>
  select(
    - Distance
  )

Mod_survival <- basta(tbl_history, 
                      studyStart = 1999, 
                      studyEnd = 2022,
                      model = "LO", 
                      shape = "simple",
                      recaptTrans = pCapt,
                      covarsStruct = "fused",   
                      formulaMort = LCov_Ls,        
                      parallel = TRUE,
                      ncpus = 8,
                      nsim = 8,
                      niter = 400000,            
                      burnin = 5000,            
                      thinning = 40,
                      updateJumps = TRUE, 
                      negSenescence = TRUE)


## III. Annual reproductive success ----------------------------------------
# i.e. for adult dispersal
# ------------------------------------------------------------------- 'Data'
tbl_RE <- 
  tbl_RE |>
  filter(Id %in% tbl_dispA$animal) |>
  left_join(
    (tbl_dispA |> select(animal,
                         Distance,
                         Propention,
                         AN)),
    join_by("Id" == "animal", "Year" == "AN")
  )

## 1. Clutch size ----------------------------------------------------------
# ~ Dispersal status
Mod_CS_Stat <- brm(Clutch_S ~ Sex + Propention,   
                   data = tbl_RE,
                   family = Gamma(link = "log"),
                   chains = 4,
                   cores = 4,
                   iter = 3000,
                   warmup = 1000,
                   thin = 5)

# ~ Dispersal distance
Mod_CS_Dist <- brm(Clutch_S ~ Sex + log(Distance),   
                   data = tbl_RE,
                   family = Gamma(link = "log"),
                   chains = 4,
                   cores = 4,
                   iter = 3000,
                   warmup = 1000,
                   thin = 5)






## 2. Clutch weight ---------------------------------------------------------
# ~ Dispersal status
Mod_CW_Stat <- brm(Clutch_W ~ Propention,   
                   data = tbl_RE,
                   family = Gamma(link = "log"),
                   chains = 4,
                   cores = 4,
                   iter = 3000,
                   warmup = 1000,
                   thin = 5)

# ~ Dispersal distance
Mod_CW_Dist <- brm(Clutch_W ~ log(Distance),   
                   data = tbl_RE,
                   family = Gamma(link = "log"),
                   chains = 4,
                   cores = 4,
                   iter = 3000,
                   warmup = 1000,
                   thin = 5)






## 3. Proportion of variables offsprings ------------------------------------
# ~ Dispersal status
Mod_VO_Stat <- brm(V | trials(V + D + A + E) ~ Propention,   # 
                   data = tbl_RE,
                   family = binomial(link = "logit"),
                   chains = 4,
                   cores = 4,
                   iter = 3000,
                   warmup = 1000,
                   thin = 5)

# ~ Dispersal distance
Mod_VO_Dist <- brm(V | trials(V + D + A + E) ~ log(Distance),   # 
                   data = x,
                   family = binomial(link = "logit"),
                   chains = 4,
                   cores = 4,
                   iter = 3000,
                   warmup = 1000,
                   thin = 5)


## IV. First reproduction event --------------------------------------------
# ------------------------------------------------------------------- 'Data'
tbl_RE <-
  tbl_RE |> 
  filter(Id %in% tbl_juv$ID, 
         Year == Year_FR)

## 1. Age at first reproduction --------------------------------------------
# ~ Statut
Mod_Age_FR_Stat <- brm(Age_FR ~ 1,   
                       data = tbl_RE,
                       family = poisson(link = "log"),
                       chains = 4,
                       cores = 4,
                       iter = 3000,
                       warmup = 1000,
                       thin = 5)

# ~ Distance
Mod_Age_FR_Dist <- brm(Age_FR ~ 1,  
                       data = tbl_RE,
                       family = poisson(link = "log"),
                       chains = 4,
                       cores = 4,
                       iter = 3000,
                       warmup = 1000,
                       thin = 5)

## 2. Clutch size at first reproduction ------------------------------------
# ~ Status
Mod_CS_FR_Stat <- brm(Clutch_S ~ Sex,   
                      data = tbl_RE,
                      family = Gamma(link = "log"),
                      chains = 4,
                      cores = 4,
                      iter = 3000,
                      warmup = 1000,
                      thin = 5)

# ~ Distance
Mod_CS_FR_Dist <- brm(Clutch_S ~  Sex,
                      family = Gamma(link = "log"),
                      chains = 4,
                      cores = 4,
                      iter = 3000,
                      warmup = 1000,
                      thin = 5)
## 3. Clutch weight at first reproduction -----------------------------------
# ~ Status
Mod_CW_FR_Stat <- brm(Clutch_W ~ 1,   
                      data = tbl_RE,
                      family = Gamma(link = "log"),
                      chains = 6,
                      cores = 4,
                      iter = 3000,
                      warmup = 1000,
                      thin = 5)

# ~ Distance
Mod_CW_FR_Dist <- brm(Clutch_W ~  1,
                      family = Gamma(link = "log"),
                      chains = 4,
                      cores = 4,
                      iter = 3000,
                      warmup = 1000,
                      thin = 5)


## 4. Proportion of viable offsprings -------------------------------------
# ~ Status
Mod_VO_FR_Stat <- brm(V | trials(V + D + A + E) ~ Propention,
                      data = tbl_RE,
                      family = binomial(link = "logit"),
                      chains = 4,
                      cores = 4,
                      iter = 3000,
                      warmup = 1000,
                      thin = 5)

# ~ Distance
Mod_VO_FR_Dist <- brm(V | trials(V + D + A + E) ~ 1,
                      data = tbl_RE,
                      family = binomial(link = "logit"),
                      chains = 4,
                      cores = 4,
                      iter = 3000,
                      warmup = 1000,
                      thin = 5)

## V. Gene dropping --------------------------------------------------------
## 1. Modelisation ---------------------------------------------------------
# According to Reid et al. 2019
# ------------------------------------------------------------------- 'Data'
tbl_capt <- rbind(tbl_dispA, tbl_dispJ)

# ------------------------------------------------------------- 'Parameters'
set_cohort <- 1999       # Focal cohort
N_iter <- 10000          # Number of iteration used in gene dropping
Yafter <- 15             # Number of year considered during gene dropping

# ------------------------------------------------------------- 'Formatting'
tbl_lab <- pedigree |>
  filter(str_detect(animal, "^[N]")) |>
  mutate(
    PARTAN = str_extract(animal, "\\d{4}"),        # Adding birth year
    Pr_Contact = NA
  )

tbl_sauv <- pedigree |>
  filter(str_detect(animal, "^[C]"))|>
  mutate(
    PARTAN = NA,
    Pr_Contact = str_extract(animal, "\\d{4}"),    # Adding first encounter year
  )

tbl_cryo <- rbind(tbl_lab, tbl_sauv)

# -------------------------- 'Saving cohorts subsequent to the focal cohort'
tbl_ancestor <- 
  tbl_cryo |>
  filter(
    PARTAN > set_cohort
  )

# --------------------------------------------------------------- 'Prunning'
# Keeping only the individuals that have reach sexual maturity - 
# i.e. had the opportunity for reproduction
tbl_rox <- 
  tbl_capt |>
  filter(
    PARTAN == set_cohort,
    Stage == "A"
  ) |>
  select(animal)

tbl_rox <- tbl_cryo |>
  filter(
    PARTAN == set_cohort,
    animal %in% tbl_rox$animal
  )

tbl_ancestor <- rbind(tbl_rox, tbl_ancestor)

# --------------------------------------------------------------- 'Founders'
# Individuals born during the focal cohort are regarded as founders ;
# Replacing their parents by 'NA'
tbl_ancestor$ID_Mom <- 
  ifelse(
  tbl_ancestor$PARTAN > set_cohort, 
  tbl_ancestor$ID_Mom, 
  NA
)

tbl_ancestor$ID_Dad <- 
  ifelse(
  tbl_ancestor$PARTAN > set_cohort, 
  tbl_ancestor$ID_Dad, 
  NA
)

# Individual found for the first time in the wild are considered as founders
tbl_lazy <- tbl_cryo |>
  filter(Pr_Contact == (set_cohort + 1))

tbl_ancestor <- rbind(tbl_lazy, tbl_ancestor)

# ------------------------------------------------------- 'Last corrections'
tbl_ancestor$PARTAN[is.na(tbl_ancestor$PARTAN)] <- 1000
tbl_ancestor <- prepPed(tbl_ancestor)
tbl_gdrop <- tbl_ancestor[order(tbl_ancestor$PARTAN, tbl_ancestor$animal),]

# Founders
Ind_cohort <- sum(tbl_gdrop$PARTAN == set_cohort)
Ind_C <- sum(tbl_gdrop$PARTAN == 1000)

Nb_Founders <- Ind_cohort + Ind_C

# -------------------------------------------------------- 'Founders matrix'
# Without inbreeding betwin founders, the matrix must be nul
l <- rep(0, times = Nb_Founders)
Fond_mat <- matrix(NA, nrow = length(l), ncol = N_iter)

for (i in 1:N_iter) {
  Fond_out <- rbinom(length(l), 1, l)
  Fond_mat[,i] <- Fond_out
}

# --------------------------------------------- 'Vector of all reaslisation'
Vec_F <- as.vector(t(Fond_mat))

# --------------------------------------------------------------- 'Function'
# Again, by Reid et al. 2019
geneDrop <- function(pedigree, N, f.outcome,
                     parallel = FALSE, ncores = getOption("mc.cores", 2L), ...){
  nPed <- numPed(pedigree)
  n <- nrow(pedigree)
  dfounders <- which(nPed[, 2] == -998)
  sfounders <- which(nPed[, 3] == -998)
  dalleles <- salleles <- vector("integer", length = n) 
  dalleles[dfounders] <- as.integer(dfounders)
  salleles[sfounders] <- as.integer(-sfounders)
  Ndalleles <- rep(dalleles, each = N)
  Nsalleles <- rep(salleles, each = N)
  f.outcome.all <- c(f.outcome, rep(0, times = N*nrow(pedigree) - length(f.outcome)))
  Nsalleles <- ifelse(f.outcome.all == 1, Ndalleles, Nsalleles)
  
  Cout <- .C("genedrop",
             as.integer(Ndalleles),
             as.integer(Nsalleles),
             as.integer(N),
             as.integer(n),
             as.integer(nPed[, 2] - 1),
             as.integer(nPed[, 3] - 1))
  
  return(list(IDs = pedigree[, 1],
              maternal = matrix(Cout[[1]], ncol = N, byrow = TRUE),
              paternal = matrix(Cout[[2]], ncol = N, byrow = TRUE),
              numericPedigree = nPed,
              dalleles.check = dalleles,
              salleles.check = salleles,
              Ndalleles.check = Ndalleles,
              Nsalleles.check = Nsalleles,
              max.dalleles = max(dalleles) # La priorité est données aux lignées maternelles, qui sont tjs complètes
  ))
}

# ----------------------------------------------------------- 'Modelisation'
Mod_gD <- geneDrop(tbl_gdrop[,1:3], 
                   N = N_iter, 
                   f.outcome = Vec_F)
## 2. Information extraction -----------------------------------------------
# ----------------------------------------------------------------- 'Output'
# Maternal allele for each individual
tbl_Mum <- as.data.frame(Mod_gD$maternal)
tbl_Mum$Id <- tbl_Ped$Id

# Paternal allele for each individual
tbl_dad <- as.data.frame(Mod_gDrop$paternal)
tbl_dad$animal <- tbl_gdrop$animal

# Founders allele
ID_allele <- 1:Mod_gDrop$max.dalleles

## Identité des allèles des individus de la cohort d'intérêt
tbl_gataca <-
  tbl_gdrop |>
  filter(PARTAN <= set_cohort)

tbl_Gat <-
  tbl_gataca |>
  filter(PARTAN == set_cohort) |>
  select(animal,
         allele)
# ------------------------------------------------- 'Founder's contribution'
# ----------------------------------------------------------------- 'Output'
# By cohort, for each known individual
# Again, from Reid et al. 2019

# Complet pedigree
for(i in 1:G_after){
  lezardeaux <- subset(tbl_Ped, tbl_Ped$GenB == (Set_gen + i))
  assign(paste("lezardeaux.", i, sep = ""), lezardeaux)
}

# Contribution to each cohort
for(i in 1:G_after){
  
  tbl_lezardeaux <- get(paste0("lezardeaux.", i))
  
  # Maternal contribution
  mum_contrib <- merge(tbl_Mum, subset(tbl_lezardeaux, select = "Id"))
  
  # Paternal contribution
  dad_contrib <- merge(tbl_Dad, subset(tbl_lezardeaux, select = "Id"))
  
  # Combine both
  all_contrib <- rbind(mum_contrib, dad_contrib)
  
  # Compute to total number of each kind of allele present in each focal cohort,
  # for each iterations
  allele_counts <- sapply(all_contrib[-1], function(x) table(x), simplify = F)
  
  # Create a dataframe with the results
  tbl_iter <- data.frame(c(-Mod_gD_P$max.dalleles:-1, 1:Mod_gD_P$max.dalleles))
  colnames(tbl_iter) <- "allele"
  
  # Add allele count to the results dataframe
  tbl_res <- 
    allele_counts |> 
    map(tibble::enframe, name = "allele", value = "count") |> 
    list_rbind(names_to = "Iter")
  
  tbl_res <- 
    tbl_res |> 
    mutate(Iter = str_replace(Iter, "V", "X")) |> 
    pivot_wider(names_from = Iter, values_from = count)
  
  tbl_res <- merge(tbl_res, tbl_iter)
  
  # Replace all NAs with zeros
  tbl_res[is.na(tbl_res)] <- 0
  
  # Keeping only positive alleles of individuals born in the focal cohort
  tbl_res <- merge(tbl_res, tbl_Gat)
  tbl_res$gen <- rep(Set_gen, times = nrow(tbl_res))
  
  tbl_analys <- tbl_res
  
  # Absolute mean and variance
  tbl_analys[["Contrib_an"]] <- apply(tbl_res |>
                                        select(
                                          - Id, 
                                          - allele,
                                        ), 1, mean)
  
  # Frequency of zero copies of each allele present in each focal individual
  tbl_analys[["Nb_0"]] <- apply(subset(tbl_res,
                                       select = -c(Id, 
                                                   allele)), 
                                1, 
                                function(x) sum(x == 0))
  
  tbl_analys[["P_ext"]] <- tbl_analys$Nb_0 / N_iter
  tbl_analys$An_capt <- rep(i, times = nrow(tbl_analys))
  
  tbl_res <- merge(tbl_res, tbl_analys)
  tbl_res <- tbl_res |> select(Id, gen, Contrib_an, P_ext, An_capt)
  
  assign(paste("tbl_res_", i, sep = ""), tbl_res)
  
  write.csv(tbl_res, paste("tbl_b_", set_cohort, "_Y", i, ".csv", sep = ""))
}








# ------------------------------------------------------------------- 'Data'

## 3. Selection gradients - natal dispersal --------------------------------
# ------------------------------------------------------------------- 'Data'
# Phenotype
tbl_GS <-
  tbl_hapl |>
  left_join((tbl_capt |>
               select(animal,
                      AN,
                      Age,
                      Distance,
                      Propention,
                      SEX)), by = c("animal" = "animal")) |>
  filter(!is.na(Distance),
         Age == 1,
         cohort < 2011) |>
  select(animal,
         cohort,
         Contrib,
         L_surv,
         AN,
         Age,
         Distance,
         Propention,
         SEX)

# Binary data
tbl_binomiale <- tbl_GS |>
  mutate(
    bin_C = case_when(
      Contrib == 0 ~ 0,
      Contrib != 0 ~ 1
    ),
    bin_LS = case_when(
      L_surv == 0 ~ 0,
      L_surv != 0 ~ 1
    ) 
  ) |>
  select(
    bin_C,
    bin_LS,
    SEX,
    Distance,
    Propention,
    cohort
  )

# Truncated data (strictly positive) 
tbl_truncated <- tbl_GS |>
  filter(
    Contrib != 0,
    L_surv != 0
  ) |>
  select(
    Contrib,
    L_surv,
    SEX,
    Distance,
    Propention,
    cohort
  )

# ------------------------------------------------------------ 'Binary data'
# Expected genetic contribution
# ~ Status
Mod_BCS <- brm(bin_C ~ Propention + SEX +
                 (1 | cohort), 
               data = tbl_binomiale,
               family = bernoulli(link = "logit"),
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)
# ~ Distance
Mod_BCD <- brm(bin_C ~ log(Distance) + SEX +
                 (1 | cohort), 
               data = tbl_binomiale,
               family = bernoulli(link = "logit"),
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)

# Lineage longevity
# ~ Status
Mod_BCS <- brm(bin_C ~ Propention + SEX +
                 (1 | cohort), 
               data = tbl_binomiale,
               family = bernoulli(link = "logit"),
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)
# ~ Distance
Mod_BCD <- brm(bin_C ~ log(Distance) + SEX +
                 (1 | cohort), 
               data = tbl_binomiale,
               family = bernoulli(link = "logit"),
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)

# --------------------------------------------------------- 'Truncated data'
# Expected genetic contribution
# ~ Status
Mod_TCS <- brm(Contrib ~ 1 +
                 (1 | cohort),
               family = Gamma(link = "log"),
               data = tbl_truncated,
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)
# ~ Distance
Mod_TCD <- brm(Contrib ~ 1 +
                 (1 | cohort),
               family = Gamma(link = "log"),
               data = tbl_truncated,
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)

# Lineage longevity
# ~ Status
Mod_TLS <- brm(L_surv ~ Propention +
                 (1 | cohort), 
               data = tbl_truncated,
               family = Gamma(link = "log"),
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)
# ~ Distance
Mod_TLD <- brm(L_surv ~ 1 +
                 (1 | cohort),
               family = Gamma(link = "log"),
               data = tbl_truncated,
               chains = 4,
               cores = 4,
               iter = 10000,
               warmup = 1000,
               thin = 10)
## 4. Selection gradients - adult dispersal --------------------------------
# ------------------------------------------------------------------- 'Data'
# Phenotype
tbl_GS <-
  tbl_hapl |>
  left_join((tbl_capt |>
               select(animal,
                      SEX) |> 
               unique()), by = c("animal" = "animal")) |>
  left_join(tbl_Propention, by = c("animal" = "animal")) |>
  filter(!is.na(BLUP),
         cohort < 2011) |>
  select(animal,
         cohort,
         Contrib,
         L_surv,
         SEX,
         BLUP, 
         Iter)

# Binary data
tbl_binomiale <- tbl_GS |>
  mutate(
    bin_C = case_when(
      Contrib == 0 ~ 0,
      Contrib != 0 ~ 1
    ),
    bin_LS = case_when(
      L_surv == 0 ~ 0,
      L_surv != 0 ~ 1
    ) 
  ) |>
  select(
    bin_C,
    bin_LS,
    SEX,
    BLUP,
    cohort,
    Iter
  )

List_binomiale <- split(tbl_binomiale, tbl_binomiale$Iter)

# Truncated data (strictly positive) 
tbl_truncated <- tbl_GS |>
  filter(
    Contrib != 0,
    L_surv != 0
  ) |>
  select(
    Contrib,
    L_surv,
    SEX,
    BLUP,
    cohort,
    Iter
  ) 

List_truncated <- split(tbl_truncated, tbl_truncated$Iter)

# ------------------------------------------------------------ 'Binary data'
# Expected genetic contribution
# ~ Propention to adult dispersal
Mod_BCM <- brm_multiple(bin_C ~ SEX, 
                        data = List_binomiale,
                        family = bernoulli(link = "logit"),
                        chains = 1,
                        cores = 4)

# --------------------------------------------------------- 'Truncated data'
# Expected genetic contribution
# ~ Propention to adult dispersal
Mod_TCM <- brm_multiple(Contrib ~ 1, 
                        data = List_truncated,
                        family = Gamma(link = "log"),
                        chains = 1)

# Lineage longevity
# ~ Propention to adult dispersal
Mod_TLM <- brm_multiple(L_surv ~ 1, 
                        data = List_truncated,
                        family = Gamma(link = "log"),
                        chains = 4)


############################################################################
# Adult dispersal information

# -------------------------------------------------- 'Dispersal probability'
Id_disp <- 
  tbl_capt |> 
  filter(
    Propention == "T"
  ) |> 
  select(animal) |>
  unique()

n_disp <- nrow(Id_disp)

Id_resi <- 
  tbl_capt |> 
  filter(
    Propention == "F",
    !(animal %in% Id_disp$animal)
  ) |> 
  select(animal) |>
  unique()

n_resi <- nrow(Id_resi)

n_disp / (n_resi + n_disp) # 0.362

# ----------------------------------------------------- 'Dispersal distance'
tbl_disp <-
  tbl_capt |>
  filter(Propention == "T")

min(tbl_disp$Distance)     # 20.1 
max(tbl_disp$Distance)     # 158.162
mean(tbl_disp$Distance)    # 42.767
sd(tbl_disp$Distance) / sqrt(length(tbl_disp$Distance)) # 1.006

# -------------------------------------------------------------- 'Sex & Age'
tbl_disp <-
  tbl_capt |>
  filter(!is.na(Age),
         !is.na(SEX),
         !is.na(Propention)) |>
  summarise(
    n_F = sum(Propention == "F", na.rm = TRUE),
    n_T = sum(Propention == "T", na.rm = TRUE),
    .by = c(Age, SEX)
  ) |>
  mutate(
    total = n_F + n_T
  )

Mod_AD_1 <- glm(cbind(n_T, total) ~ Age + SEX,
                family = binomial,
                tbl_disp)
# No significant effect of Sex

Mod_AD_1 <- glm(cbind(n_T, total) ~ Age,
                family = binomial,
                tbl_disp)









