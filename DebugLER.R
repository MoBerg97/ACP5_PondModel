# library(gotmtools)
# library(LakeEnsemblR)
# getwd()
# 
# detach("package:LakeEnsemblR", unload=T)
# remove.packages("LakeEnsemblR")
# remotes::install_github("aemon-j/LakeEnsemblR")
# library(LakeEnsemblR)

r <- unclass(lsf.str(envir = asNamespace("LakeEnsemblR"), all = T))
# create functions in the Global Env. with the same name
for(name in r) eval(parse(text=paste0(name, '<-LakeEnsemblR:::', name)))

# cali_ensemble ####

cali_ensemble2 = function (config_file, num = NULL, param_file = NULL, cmethod = "LHC", 
                          qualfun = qual_fun, parallel = FALSE, job_name, model = c("FLake", 
                                                                                    "GLM", "GOTM", "Simstrat", "MyLake"), folder = ".", spin_up = NULL, 
                          out_f = "cali", ncores = NULL, tmp_dir = NULL) 
{
  if (!missing(job_name)) {
    if (make.names(job_name) != job_name) {
      stop("job_name '", job_name, "' is not a syntactically valid variable name.")
    }
    call <- match.call()
    call$config_file <- config_file
    call$num <- num
    call$param_file <- param_file
    call$cmethod <- cmethod
    call$parallel <- parallel
    call$model <- model
    call$folder <- folder
    call$spin_up <- spin_up
    call$out_f <- out_f
    nout_fun <- length(qualfun(c(1, 1, 1, 1), c(1.1, 0.9, 
                                                1, 1.2)))
    call_list <- lapply(call, eval)
    call[names(call_list)[-1]] <- call_list[-1]
    script <- make_script(call = call, name = job_name)
    if (!requireNamespace("rstudioapi", quietly = TRUE)) {
      stop("Jobs are only supported in RStudio.")
    }
    if (!rstudioapi::isAvailable("1.2")) {
      stop("Need at least version 1.2 of RStudio to use jobs. Currently running ", 
           rstudioapi::versionInfo()$version, ".")
    }
    job <- rstudioapi::jobRunScript(path = script, name = job_name, 
                                    exportEnv = "R_GlobalEnv")
    return(invisible(job))
  }
  if (!cmethod %in% c("modFit", "LHC", "LHC_old", "MCMC")) {
    stop(paste0("Method ", cmethod, " not allowed. Use one of: modFit, LHC,\n                LHC_old, or MCMC"))
  }
  model <- check_models(model, check_package_install = TRUE)
  check_master_config(config_file, model)
  original_tz <- Sys.getenv("TZ")
  Sys.setenv(TZ = "UTC")
  tz <- "UTC"
  nout_fun <- length(qualfun(c(1, 1, 1, 1), c(1.1, 0.9, 1, 
                                              1.2)))
  oldwd <- getwd()
  on.exit({
    setwd(oldwd)
    Sys.setenv(TZ = original_tz)
    if (file.exists(file.path(folder, "LER_CNFG_TMP.yaml"))) {
      file.remove(file.path(folder, "LER_CNFG_TMP.yaml"))
    }
  })
  yaml <- file.path(folder, config_file)
  start <- gotmtools::get_yaml_value(yaml, label = "time", 
                                     key = "start")
  stop <- gotmtools::get_yaml_value(yaml, label = "location", 
                                    key = "stop")
  obs_file <- gotmtools::get_yaml_value(file = yaml, label = "temperature", 
                                        key = "file")
  time_unit <- gotmtools::get_yaml_value(yaml, "output", "time_unit")
  time_step <- gotmtools::get_yaml_value(yaml, "output", "time_step")
  cnfg_l <- lapply(model, function(m) gotmtools::get_yaml_value(yaml, 
                                                                "config_files", m))
  names(cnfg_l) <- model
  met_timestep <- get_meteo_time_step(file.path(folder, gotmtools::get_yaml_value(yaml, 
                                                                                  "meteo", "file")))
  if (is.null(spin_up)) {
    out_time <- seq.POSIXt(as.POSIXct(start, tz = tz), as.POSIXct(stop, 
                                                                  tz = tz), by = paste(time_step, time_unit))
  }
  else {
    start <- as.POSIXct(start, tz = tz) + spin_up * 24 * 
      60 * 60
    stop <- as.POSIXct(stop, tz = tz)
    out_time <- seq.POSIXt(as.POSIXct(start, tz = tz), as.POSIXct(stop, 
                                                                  tz = tz), by = paste(time_step, time_unit))
  }
  out_time <- data.frame(datetime = out_time)
  if (met_timestep == 86400) {
    out_hour <- lubridate::hour(start)
  }
  else {
    out_hour <- 0
  }
  message("Loading observed wtemp data...")
  obs <- read.csv(file.path(folder, obs_file), stringsAsFactors = FALSE)
  obs$datetime <- as.POSIXct(obs$datetime, tz = tz)
  obs <- obs[obs$datetime %in% out_time$datetime, ]
  obs_deps <- sort(unique(obs$Depth_meter))
  if (any(duplicated(paste0(obs$datetime, obs$Depth_meter)))) {
    warning(paste0("There are non-unique observations in the observed", 
                   " water temperature file ", obs_file, "! Non-unique ", 
                   "observations are averaged."))
  }
  obs_out <- reshape2::dcast(obs, datetime ~ Depth_meter, value.var = "Water_Temperature_celsius", 
                             fun.aggregate = mean, na.rm = TRUE)
  str_depths <- colnames(obs_out)[2:ncol(obs_out)]
  colnames(obs_out) <- c("datetime", paste("wtr_", str_depths, 
                                           sep = ""))
  obs_out$datetime <- as.POSIXct(obs_out$datetime)
  message("Finished!")
  dir.create(file.path(folder, out_f), showWarnings = FALSE)
  configr_master_config <- configr::read.config(yaml)
  cal_section <- configr_master_config[["calibration"]][["met"]]
  params_met <- sapply(names(cal_section), function(n) cal_section[[n]]$initial)
  p_lower_met <- sapply(names(cal_section), function(n) cal_section[[n]]$lower)
  p_upper_met <- sapply(names(cal_section), function(n) cal_section[[n]]$upper)
  p_log_met <- sapply(names(cal_section), function(n) cal_section[[n]]$log)
  cal_section <- configr_master_config[["calibration"]][["Kw"]]
  params_kw <- c(Kw = cal_section$initial)
  p_lower_kw <- c(Kw = cal_section$lower)
  p_upper_kw <- c(Kw = cal_section$upper)
  p_log_kw <- c(Kw = cal_section$log)
  model_p <- model[model %in% names(configr_master_config[["calibration"]])]
  cal_section <- lapply(model_p, function(m) configr_master_config[["calibration"]][[m]])
  names(cal_section) <- model_p
  params_mod <- lapply(model_p, function(m) {
    sapply(names(cal_section[[m]]), function(n) as.numeric(cal_section[[m]][[n]]$initial))
  })
  names(params_mod) <- model_p
  p_lower_mod <- lapply(model_p, function(m) {
    sapply(names(cal_section[[m]]), function(n) as.numeric(cal_section[[m]][[n]]$lower))
  })
  names(p_lower_mod) <- model_p
  p_upper_mod <- lapply(model_p, function(m) {
    sapply(names(cal_section[[m]]), function(n) as.numeric(cal_section[[m]][[n]]$upper))
  })
  names(p_upper_mod) <- model_p
  log_mod <- lapply(model_p, function(m) {
    sapply(names(cal_section[[m]]), function(n) as.logical(cal_section[[m]][[n]]$log))
  })
  names(log_mod) <- model_p
  pars_l <- lapply(model, function(m) {
    df <- data.frame(pars = c(params_met, params_kw, params_mod[[m]], 
                              recursive = TRUE), name = c(names(params_met), names(params_kw), 
                                                          names(params_mod[[m]]), recursive = TRUE), upper = c(p_upper_met, 
                                                                                                               p_upper_kw, p_upper_mod[[m]], recursive = TRUE), 
                     lower = c(p_lower_met, p_lower_kw, p_lower_mod[[m]], 
                               recursive = TRUE), type = c(rep("met", length(params_met)), 
                                                           rep("kw", length(params_kw)), rep("model", length(params_mod[[m]])), 
                                                           recursive = TRUE), log = c(p_log_met, p_log_kw, 
                                                                                      log_mod[[m]], recursive = TRUE), stringsAsFactors = FALSE)
    colnames(df) <- c("pars", "name", "upper", "lower", "type", 
                      "log")
    return(df)
  })
  names(pars_l) <- model
  par_sets <- setNames(sapply(model, function(m) length(pars_l[[m]]$pars)), 
                       model)
  outf_n <- paste0(cmethod, "_", format(Sys.time(), "%Y%m%d%H%M"))
  if (cmethod %in% c("LHC_old", "LHC")) {
    if (!is.null(param_file)) {
      if (length(unique(par_sets)) > 1) {
        stop(paste0("The calibration configuration in the master config file ", 
                    config_file, "results in ", length(unique(par_sets)), 
                    " In this case providing own calibration file is not supported."))
      }
      outf_n <- gsub("_params_", "", basename(param_file))
    }
    if (is.null(param_file)) {
      pars_lhc <- list()
      for (m in model) {
        prange <- matrix(c(pars_l[[m]]$lower, pars_l[[m]]$upper), 
                         ncol = 2)
        prange[pars_l[[m]]$log, ] <- log10(prange[pars_l[[m]]$log, 
        ])
        pars_lhc[[m]] <- FME::Latinhyper(parRange = prange, 
                                         num = num)
        pars_lhc[[m]][, pars_l[[m]]$log] <- 10^pars_lhc[[m]][, 
                                                             pars_l[[m]]$log]
        pars_lhc[[m]] <- signif(pars_lhc[[m]], 5)
        colnames(pars_lhc[[m]]) <- pars_l[[m]]$name
        pars_lhc[[m]] <- as.data.frame(pars_lhc[[m]])
        pars_lhc[[m]]$par_id <- paste0("p", formatC(seq_len(num), 
                                                    width = round(log10(num)) + 1, format = "d", 
                                                    flag = "0"))
        write.table(pars_lhc[[m]], file = file.path(folder, 
                                                    out_f, paste0("params_", m, "_", outf_n, ".csv")), 
                    quote = FALSE, row.names = FALSE, sep = ",")
      }
    }
    else {
      pars_lhc <- lapply(model, function(m) read.csv(param_file, 
                                                     stringsAsFactors = FALSE))
      names(pars_lhc) <- model
      if ((ncol(pars_lhc[[1]]) - 1) != unique(par_sets)) {
        stop(paste0("Number of parameters in file ", 
                    param_file, " (", (ncol(pars_lhc[[1]]) - 1), 
                    ") ", "and number of parameters to calibrate in master config file (", 
                    unique(par_sets), ") do not match!"))
      }
      num <- nrow(pars_lhc[[1]])
    }
  }
  else {
    pars_lhc <- NULL
  }
  if (file.exists(file.path(folder, "LER_CNFG_TMP.yaml"))) {
    warning(strwrap("The file 'LER_CNFG_TMP.yaml' exists in your folder which\n                    is a reserved file name. This will be overwritten."))
    unlink(file.path(folder, "LER_CNFG_TMP.yaml"))
  }
  file.copy(yaml, file.path(folder, "LER_CNFG_TMP.yaml"))
  lst_config_tmp <- configr::read.config(file.path(folder, 
                                                   "LER_CNFG_TMP.yaml"))
  scfctrs_to_calibrate <- names(lst_config_tmp[["calibration"]][["met"]])
  names_scale_section <- names(lst_config_tmp[["scaling_factors"]])
  for (i in names_scale_section) {
    if (i == "all") {
      for (j in scfctrs_to_calibrate) {
        lst_config_tmp[["scaling_factors"]][["all"]][[j]] <- 1
      }
    }
    else {
      for (j in scfctrs_to_calibrate) {
        if (!is.null(lst_config_tmp[["scaling_factors"]][[i]][[j]])) {
          lst_config_tmp[["scaling_factors"]][[i]][[j]] <- 1
        }
      }
    }
  }
  configr::write.config(lst_config_tmp, file.path = file.path(folder, 
                                                              "LER_CNFG_TMP.yaml"), write.type = "yaml", indent = 3)
  export_meteo(config_file = "LER_CNFG_TMP.yaml", model = model, 
               folder = folder)
  met_l <- lapply(model, function(m) {
    met_name <- get_model_met_name(m, cnfg_l[[m]])
    l_names <- as.list(met_var_dic$standard_name)
    names(l_names) <- met_var_dic$short_name
    if (m == "MyLake") {
      met_m <- read.table(file.path(folder, m, met_name), 
                          sep = "\t", header = FALSE)
      colnames(met_m) <- c(l_names$time, l_names$swr, l_names$cc, 
                           l_names$airt, l_names$relh, l_names$p_surf, l_names$wind_speed, 
                           l_names$precip)
    }
    else if (m == "GLM") {
      met_m <- read.table(file.path(folder, m, met_name), 
                          sep = ",", header = TRUE)
    }
    else if (m == "FLake") {
      met_m <- read.table(file.path(folder, m, met_name), 
                          sep = "\t", header = FALSE)
      colnames(met_m) <- c("!Shortwave_Radiation_Downwelling_wattPerMeterSquared", 
                           "Air_Temperature_celsius", "Vapour_Pressure_milliBar", 
                           "Ten_Meter_Elevation_Wind_Speed_meterPerSecond", 
                           "Cloud_Cover_decimalFraction", "datetime")
    }
    else if (m == "GOTM") {
      met_m <- read.table(file.path(folder, m, met_name), 
                          sep = "\t", header = TRUE)
      colnames(met_m)[1] <- "!datetime"
    }
    else if (m == "Simstrat") {
      met_m <- read.table(file.path(folder, m, met_name), 
                          sep = "\t", header = TRUE)
    }
    return(met_m)
  })
  names(met_l) <- model
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl = cl, unclass(lsf.str(envir = asNamespace("LakeEnsemblR"), 
                                                     all = T)), envir = as.environment(asNamespace("LakeEnsemblR")))
    if (cmethod == "LHC_old") {
      parallel::clusterExport(cl, varlist = list("pars_lhc", 
                                                 "pars_l", "model", "config_file", "met_l", "folder", 
                                                 "out_f", "cnfg_l", "obs_deps", "obs_out", "out_hour", 
                                                 "qualfun", "outf_n"), envir = environment())
      message("\nStarted parallel LHC [", Sys.time(), "]\n")
      model_out <- setNames(parLapply(cl, model, function(m) LHC_model(pars = pars_lhc[[m]], 
                                                                       type = pars_l[[m]]$type, model = m, var = "temp", 
                                                                       config_file = config_file, met = met_l[[m]], 
                                                                       folder = folder, out_f = out_f, config_f = cnfg_l[[m]], 
                                                                       obs_deps = obs_deps, obs_out = obs_out, out_hour = out_hour, 
                                                                       qualfun = qualfun, nout_fun = nout_fun, outf_n = outf_n)), 
                            model)
      message("\nFinished parallel LHC [", Sys.time(), 
              "]\n")
    }
    if (cmethod == "LHC") {
      browser()
      model_out <- setNames(lapply(model, function(m) {
        temp_dirs <- make_temp_dir(model = m, folder = folder, 
                                   n = ncores, tmp_dir = tmp_dir)
        param_list <- split(pars_lhc[[m]], rep(1:ncores))
        type <- pars_l[[m]]$type
        met <- met_l[[m]]
        config_f <- cnfg_l[[m]]
        varlist = list("config_file", "m", "temp_dirs", 
                       "type", "met", "obs_out", "out_hour", "config_f", 
                       "nout_fun", "qualfun", "folder", "out_f", "obs_deps", 
                       "outf_n")
        parallel::clusterExport(cl, varlist = varlist, 
                                envir = environment())
        message(m, ": Starting LHC calibration with ", 
                num, " parameters using ", ncores, " cores. [", 
                Sys.time(), "]")
        model_out <- parallel::parLapply(cl, seq_along(param_list), 
                                         function(pars, i) {
                                           temp_dir <- temp_dirs[i]
                                           names_out_qfun <- colnames(qualfun(c(1, 1), 
                                                                              c(0.9, 0.8)))
                                           out_i <- as.data.frame(matrix(NA, nrow = nrow(pars[[i]]), 
                                                                         ncol = nout_fun + 1))
                                           names(out_i) <- c(names_out_qfun, "par_id")
                                           out_i$par_id <- pars[[i]]$par_id
                                           for (p in seq_len(nrow(pars[[i]]))) {
                                             change_pars(config_file = config_file, 
                                                         model = m, pars = pars[[i]][p, -ncol(pars[[i]]), 
                                                                                     drop = FALSE], type = type, met = met, 
                                                         folder = temp_dir)
                                             qual_i <- cost_model(config_file = config_file, 
                                                                  model = m, var = "temp", folder = temp_dir, 
                                                                  obs_deps = obs_deps, obs_out = obs_out, 
                                                                  out_hour = out_hour, qualfun = qualfun, 
                                                                  config_f = config_f)
                                             if (any(is.na(qual_i))) {
                                               qual_i <- setNames(rep(NA, nout_fun), 
                                                                  names_out_qfun)
                                             }
                                             out_i[p, -ncol(out_i)] <- qual_i
                                           }
                                           return(out_i)
                                         }, pars = param_list)
        message(m, ": Finished LHC calibration. [", Sys.time(), 
                "]")
        g1 <- do.call(rbind, model_out)
        g1 <- g1[order(g1$par_id), ]
        out_name <- paste0(m, "_", outf_n, ".csv")
        flsw <- file.exists(file.path(oldwd, out_f, out_name))
        write.table(x = g1, file = file.path(oldwd, out_f, 
                                             out_name), append = ifelse(flsw, TRUE, FALSE), 
                    sep = ",", row.names = FALSE, col.names = ifelse(flsw, 
                                                                     FALSE, TRUE), quote = FALSE)
        return(g1)
      }), model)
    }
    on.exit({
      if (!is.null(tmp_dir)) {
        unlink(tmp_dir, recursive = TRUE, force = TRUE)
      }
    })
    if (cmethod == "MCMC") {
      parallel::clusterExport(cl, varlist = list("pars_lhc", 
                                                 "pars_l", "model", "config_file", "met_l", "folder", 
                                                 "out_f", "cnfg_l", "obs_deps", "obs_out", "out_hour", 
                                                 "qualfun", "outf_n"), envir = environment())
      message("\nStarted parallel MCMC\n")
      model_out <- setNames(parLapply(cl, model, function(m) {
        FME::modMCMC(f = wrap_model, p = setNames(pars_l[[m]]$pars, 
                                                  pars_l[[m]]$name), type = pars_l[[m]]$type, 
                     model = m, var = "temp", config_file = config_file, 
                     met = met_l[[m]], folder = folder, config_f = cnfg_l[[m]], 
                     out_f = out_f, obs_deps = obs_deps, obs_out = obs_out, 
                     out_hour = out_hour, qualfun = function(O, 
                                                             P) {
                       ssr = sum((as.matrix(O[, -1]) - as.matrix(P[, 
                                                                   -1]))^2, na.rm = TRUE)
                     }, outf_n = outf_n, niter = num, lower = setNames(pars_l[[m]]$lower, 
                                                                       pars_l[[m]]$name), upper = setNames(pars_l[[m]]$upper, 
                                                                                                           pars_l[[m]]$name), ...)
      }), model)
      message("\nFinished parallel MCMC\n")
    }
    if (cmethod == "modFit") {
      parallel::clusterExport(cl, varlist = list("pars_lhc", 
                                                 "pars_l", "model", "config_file", "met_l", "folder", 
                                                 "out_f", "cnfg_l", "obs_deps", "obs_out", "out_hour", 
                                                 "qualfun", "outf_n"), envir = environment())
      message("\nStarted parallel modFit\n")
      model_out <- setNames(parLapply(cl, model, function(m) {
        FME::modFit(f = wrap_model, p = setNames(pars_l[[m]]$pars, 
                                                 pars_l[[m]]$name), type = pars_l[[m]]$type, 
                    model = m, var = "temp", config_file = config_file, 
                    met = met_l[[m]], folder = folder, config_f = cnfg_l[[m]], 
                    out_f = out_f, obs_deps = obs_deps, obs_out = obs_out, 
                    out_hour = out_hour, qualfun = function(O, 
                                                            P) {
                      res = na.exclude(as.vector(as.matrix(O[, 
                                                             -1]) - as.matrix(P[, -1])))
                    }, outf_n = "", write = FALSE, lower = setNames(pars_l[[m]]$lower, 
                                                                    pars_l[[m]]$name), upper = setNames(pars_l[[m]]$upper, 
                                                                                                        pars_l[[m]]$name), ...)
      }), model)
      message("\nFinished parallel modFit\n")
    }
  }
  else {
    if (cmethod == "LHC") {
      browser()
      model_out <- setNames(lapply(model, function(m) LHC_model(pars = pars_lhc[[m]], 
                                                                type = pars_l[[m]]$type, model = m, var = "temp", 
                                                                config_file = config_file, met = met_l[[m]], 
                                                                folder = folder, out_f = out_f, config_f = cnfg_l[[m]], 
                                                                obs_deps = obs_deps, obs_out = obs_out, out_hour = out_hour, 
                                                                qualfun = qualfun, nout_fun = nout_fun, outf_n = outf_n)), 
                            model)
    }
    if (cmethod == "MCMC") {
      model_out <- setNames(lapply(model, function(m) {
        message(paste0("\nStarted MCMC for model ", m, 
                       "\n"))
        res <- FME::modMCMC(f = wrap_model, p = setNames(pars_l[[m]]$pars, 
                                                         pars_l[[m]]$name), type = pars_l[[m]]$type, 
                            model = m, var = "temp", config_file = config_file, 
                            met = met_l[[m]], folder = folder, config_f = cnfg_l[[m]], 
                            out_f = out_f, obs_deps = obs_deps, obs_out = obs_out, 
                            out_hour = out_hour, qualfun = function(O, 
                                                                    P) {
                              ssr = sum((as.matrix(O[, -1]) - as.matrix(P[, 
                                                                          -1]))^2, na.rm = TRUE)
                            }, outf_n = outf_n, niter = num, lower = setNames(pars_l[[m]]$lower, 
                                                                              pars_l[[m]]$name), upper = setNames(pars_l[[m]]$upper, 
                                                                                                                  pars_l[[m]]$name), ...)
        message(paste0("\nFinished MCMC for model ", 
                       m, "\n"))
        return(res)
      }), model)
    }
    if (cmethod == "modFit") {
      browser()
      model_out <- setNames(lapply(model, function(m) {
        message(paste0("\nStarted fitting of model ", 
                       m, "\n"))
        res <- FME::modFit(f = wrap_model, p = setNames(pars_l[[m]]$pars, 
                                                        pars_l[[m]]$name), type = pars_l[[m]]$type, 
                           model = m, var = "temp", config_file = config_file, 
                           met = met_l[[m]], folder = folder, config_f = cnfg_l[[m]], 
                           out_f = out_f, obs_deps = obs_deps, obs_out = obs_out, 
                           out_hour = out_hour, qualfun = function(O, 
                                                                   P) {
                             res = na.exclude(as.vector(as.matrix(O[, 
                                                                    -1]) - as.matrix(P[, -1])))
                           }, outf_n = "", write = FALSE, lower = setNames(pars_l[[m]]$lower, 
                                                                           pars_l[[m]]$name), upper = setNames(pars_l[[m]]$upper, 
                                                                                                               pars_l[[m]]$name), ...)
        message(paste0("\nFinished fitting of model ", 
                       m, "\n"))
        return(res)
      }), model)
    }
  }
  return(model_out)
}

# cali_res_FLake_test <- cali_ensemble2(config_file = config_file, num = 1, cmethod = "LHC",
#                                  parallel = T, model = "FLake", ncores = 1)


# export_location ####

export_location2 = function (config_file, model = c("GOTM", "GLM", "Simstrat", "FLake", 
                                 "MyLake"), folder = ".") 
{
  oldwd <- getwd()
  setwd(folder)
  on.exit({
    setwd(oldwd)
  })
  model <- check_models(model)
  lat <- get_yaml_value(config_file, "location", "latitude")
  lon <- get_yaml_value(config_file, "location", "longitude")
  elev <- get_yaml_value(config_file, "location", "elevation")
  max_depth <- get_yaml_value(config_file, "location", "depth")
  init_depth <- get_yaml_value(config_file, "location", "init_depth")
  hyp_file <- get_yaml_value(config_file, "location", "hypsograph")
  if (!file.exists(hyp_file)) {
    stop(hyp_file, " does not exist. Check filepath in ", 
         config_file)
  }
  hyp <- read.csv(hyp_file)
  use_ice <- get_yaml_value(config_file, "ice", "use")
  output_depths <- get_yaml_value(config_file, "output", "depths")
  if ("FLake" %in% model) {
    browser()
    fla_fil <- file.path(folder, get_yaml_value(config_file, 
                                                "config_files", "FLake"))
    bth_area <- hyp$Area_meterSquared
    bth_depth <- hyp$Depth_meter
    top <- min(bth_depth)
    bottom <- max(bth_depth)
    layer_d <- seq(top, bottom, 0.1)
    layer_a <- stats::approx(bth_depth, bth_area, layer_d)$y
    if (init_depth < max_depth) {
      layer_a <- layer_a[(init_depth + (layer_d - max_depth)) >= 
                           0]
      layer_d <- layer_d[(init_depth + (layer_d - max_depth)) >= 
                           0]
      layer_d <- layer_d - min(layer_d)
    }
    vols <- c()
    for (i in 2:length(layer_d)) {
      h <- layer_d[i] - layer_d[i - 1]
      cal_v <- (h/3) * (layer_a[i] + layer_a[i - 1] + sqrt(layer_a[i] * 
                                                             layer_a[i - 1]))
      vols <- c(vols, cal_v)
    }
    vol <- sum(vols)
    # mean_depth <- signif((vol/layer_a[1]), 4)
    mean_depth <- 0.301
    input_nml(fla_fil, label = "SIMULATION_PARAMS", key = "h_ML_in", 
              mean_depth)
    input_nml(fla_fil, label = "LAKE_PARAMS", key = "depth_w_lk", 
              mean_depth)
    input_nml(fla_fil, label = "LAKE_PARAMS", key = "latitude_lk", 
              lat)
  }
  if ("GLM" %in% model) {
    browser()
    glm_nml <- file.path(folder, get_yaml_value(config_file, 
                                                "config_files", "GLM"))
    nml <- read_nml(glm_nml)
    glm_hyp <- hyp
    glm_hyp[, 1] <- elev - glm_hyp[, 1]
    Ao <- max(glm_hyp[, 2])
    bsn_wid <- sqrt((2 * Ao)/pi)
    bsn_len <- 2 * bsn_wid
    min_layer_thick <- get_nml_value(nml, "min_layer_thick")
    max_layers <- round(max_depth/min_layer_thick) + nrow(hyp)
    if (!("crest_elev" %in% names(nml[["morphometry"]]))) {
      nml[["morphometry"]][["crest_elev"]] <- 0
    }
    inp_list <- list(lake_name = get_yaml_value(config_file, 
                                                "location", "name"), latitude = lat, longitude = lon, 
                     lake_depth = max_depth, crest_elev = max((glm_hyp[, 
                                                                       1])), bsn_vals = length(glm_hyp[, 1]), H = rev(glm_hyp[, 
                                                                                                                              1]), A = rev(glm_hyp[, 2]), bsn_len = bsn_len, 
                     bsn_wid = bsn_wid, max_layers = max_layers, max_layer_thick = 1, 
                     lake_depth = init_depth)
    nml <- glmtools::set_nml(nml, arg_list = inp_list)
    write_nml(nml, glm_nml)
  }
  if ("GOTM" %in% model) {
    browser()
    got_yaml <- file.path(folder, get_yaml_value(config_file, 
                                                 "config_files", "GOTM"))
    input_yaml(got_yaml, "location", "name", get_yaml_value(config_file, 
                                                            "location", "name"))
    input_yaml(got_yaml, "location", "latitude", lat)
    input_yaml(got_yaml, "location", "longitude", lon)
    input_yaml(got_yaml, "location", "depth", max_depth)
    input_yaml(got_yaml, "grid", "nlev", round(max_depth/0.5))
    if (use_ice) {
      input_yaml(got_yaml, "ice", "model", 2)
    }
    else {
      input_yaml(got_yaml, "ice", "model", 0)
    }
    ndeps <- nrow(hyp)
    got_hyp <- hyp
    got_hyp[, 1] <- -got_hyp[, 1]
    if (init_depth < max_depth) {
      got_hyp$Depth_meter <- got_hyp$Depth_meter + (max_depth - 
                                                      init_depth)
    }
    colnames(got_hyp) <- c(as.character(ndeps), "2")
    write.table(got_hyp, "GOTM/hypsograph.dat", quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = TRUE)
    input_yaml(got_yaml, "location", "hypsograph", "hypsograph.dat")
  }
  if ("Simstrat" %in% model) {
    browser()
    sim_par <- file.path(folder, get_yaml_value(config_file, 
                                                "config_files", "Simstrat"))
    sim_hyp <- hyp
    sim_hyp[, 1] <- -sim_hyp[, 1]
    if (init_depth < max_depth) {
      sim_hyp$Depth_meter <- sim_hyp$Depth_meter + (max_depth - 
                                                      init_depth)
    }
    colnames(sim_hyp) <- c("Depth [m]", "Area [m^2]")
    write.table(sim_hyp, "Simstrat/hypsograph.dat", quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = TRUE)
    input_json(sim_par, "Input", "Grid", round(max_depth/output_depths))
    input_json(sim_par, "Input", "Morphology", "\"hypsograph.dat\"")
    input_json(sim_par, "ModelParameters", "lat", lat)
    if (use_ice) {
      input_json(sim_par, "ModelConfig", "IceModel", 1)
    }
    else {
      input_json(sim_par, "ModelConfig", "IceModel", 0)
      input_json(sim_par, "ModelConfig", "SnowModel", 0)
    }
    surf_area <- max(sim_hyp[, 2])/1e+06
    a_seiche <- 10^(-2.8591 + 0.7029 * log10(surf_area))
    input_json(sim_par, "ModelParameters", "a_seiche", a_seiche)
  }
  if ("MyLake" %in% model) {
    browser()
    load(get_yaml_value(config_file, "config_files", "MyLake"))
    c_shelter <- 1 - exp(-0.3 * (hyp$Area_meterSquared[1] * 
                                   1e-06))
    mylake_config[["Phys.par"]][5] <- c_shelter
    mylake_config[["Phys.par"]][6] <- lat
    mylake_config[["Phys.par"]][7] <- lon
    if (init_depth < max_depth) {
      myl_hyp <- hyp[(hyp$Depth_meter - max_depth) >= -init_depth, 
      ]
      myl_hyp$Depth_meter <- myl_hyp$Depth_meter - min(myl_hyp$Depth_meter)
      mylake_config[["In.Az"]] <- as.matrix(myl_hyp$Area_meterSquared)
      mylake_config[["In.Z"]] <- as.matrix(myl_hyp$Depth_meter)
    }
    else {
      myl_hyp <- hyp
      mylake_config[["In.Az"]] <- as.matrix(hyp$Area_meterSquared)
      mylake_config[["In.Z"]] <- as.matrix(hyp$Depth_meter)
    }
    mylake_config[["In.FIM"]] <- matrix(rep(0.92, nrow(myl_hyp)), 
                                        ncol = 1)
    mylake_config[["In.Chlz.sed"]] <- matrix(rep(196747, 
                                                 nrow(myl_hyp)), ncol = 1)
    mylake_config[["In.TPz.sed"]] <- matrix(rep(756732, nrow(myl_hyp)), 
                                            ncol = 1)
    mylake_config[["In.DOCz"]] <- matrix(rep(3000, nrow(myl_hyp)), 
                                         ncol = 1)
    mylake_config[["In.Chlz"]] <- matrix(rep(7, nrow(myl_hyp)), 
                                         ncol = 1)
    mylake_config[["In.DOPz"]] <- matrix(rep(7, nrow(myl_hyp)), 
                                         ncol = 1)
    mylake_config[["In.TPz"]] <- matrix(rep(21, nrow(myl_hyp)), 
                                        ncol = 1)
    mylake_config[["In.Sz"]] <- matrix(rep(0, nrow(myl_hyp)), 
                                       ncol = 1)
    mylake_config[["In.Cz"]] <- matrix(rep(0, nrow(myl_hyp)), 
                                       ncol = 1)
    temp_fil <- gsub(".*/", "", get_yaml_value(config_file, 
                                               "config_files", "MyLake"))
    save(mylake_config, file = file.path(folder, "MyLake", 
                                         temp_fil))
  }
  message("export_location complete!")
}

# export_location2(config_file=config_file, model = c("FLake","GLM","GOTM","Simstrat","MyLake"))


# plot_LHC ####

plot_LHC2 = function (config_file, model, res_files, qual_met = "rmse", best_quant = 0.1, 
          best = "low") 
{
  model <- check_models(model)
  if (!best %in% c("low", "high")) {
    stop("best must be either low or high")
  }
  if (length(res_files) != 2 * length(model)) {
    stop(paste0("The number of models (", length(model),
                ") and the number of res_files (", length(res_files),
                ") does not fit. There should be 2 files (results ",
                "and parameter sets) per model."))
  }
  configr_master_config <- configr::read.config(file.path(config_file))
  met_pars <- names(configr_master_config[["calibration"]][["met"]])
  if ("Kw" %in% names(configr_master_config[["calibration"]])) {
    kw_pars <- "Kw"
  }
  else {
    kw_pars <- NULL
  }
  model_p <- model[model %in% names(configr_master_config[["calibration"]])]
  model_pars <- lapply(model_p, function(m) {
    gsub("/", ".", names(configr_master_config[["calibration"]][[m]]))
  })
  names(model_pars) <- model_p
  res <- lapply(res_files, function(f) na.exclude(read.csv(f)))
  names(res) <- basename(gsub("_LHC_.*", "", res_files))
  res <- lapply(model, function(m) merge(res[[m]], res[[paste0("params_", 
                                                               m)]]))
  names(res) <- model
  if (best == "low") {
    best_l <- lapply(model, function(m) {
      subset(res[[m]], res[[m]][[qual_met]] < quantile(res[[m]][[qual_met]], 
                                                       best_quant))
    })
    names(best_l) <- model
  }
  else {
    best_l <- lapply(model, function(m) {
      subset(res[[m]], res[[m]][[qual_met]] > quantile(res[[m]][[qual_met]], 
                                                       (1 - best_quant)))
    })
    names(best_l) <- model
  }
  ret_l <- list()
  for (m in model) {
    if (!(qual_met %in% colnames(res[[m]]))) {
      av_met <- colnames(res[[m]])[!(colnames(res[[m]])) %in% 
                                     c("par_id", met_pars, kw_pars, model_pars[[m]])]
      stop(paste0("Model performance metric ", qual_met, 
                  " not available for model ", m, " available metrics: ", 
                  paste0(av_met, collapse = ", ")))
    }
    for (p in c(met_pars, kw_pars, model_pars[[m]])) {
      ret_l[[m]][[p]] <- ggplot(res[[m]]) + geom_point(aes_string(x = p, 
                                                                  y = qual_met)) + scale_color_gradient(low = "green", 
                                                                                                        high = "red") + ggtitle(paste0("Scatterplot of ", 
                                                                                                                                       p, " for model ", m))
      if (best == "low") {
        best_par <- res[[m]][[p]][res[[m]][[qual_met]] == 
                                    min(res[[m]][[qual_met]])]
        best_q <- res[[m]][[qual_met]][res[[m]][[qual_met]] == 
                                         min(res[[m]][[qual_met]])]
      }
      else {
        best_par <- res[[m]][[p]][res[[m]][[qual_met]] == 
                                    max(res[[m]][[qual_met]])]
        best_q <- res[[m]][[qual_met]][res[[m]][[qual_met]] == 
                                         max(res[[m]][[qual_met]])]
      }
      browser()
      annotations <- data.frame(xpos = -Inf, ypos = Inf, 
                                annotateText = paste0("Best: ", p, " = ", signif(best_par, 
                                                                                 4), "; ", qual_met, " = ", signif(best_q, 4)), 
                                hjustvar = 0, vjustvar = 1)
      ret_l[[m]][[paste0("dist_", p)]] <- ggplot(best_l[[m]]) + 
        geom_histogram(aes_string(x = p, y = "..density.."), 
                       color = "black", fill = "white") + geom_density(aes_string(x = p), 
                                                                       alpha = 0.2, fill = "#FF6666") + xlim(c(min(res[[m]][[p]]), 
                                                                                                               max(res[[m]][[p]]))) + geom_vline(aes_string(xintercept = best_par), 
                                                                                                                                                 color = "blue", linetype = "dashed", size = 1, 
                                                                                                                                                 show.legend = TRUE) + geom_label(data = annotations, 
                                                                                                                                                                                  aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, 
                                                                                                                                                                                      label = annotateText), color = "blue") + ggtitle(paste0("Distribution of best ", 
                                                                                                                                                                                                                                              round(best_quant * 100, 1), "% of ", p, " for model ", 
                                                                                                                                                                                                                                              m))
    }
  }
  return(ret_l)
}

# plot_LHC2(config_file = config_file, model = "GOTM", res_files = cali_res_GOTM$GOTM,
#          qual_met = "nse", best = "high")
 
# export_config ####

export_config2 <- function (config_file, model = c("GOTM", "GLM", "Simstrat", "FLake", 
                                                  "MyLake"), dirs = TRUE, time = TRUE, location = TRUE, output_settings = TRUE, 
                           meteo = TRUE, init_cond = TRUE, extinction = TRUE, flow = TRUE, 
                           model_parameters = TRUE, folder = ".") 
{
  oldwd <- getwd()
  setwd(folder)
  original_tz <- Sys.getenv("TZ")
on.exit({
    setwd(oldwd)
    Sys.setenv(TZ = original_tz)
  })
  Sys.setenv(TZ = "GMT")
  if (!file.exists(config_file)) {
    stop(config_file, " does not exist.")
  }
  check_master_config(config_file, exp_cnf = TRUE)
  model <- check_models(model)
  if (dirs) {
    export_dirs(config_file = config_file, model = model, 
                folder = folder)
  }
  if (time) {
    export_time(config_file = config_file, model = model, 
                folder = folder)
  }
  if (location) {
    export_location(config_file = config_file, model = model, 
                    folder = folder)
  }
  if (output_settings) {
    export_output_settings(config_file = config_file, model = model, 
                           folder = folder)
  }
  if (meteo) {
    export_meteo(config_file = config_file, model = model, 
                 folder = folder)
  }
  if (init_cond) {
    export_init_cond(config_file = config_file, model = model, 
                     print = TRUE, folder = folder)
  }
  if (extinction) {
    export_extinction(config_file = config_file, model = model, 
                      folder = folder)
  }
  if (flow) {
    export_flow(config_file = config_file, model = model, 
                folder = folder)
  }
  if (model_parameters) {
    export_model_parameters(config_file = config_file, model = model, 
                            folder = folder)
  }
}

# plot_heatmap ####

plot_heatmap2 = function (ncdf = NULL, var = "temp", dim = "model", dim_index = 1, 
          var_list = NULL, model = NULL, tile_width = NULL, tile_height = NULL) 
{
  model <- check_models(model)
  if (!is.null(ncdf)) {
    if (!file.exists(ncdf)) {
      stop("File '", ncdf, "' does not exist. Check you have the correct filepath.")
    }
    vars <- gotmtools::list_vars(ncdf)
    if (!(var %in% vars)) {
      stop("Variable '", var, "' is not present in the netCDF file '", 
           ncdf, "'")
    }
    var_list <- load_var(ncdf, var = var, return = "list", 
                         dim = dim, dim_index = dim_index)
  }
  else {
    var_list <- var_list
  }
  if (!is.null(model)) {
    var_list <- var_list[c(model, "Obs")]
  }
  mod_names <- names(var_list)
  data <- var_list %>% reshape2::melt(id.vars = "datetime") %>% 
    dplyr::group_by(datetime)
  colnames(data) <- c("datetime", "Depth", "value", "Model")
  data$depth <- -as.numeric(gsub("wtr_", "", data$Depth))
  data <- as.data.frame(data)
  data$Model <- factor(data$Model)
  data$Model <- factor(data$Model, levels = mod_names)
  if (is.null(tile_width)) {
    tile_width <- as.numeric(difftime(data$datetime[2], data$datetime[1], 
                                      units = "secs"))
  }
  if (is.null(tile_height)) {
    the_depths <- unique(data$depth)
    tile_height <- abs(min(diff(the_depths)))
  }
  spec <- RColorBrewer::brewer.pal(11, "Spectral")
  data <- data[!is.na(data$value), ]
  if (nrow(data) == 0) {
    stop("Modelled  and observed data is all NAs.\n         Please inspect the model output and re-run 'run_ensemble()' if necessary.")
  }
  p1 <- ggplot(data) + geom_tile(aes(datetime, depth, fill = value), 
                                 width = tile_width, height = tile_height) + scale_fill_gradientn(colours = rev(spec)) + 
    facet_wrap(~Model, ncol = 2)
  return(p1)
}

# plot_heatmap2(ncdf)


# load_var ####

load_var2 = function (ncdf, var, return = "list", dim = "model", dim_index = 1, 
                      print = TRUE) 
{
  match.arg(return, c("list", "array"))
  match.arg(dim, c("model", "member"))
  if (!file.exists(ncdf)) {
    stop(ncdf, " does not exist. Check the filepath is correct.")
  }
  if (!var %in% lake_var_dic$short_name) {
    stop(paste0("Variabel '", var, "' unknown. Allowed names for var: ", 
                paste0(lake_var_dic$short_name, collapse = ", ")))
  }
  tryCatch({
    fid <- ncdf4::nc_open(ncdf)
    tim <- ncdf4::ncvar_get(fid, "time")
    tunits <- ncdf4::ncatt_get(fid, "time")
    tustr <- strsplit(tunits$units, " ")
    tdstr <- strsplit(unlist(tustr)[3], "-")
    tmonth <- as.integer(unlist(tdstr)[2])
    tday <- as.integer(unlist(tdstr)[3])
    tyear <- as.integer(unlist(tdstr)[1])
    tdstr <- strsplit(unlist(tustr)[4], ":")
    thour <- as.integer(unlist(tdstr)[1])
    tmin <- as.integer(unlist(tdstr)[2])
    origin <- as.POSIXct(paste0(tyear, "-", tmonth, "-", 
                                tday, " ", thour, ":", tmin), format = "%Y-%m-%d %H:%M", 
                         tz = "UTC")
    time <- as.POSIXct(tim, origin = origin, tz = "UTC")
    mod_names <- ncdf4::ncatt_get(fid, "model", "Model")$value
    mod_names <- strsplit(mod_names, ", ")[[1]]
    mod_names <- substring(mod_names, 5)
    mem <- ncdf4::ncvar_get(fid, "member")
    var1 <- ncdf4::ncvar_get(fid, var)
    tunits <- ncdf4::ncatt_get(fid, var)
    miss_val <- tunits$missing_value
    var1[var1 >= miss_val] <- NA
    var_dim <- strsplit(tunits$coordinates, " ")[[1]]
    if (length(dim(var1)) > 2) {
      z <- ncvar_get(fid, "z")
    }
  }, warning = function(w) {
    return_val <- "Warning"
  }, error = function(e) {
    return_val <- "Error"
    warning("Error creating netCDF file!")
  }, finally = {
    ncdf4::nc_close(fid)
  })
  mat <- matrix(data = c(var, tunits$units), 
                dimnames = list(c("short_name", "units"), 
                                c()))
  if (print == TRUE) {
    message("Extracted ", var, " from ", ncdf)
    print(mat)
  }
  if (length(dim(var1)) == 4) {
    n_vals <- dim(var1)[2] * (dim(var1)[3]) * (dim(var1)[4])
    for (m in seq_len(dim(var1)[1])) {
      nas <- sum(is.na(var1[m, , , ]))
      if (nas == n_vals) {
        break
      }
    }
    if (m != dim(var1)[1]) {
      var1 <- var1[(seq_len(m - 1)), , , ]
    }
  }
  else if (length(dim(var1)) == 3) {
    n_vals <- dim(var1)[2] * (dim(var1)[3])
    for (m in seq_len(dim(var1)[1])) {
      nas <- sum(is.na(var1[m, , ]))
      if (nas == n_vals) {
        break
      }
    }
    if (m != dim(var1)[1]) {
      var1 <- var1[(seq_len(m - 1)), , ]
    }
  }
  if (return == "array") {
    if (length(dim(var1)) == 4) {
      dimnames(var1) <- list(paste0("member_", seq_len(dim(var1)[1])), 
                             mod_names, as.character(time), z)
    }
    if (length(dim(var1)) == 3 & var == "temp") {
      dimnames(var1) <- list(mod_names, as.character(time), 
                             z)
    }
    if (length(dim(var1)) == 3 & var == "ice_height") {
      dimnames(var1) <- list(paste0("member_", seq_len(dim(var1)[1])), 
                             mod_names, as.character(time))
    }
    if (length(dim(var1)) == 2) {
      dimnames(var1) <- list(mod_names, as.character(time))
    }
    return(var1)
  }
  if (return == "list") {
    if ("z" %in% var_dim) {
      if (length(dim(var1)) == 4) {
        if (dim == "model") {
          if (dim_index > dim(var1)[1]) {
            stop("Dimension index ", dim_index, " out of bounds!\nAvailable dimensions: ", 
                 paste(seq_len(dim(var1)[1]), collapse = ","))
          }
          var_list <- lapply(seq(dim(var1)[2]), function(x) var1[dim_index, 
                                                                 x, , ])
          names(var_list) <- mod_names
        }
        else if (dim == "member") {
          if (dim_index > dim(var1)[2]) {
            stop("Dimension index ", dim_index, " out of bounds!\nAvailable dimensions: ", 
                 paste(seq_len(dim(var1)[2]), collapse = ","))
          }
          var_list <- lapply(seq(dim(var1)[1]), function(x) var1[x, 
                                                                 dim_index, , ])
          n_vals <- dim(var_list[[1]])[1] * (dim(var_list[[1]])[2])
          for (m in seq_len(length(var_list))) {
            nas <- sum(is.na(var_list[[m]]))
            if (nas == n_vals) {
              break
            }
          }
          if (m != length(var_list)) {
            var_list <- var_list[(seq_len(m - 1))]
          }
          names(var_list) <- paste0(mod_names[dim_index], 
                                    "_member_", seq_len(length(var_list)))
        }
        var_list <- lapply(var_list, function(x) {
          x <- as.data.frame(x)
          x <- cbind(time, x)
          colnames(x) <- c("datetime", paste0("wtr_", 
                                              abs(z)))
          return(x)
        })
      }
      else if (length(dim(var1)) == 3) {
        var_list <- lapply(seq(dim(var1)[1]), function(x) var1[x, 
                                                               , ])
        names(var_list) <- mod_names
        var_list <- lapply(var_list, function(x) {
          x <- as.data.frame(x)
          x <- cbind(time, x)
          colnames(x) <- c("datetime", paste0("wtr_", 
                                              abs(z)))
          return(x)
        })
      }
    }
    else {
      if (length(dim(var1)) == 2) {
        var_list <- lapply(seq(dim(var1)[1]), function(x) var1[x, 
        ])
        names(var_list) <- mod_names
        var_list <- lapply(var_list, function(x) {
          x <- as.data.frame(x)
          x <- cbind(time, x)
          colnames(x)[2] <- var
          return(x)
        })
      }
      else if (length(dim(var1)) == 3) {
        if (dim == "model") {
          if (dim_index > dim(var1)[2]) {
            stop("Dimension index ", dim_index, " out of bounds!\nAvailable dimensions: ", 
                 paste(seq_len(dim(var1)[2]), collapse = ","))
          }
          var_list <- lapply(seq(dim(var1)[2]), function(x) var1[dim_index, 
                                                                 x, ])
          names(var_list) <- mod_names
        }
        else {
          var_list <- lapply(seq(dim(var1)[1]), function(x) var1[x, 
                                                                 dim_index, ])
          n_vals <- length(var_list[[1]])
          for (m in seq_len(length(var_list))) {
            nas <- sum(is.na(var_list[[m]]))
            if (nas == n_vals) {
              break
            }
          }
          if (m != length(var_list)) {
            var_list <- var_list[(seq_len(m - 1))]
          }
          names(var_list) <- paste0(mod_names[dim_index], 
                                    "_member_", seq_len(length(var_list)))
        }
        var_list <- lapply(var_list, function(x) {
          x <- as.data.frame(x)
          x <- cbind(time, x)
          colnames(x)[2] <- var
          return(x)
        })
      }
    }
  }
  return(var_list)
}

load_var2(ncdf = "output/ensemble_output.nc", var = "temp", return = "list", 
          dim = "model", dim_index = 1)
