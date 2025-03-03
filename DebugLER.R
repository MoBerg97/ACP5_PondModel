options(error=browser)

cali_res <- cali_ensemble(config_file = config_file, num = 10, cmethod = "LHC",
                          parallel = TRUE, model = model)
