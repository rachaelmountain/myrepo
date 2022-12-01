#' Returns results of validation tests for overdiagnosis
#' @param n_sim number of agents
#' @return validation test results
#' @export
validate_overdiagnosis <- function(n_sim = 1e+04) {
  
  n_sim=150000

  settings <- get_default_settings()
  settings$record_mode <- record_mode["record_mode_none"]
  settings$agent_stack_size <- 0
  settings$n_base_agents <- n_sim
  settings$event_stack_size <- 0
  init_session(settings = settings)
  
  model_input <- get_input()
  input <- model_input$values
  input$diagnosis$p_correct_overdiagnosis <- 0.65
  
  res <- run(input = input)
  if (res < 0)
    stop("Execution stopped.\n")
  
  inputs <- Cget_inputs()
  output_ex <- Cget_output_ex()
  
  message("Here are the proportion of non-COPD subjects overdiagnosed over model time: \n")
  
  overdiag <- data.frame(Year=1:inputs$global_parameters$time_horizon,
                         NonCOPD=output_ex$n_COPD_by_ctime_severity[,1],
                         Overdiagnosed=rowSums(output_ex$n_Overdiagnosed_by_ctime_sex))
  
  overdiag$Proportion <- overdiag$Overdiagnosed/overdiag$NonCOPD
  
  print(overdiag)
  
  message("The average proportion overdiagnosed from year", round(length(overdiag$Proportion)/2,0), "to", length(overdiag$Proportion), "is",
          mean(overdiag$Proportion[(round(length(overdiag$Proportion)/2,0)):(length(overdiag$Proportion))]),"\n")
  
  overdiag.plot <- tidyr::gather(data=overdiag, key="Variable", value="Number", c(NonCOPD, Overdiagnosed))
  
  overdiag.plotted <- ggplot2::ggplot(overdiag.plot, aes(x=Year, y=Number, col=Variable)) +
    geom_line() + geom_point() + expand_limits(y = 0) +
    theme_bw() + ylab("Number of non-COPD subjects") + xlab("Years")
  
  plot(overdiag.plotted)
  
  
  
  overdiag.plot <- overdiag.plot %>% 
    mutate(Proportion=ifelse(Variable=="NonCOPD",1-Proportion,Proportion))
  
  
  overdiag.plotted2 <- ggplot2::ggplot(overdiag.plot, aes(x=Year, y=Proportion, col=Variable)) +
    geom_line() + geom_point() + expand_limits(y = 0) +
    theme_bw() + ylab("Proportion of non-COPD subjects") + xlab("Years")
  
  
  
  
  
  
  
  
  
  
  
  message("\n")
  
  terminate_session()
  
}

validate_overdiagnosis(n_sim = 100000)
