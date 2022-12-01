validate_diagnosis <- function(n_sim = 1e+04) {

  n_sim = 150000
  
  settings <- get_default_settings()
  settings$record_mode <- record_mode["record_mode_none"]
  settings$agent_stack_size <- 0
  settings$n_base_agents <- n_sim
  settings$event_stack_size <- 0
  init_session(settings = settings)
  
  input <- get_input()
  
  res <- run(input = input$values)

  
  inputs <- Cget_inputs()
  output_ex <- Cget_output_ex()
  
  message("Here are the proportion of COPD patients diagnosed over model time: \n")
  
  diag <- data.frame(Year=1:inputs$global_parameters$time_horizon,
                     COPD=rowSums(output_ex$n_COPD_by_ctime_sex),
                     Diagnosed=rowSums(output_ex$n_Diagnosed_by_ctime_sex))
  
  diag$Proportion <- round(diag$Diagnosed/diag$COPD,2)
  
  print(diag)
  
  message("The average proportion diagnosed from year", round(length(diag$Proportion)/2,0), "to", length(diag$Proportion), "is",
          mean(diag$Proportion[(round(length(diag$Proportion)/2,0)):(length(diag$Proportion))]),"\n")
  
  diag.plot <- tidyr::gather(data=diag, key="Variable", value="Number", c(COPD,Diagnosed))
  
  diag.plotted <- ggplot2::ggplot(diag.plot, aes(x=Year, y=Number, col=Variable)) +
    geom_line() + geom_point() + expand_limits(y = 0) +
    theme_bw() + ylab("Number of COPD patients") + xlab("Years")
  
  plot(diag.plotted)
  
  
  
  
  diag.plot <- diag.plot %>% 
    mutate(Proportion=ifelse(Variable=="COPD",1-Proportion,Proportion)) %>% 
    mutate(Variable=ifelse(Variable=="COPD","Undiagnosed",Variable))
  
  
  diag.plotted2 <- ggplot2::ggplot(diag.plot, aes(x=Year, y=Proportion, col=Variable)) +
    geom_line() + geom_point() + expand_limits(y = 0) +
    theme_bw() + ylab("Proportion of COPD subjects") + xlab("Years")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  message("\n")
  message("Now let's look at the proportion diagnosed by COPD severity.\n")
  
  prop <- data.frame(Year=1:inputs$global_parameters$time_horizon,
                     output_ex$n_Diagnosed_by_ctime_severity/output_ex$n_COPD_by_ctime_severity)[,c(1,3,4,5,6)]
  
  names(prop) <- c("Year","GOLD1","GOLD2","GOLD3","GOLD4")
  prop <- prop[-1,]
  print(prop)
  
  message("The average proportion of GOLD 1 and 2 that are diagnosed from year", round(nrow(prop)/2,0), "to", max(prop$Year), "is",
          (mean(prop$GOLD1[round((nrow(prop)/2),0):nrow(prop)]) + mean(prop$GOLD2[round((nrow(prop)/2),0):nrow(prop)]))/2,"\n")
  
  prop.plot <- tidyr::gather(data=prop, key="GOLD", value="Proportion", c(GOLD1:GOLD4))
  
  prop.plotted <- ggplot2::ggplot(prop.plot, aes(x=Year, y=Proportion, col=GOLD)) +
    geom_line() + geom_point() + expand_limits(y = 0) +
    theme_bw() + ylab("Proportion diagnosed") + xlab("Years")
  
  plot(prop.plotted)
  
  terminate_session()
}






