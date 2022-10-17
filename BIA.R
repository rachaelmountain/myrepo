## BIA draft code

detach("package:epicR", unload = TRUE)
install.packages("C:/Users/rcm.stu/OneDrive - UBC/Documents/epicR",
                 repos = NULL,
                 type = "source")


# Set-up ------------------------------------------------------------------

{
  # * load required packages ------------------------------------------------
  library(epicR)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)

  # * set seed --------------------------------------------------------------
  #set.seed(444)

  # * settings --------------------------------------------------------------
  record_mode <- c(record_mode_none=0,
                   record_mode_agent=1,
                   record_mode_event=2,
                   record_mode_some_event=3)

  settings <- get_default_settings()
  settings$record_mode <- record_mode["record_mode_none"] # stores the events matrix, switch to "record_mode_none" for faster code
  settings$n_base_agents <- 18000000 #18e+6 # change number of agents - 18mil=est Canadian pop >= 40, 1st Jan 2015 -- https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000501
  settings$n_base_agents <- settings$n_base_agents - 1 # since it starts at id=0 ; ask Amin about this
  settings$update_continuous_outcomes_mode <- 1


  # * common input values ---------------------------------------------------
  # these input values all remain constant across case detection scenarios
  time_horizon <- 13
  yrs_btw_CD <- 20 # ensures no individual tested twice
  CD_start_year <- 8 #2022
  discount_rate <- 0 # the default
  med_adherence <- 0.7 # the default, adjust for sensitivity analysis
  smoking_adherence <- 0.7 # the default, adjust for sensitivity analysis
  p_case_detection <- c(0.1,0.2,0.3,0.4,0.5)
  p_correct_overdiagnosis <- 0.5
  price_yr <- 1.03761




  # Functions ---------------------------------------------------------------


  # * Simulation ------------------------------------------------------------

  # CD_method : "None" "CDQ17" "CDQ195" "CDQ165" "FlowMeter" "FlowMeter_CDQ"
  # eligibility_criteria : "All patients" "Symptomatic" "Eversmokers"

  BIA_simul <- function(CD_method, eligibility_criteria){

    #CD_method <- "FlowMeter" ; eligibility_criteria <- "All patients"

    if(!CD_method %in% c("None", "CDQ17", "CDQ195", "CDQ165", "FlowMeter", "FlowMeter_CDQ")){
      print("Unknown case detection method")
    }

    ### run simulation
    init_session(settings=settings)
    input <- get_input()

    if(eligibility_criteria=="All patients"){
      min_age <- 40
      min_pack_years <- 0
      num_symptoms <- 0
      CD_values <- input$values$diagnosis$case_detection_methods[,CD_method]
    }else if (eligibility_criteria=="Symptomatic"){
      min_age <- 40
      min_pack_years <- 0
      num_symptoms <- 1
      CD_values <- input$values$diagnosis$case_detection_methods_symptomatic[,CD_method]
    }else if (eligibility_criteria=="Eversmokers"){
      min_age <- 50
      min_pack_years <- 0.001
      num_symptoms <- 0
      CD_values <- input$values$diagnosis$case_detection_methods_eversmokers[,CD_method]
    }else{
      terminate_session()
      return("Unknown eligibility criteria")
    }

    input$values$global_parameters$time_horizon <- time_horizon #time_horizon + (CD_start_year-2015) + 1
    input$values$diagnosis$case_detection_start_yr <- CD_start_year
    input$values$diagnosis$years_btw_case_detection <- yrs_btw_CD
    input$values$global_parameters$discount_cost <- discount_rate
    input$values$medication$medication_adherence <- med_adherence
    input$values$smoking$smoking_cessation_adherence <- smoking_adherence
    if(CD_method!="None"){
      #input$values$diagnosis$p_case_detection[(CD_start_year-2015+1):(CD_start_year-2015+time_horizon)] <- p_case_detection
      input$values$diagnosis$p_case_detection[CD_start_year:(CD_start_year+4)] <- p_case_detection
    }
    input$values$diagnosis$p_correct_overdiagnosis <- p_correct_overdiagnosis

    input$values$diagnosis$min_cd_age <- min_age
    input$values$diagnosis$min_cd_pack_years <- min_pack_years
    input$values$diagnosis$min_cd_symptoms <- num_symptoms

    input$values$cost$cost_case_detection <- CD_values[3]*price_yr
    input$values$diagnosis$logit_p_prevalent_diagnosis_by_sex["case_detection.None",] <- CD_values[1]
    input$values$diagnosis$logit_p_diagnosis_by_sex["case_detection.None",] <- CD_values[1]
    input$values$diagnosis$logit_p_overdiagnosis_by_sex["case_detection.None",] <- CD_values[2]

    input$values$cost$bg_cost_by_stage <- input$values$cost$bg_cost_by_stage*price_yr
    input$values$cost$cost_smoking_cessation <- input$values$cost$cost_smoking_cessation*price_yr
    input$values$cost$exac_dcost <- input$values$cost$exac_dcost*price_yr
    input$values$cost$cost_outpatient_diagnosis <- input$values$cost$cost_outpatient_diagnosis*price_yr
    input$values$cost$cost_gp_visit <- input$values$cost$cost_gp_visit*price_yr
    input$values$medication$medication_costs <- input$values$medication$medication_costs*price_yr

    set.seed(123)
    run(input = input$values)

    #inputs <- Cget_inputs()
    output <- Cget_output()
    output_ex <- Cget_output_ex()
    # eventmat <- as.data.frame(Cget_all_events_matrix())
    # eventmat <- eventmat %>%
    #   mutate(ctime=local_time+time_at_creation,.after=local_time)


    terminate_session()


    # number alive
    base <- data.frame(year=1:(time_horizon),
                       criteria=eligibility_criteria,
                       CD_method=CD_method,
                       alive=rowSums(output_ex$n_alive_by_ctime_sex),
                       sex=output_ex$n_alive_by_ctime_sex[,2],
                       smoke=output_ex$n_smoking_status_by_ctime[,2]
    )

    # number with COPD and diagnosed
    diags <- data.frame(copd=rowSums(output_ex$n_COPD_by_ctime_severity[,2:5]),
                        copd_sev=output_ex$n_COPD_by_ctime_severity[,4],
                        diags_true=rowSums(output_ex$n_Diagnosed_by_ctime_severity[,2:5]), # diagnosed here is true diagnosed only
                        diags_false=rowSums(output_ex$n_Overdiagnosed_by_ctime_sex),
                        diags_total=rowSums(output_ex$n_Diagnosed_by_ctime_severity[,2:5])+rowSums(output_ex$n_Overdiagnosed_by_ctime_sex)
    )

    # health outcomes
    healthOC <- data.frame(exacs=rowSums(output_ex$n_exac_by_ctime_severity),
                           exacs_mild=output_ex$n_exac_by_ctime_severity[,1],
                           exacs_mod=output_ex$n_exac_by_ctime_severity[,2],
                           exacs_sev=output_ex$n_exac_by_ctime_severity[,3],
                           exacs_vsev=output_ex$n_exac_by_ctime_severity[,4],
                           deaths=rowSums(output_ex$n_exac_death_by_ctime_severity))

    # number of case detections
    case_detection <- data.frame(case_detection=rowSums(output_ex$n_case_detection_by_ctime),
                                 CD_pos=output_ex$n_case_detection_by_ctime[,1],
                                 CD_neg=output_ex$n_case_detection_by_ctime[,2])

    # costs
    costs <- output_ex$annual_cost_ctime

    hosps <- output_ex$n_exac_by_ctime_severity[,3:4]
    hosps_cost <- input$values$cost$exac_dcost[3:4]
    hosps_cost_yr <- rowSums(t(c(hosps_cost)*t(hosps)))

    copd <- output_ex$n_COPD_by_ctime_severity[,2:5] # column 1 is gold grade 0 i.e. no copd
    main_cost <- input$values$cost$bg_cost_by_stage[2:5]
    copd_cost_yr <- rowSums(t(c(main_cost)*t(copd)))

    caseds <- output_ex$n_case_detection_by_ctime
    caseds_cost <- c(input$values$cost$cost_case_detection,
                     input$values$cost$cost_case_detection+input$values$cost$cost_outpatient_diagnosis)
    caseds_cost_yr <- rowSums(t(c(caseds_cost)*t(caseds)))

    medyrs <- output_ex$medication_time_by_ctime_class[,1:4]
    meds_cost <- input$values$medication$medication_costs[c(2,5,7,15)]
    meds_cost_yr <- rowSums(t(c(meds_cost)*t(medyrs)))

    costs_subcs <- data.frame(cost_total=costs,
                              cost_case_detection=caseds_cost_yr,
                              cost_treat=meds_cost_yr,
                              cost_hosp=hosps_cost_yr,
                              cost_maint=copd_cost_yr)

    results <- cbind(base,diags,healthOC,case_detection,costs_subcs)


    return(results)

  }

}





# Run simulations ---------------------------------------------------------




# * S1: All patients ------------------------------------------------------

s10 <- BIA_simul(CD_method = "None", eligibility_criteria = "All patients")
saveRDS(s10,"s10_2.Rda")
s1a <- BIA_simul(CD_method = "CDQ17", eligibility_criteria = "All patients")
saveRDS(s1a,"s1a_2.Rda")
s1b <- BIA_simul(CD_method = "FlowMeter", eligibility_criteria = "All patients")
saveRDS(s1b,"s1b_2.Rda")
s1c <- BIA_simul(CD_method = "FlowMeter_CDQ", eligibility_criteria = "All patients")
saveRDS(s1c,"s1c_2.Rda")


# * S2: Symptomatic patients ----------------------------------------------

#s20 <- BIA_simul(CD_method = "None", eligibility_criteria = "Symptomatic")
s2a <- BIA_simul(CD_method = "FlowMeter", eligibility_criteria = "Symptomatic")
saveRDS(s2a,"s2a_2.Rda")


# * S3: Eversmokers + age >= 50 -------------------------------------------

#s30 <- BIA_simul(CD_method = "None", eligibility_criteria = "Eversmokers")
s3a <- BIA_simul(CD_method = "CDQ195", eligibility_criteria = "Eversmokers")
saveRDS(s3a,"s3a_2.Rda")
s3b <- BIA_simul(CD_method = "CDQ165", eligibility_criteria = "Eversmokers")
saveRDS(s3b,"s3b_2.Rda")
s3c <- BIA_simul(CD_method = "FlowMeter", eligibility_criteria = "Eversmokers")
saveRDS(s3c,"s3c_2.Rda")
s3d <- BIA_simul(CD_method = "FlowMeter_CDQ", eligibility_criteria = "Eversmokers")
saveRDS(s3d,"s3d_2.Rda")





s10 <- readRDS("s10_2.Rda")
s1a <- readRDS("s1a_2.Rda")
s1b <- readRDS("s1b_2.Rda")
s1c <- readRDS("s1c_2.Rda")
s2a <- readRDS("s2a_2.Rda")
s3a <- readRDS("s3a_2.Rda")
s3b <- readRDS("s3b_2.Rda")
s3c <- readRDS("s3c_2.Rda")
s3d <- readRDS("s3d_2.Rda")






# Results -----------------------------------------------------------------

baseline_yr <- CD_start_year-1
CD_yrs <- CD_start_year:(CD_start_year+4)

s10 <- s10 %>% mutate(scenario="s10", .after=year)
s1a <- s1a %>% mutate(scenario="s1a", .after=year)
s1b <- s1b %>% mutate(scenario="s1b", .after=year)
s1c <- s1c %>% mutate(scenario="s1c", .after=year)

#s20 <- s20 %>% mutate(scenario="s20", .after=year)
s2a <- s2a %>% mutate(scenario="s2a", .after=year)

#s30 <- s30 %>% mutate(scenario="s30", .after=year)
s3a <- s3a %>% mutate(scenario="s3a", .after=year)
s3b <- s3b %>% mutate(scenario="s3b", .after=year)
s3c <- s3c %>% mutate(scenario="s3c", .after=year)
s3d <- s3d %>% mutate(scenario="s3d", .after=year)


all_results <- bind_rows(s10,s1a,s1b,s1c,s2a,s3a,s3b,s3c,s3d) %>% #
  filter(year %in% c(baseline_yr,CD_yrs)) %>% # remove unnecessary data
  mutate(year=year-baseline_yr) %>% # renumber so baseline year = 0 then 1:5 case detection years
  mutate(copd_prev=copd/alive*100, .after=copd) %>%
  mutate(diags_true_prev=diags_true/alive*100, .after=diags_true) %>%
  mutate(diags_false_prev=diags_false/(alive-copd)*100, .after=diags_false) %>%
  mutate(diags_total_prev=diags_total/alive*100, .after=diags_total) %>%
  mutate(cost_total_pp=cost_total/copd) %>%
  mutate(cost_case_detection_pp=cost_case_detection/copd) %>%
  mutate(cost_treat_pp=cost_treat/copd) %>%
  mutate(cost_hosp_pp=cost_hosp/copd) %>%
  mutate(cost_maint_pp=cost_maint/copd)



baseline_compar_raw <- all_results %>%
  filter(year==0)

baseline_compar_percent <- baseline_compar_raw %>%
  mutate(across(alive:cost_maint_pp, function(x) (x - mean(x))/mean(x)*100))


all_results_costs <- all_results %>%
  select(year,scenario,criteria,CD_method,cost_total:cost_maint) %>%
  pivot_longer(cols=cost_total:cost_maint,names_to="cost_group",values_to="CAD") %>%
  pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)



bia_table <- function(table,base_scenario,alt_scenario){

  base_table <- table %>%
    filter(scenario==base_scenario) %>%
    select(year,scenario,criteria,CD_method,cost_total:cost_maint)

  alt_table <- table %>%
    filter(scenario==alt_scenario) %>%
    select(year,scenario,criteria,CD_method,cost_total:cost_maint)

  bia_table <- data.frame(year=0:5,
                          scenario=NA,
                          criteria=NA,
                          CD_method=NA,
                          bia_total=base_table$cost_total-alt_table$cost_total,
                          bia_case_detection=base_table$cost_case_detection-alt_table$cost_case_detection,
                          bia_treat=base_table$cost_treat-alt_table$cost_treat,
                          bia_hosp=base_table$cost_hosp-alt_table$cost_hosp,
                          bia_maint=base_table$cost_maint-alt_table$cost_maint)

  table <- bind_rows(
    base_table %>%
      pivot_longer(cols=cost_total:cost_maint,names_to="cost_group",values_to="CAD") %>%
      pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)
    ,
    alt_table %>%
      pivot_longer(cols=cost_total:cost_maint,names_to="cost_group",values_to="CAD") %>%
      pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)
    ,
    bia_table %>%
      pivot_longer(cols=bia_total:bia_maint,names_to="cost_group",values_to="CAD") %>%
      pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)
  )

}

bia_s10_s1a <- bia_table(all_results,"s10","s1a")
bia_s10_s1b <- bia_table(all_results,"s10","s1b")
bia_s10_s1c <- bia_table(all_results,"s10","s1c")




all_results_costs_long <- all_results_costs %>%
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="CAD")
ggplot(all_results_costs_long %>% filter(cost_group=="cost_total"),
       aes(x=year,y=CAD,group=scenario,col=scenario)) +
  geom_point() +
  geom_line() +
  theme_bw()


all_results_costs_long <- all_results_costs %>%
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="CAD")
ggplot(all_results_costs_long %>% filter(cost_group=="cost_case_detection"),
       aes(x=year,y=CAD,group=scenario,col=scenario)) +
  geom_point() +
  geom_line() +
  theme_bw()





all_results_ho <- all_results %>%
  select(year,scenario,criteria,CD_method,copd_sev,diags_true,diags_false,exacs,exacs_mild,exacs_mod,exacs_sev,exacs_vsev,deaths) %>%
  pivot_longer(cols=copd_sev:deaths,names_to="ho_group",values_to="count") %>%
  pivot_wider(names_from=year,names_prefix="year_",values_from=count)

all_results_ho_long <- all_results_ho %>%
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="count")
p1 <- ggplot(all_results_ho_long %>% filter(ho_group=="copd_sev"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p2 <- ggplot(all_results_ho_long %>% filter(ho_group=="diags_true"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p3 <- ggplot(all_results_ho_long %>% filter(ho_group=="diags_false"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p4 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p5 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_mild"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p6 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_mod"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p7 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_sev"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p8 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_vsev"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p9 <- ggplot(all_results_ho_long %>% filter(ho_group=="deaths"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
grid.arrange(p2,p3,nrow=2)
grid.arrange(p5,p6,p7,p8,nrow=2)








# Per COPD patient results ------------------------------------------------

all_results_costs_pp <- all_results %>%
  select(year,scenario,criteria,CD_method,cost_total_pp:cost_maint_pp) %>%
  pivot_longer(cols=cost_total_pp:cost_maint_pp,names_to="cost_group",values_to="CAD") %>%
  pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)



bia_table_pp <- function(table,base_scenario,alt_scenario){

  base_table <- table %>%
    filter(scenario==base_scenario) %>%
    select(year,scenario,criteria,CD_method,cost_total_pp:cost_maint_pp)

  alt_table <- table %>%
    filter(scenario==alt_scenario) %>%
    select(year,scenario,criteria,CD_method,cost_total_pp:cost_maint_pp)

  bia_table <- data.frame(year=0:5,
                          scenario=NA,
                          criteria=NA,
                          CD_method=NA,
                          bia_total_pp=base_table$cost_total_pp-alt_table$cost_total_pp,
                          bia_case_detection_pp=base_table$cost_case_detection_pp-alt_table$cost_case_detection_pp,
                          bia_treat_pp=base_table$cost_treat_pp-alt_table$cost_treat_pp,
                          bia_hosp_pp=base_table$cost_hosp_pp-alt_table$cost_hosp_pp,
                          bia_maint_pp=base_table$cost_maint_pp-alt_table$cost_maint_pp)

  table <- bind_rows(
    base_table %>%
      pivot_longer(cols=cost_total_pp:cost_maint_pp,names_to="cost_group",values_to="CAD") %>%
      pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)
    ,
    alt_table %>%
      pivot_longer(cols=cost_total_pp:cost_maint_pp,names_to="cost_group",values_to="CAD") %>%
      pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)
    ,
    bia_table %>%
      pivot_longer(cols=bia_total_pp:bia_maint_pp,names_to="cost_group",values_to="CAD") %>%
      pivot_wider(names_from=year,names_prefix="year_",values_from=CAD)
  )

  return(table)

}

bia_pp_s10_s1a <- bia_table_pp(all_results,"s10","s1a")
bia_pp_s10_s1b <- bia_table_pp(all_results,"s10","s1b")
bia_pp_s10_s1c <- bia_table_pp(all_results,"s10","s1c")




all_results_costs_long_pp <- all_results_costs_pp %>%
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="CAD") %>%
  mutate(year=factor(year)) %>%
  mutate(year=fct_recode(year,"2021"="year_0","2022"="year_1","2023"="year_2","2024"="year_3","2025"="year_4","2026"="year_5"))
ggplot(all_results_costs_long_pp %>% filter(cost_group=="cost_total_pp"),
       aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Total cost per COPD patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()

all_results_costs_long_pp <- all_results_costs_pp %>%
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="CAD")
p1 <- ggplot(all_results_costs_long_pp %>% filter(cost_group=="cost_case_detection_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p2 <- ggplot(all_results_costs_long_pp %>% filter(cost_group=="cost_treat_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Treatment cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p3 <- ggplot(all_results_costs_long_pp %>% filter(cost_group=="cost_hosp_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Hospitalisation cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
p4 <- ggplot(all_results_costs_long_pp %>% filter(cost_group=="cost_maint_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Maintenance cost (CAD)") +
  geom_point() +
  geom_line() +
  theme_bw()
grid.arrange(p1,p2,p3,p4,nrow=2)
