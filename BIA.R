## BIA draft code

detach("package:epicR", unload = TRUE)
install.packages("C:/Users/rcm.stu/OneDrive - UBC/Documents/epicR",
                 repos = NULL,
                 type = "source")

remotes::install_github('resplab/epicR')


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
  settings$record_mode <- record_mode["record_mode_agent"] # stores the events matrix, switch to "record_mode_none" for faster code
  settings$n_base_agents <- 1000 #18e+6 # change number of agents - 18mil=est Canadian pop >= 40, 1st Jan 2015 -- https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000501
  settings$n_base_agents <- settings$n_base_agents - 1 # since it starts at id=0 ; ask Amin about this
  settings$update_continuous_outcomes_mode <- 1
  settings$random_number_agent_refill <- 1


  # * common input values ---------------------------------------------------
  # these input values all remain constant across case detection scenarios
  time_horizon <- 14
  yrs_btw_CD <- 20 # ensures no individual tested twice
  CD_start_cal_year <- 2022
  CD_start_year <- CD_start_cal_year - 2015
  discount_rate <- 0 # the default
  med_adherence <- 0.7 # the default, adjust for sensitivity analysis
  smoking_adherence <- 0.7 # the default, adjust for sensitivity analysis
  p_case_detection <- c(0.1,0.2,0.3,0.4,0.5)
  p_correct_overdiagnosis <- 0.65 # did my own validation on this
  price_yr <- 1.03761 # Jan 2022
  
  
  
  settings$runif_buffer_size <- time_horizon*80
  settings$rexp_buffer_size <- time_horizon*150
  settings$rnorm_buffer_size <- time_horizon*20
  settings$rgamma_buffer_size <- time_horizon + 1




# Functions ---------------------------------------------------------------


  # * EPIC simulation ------------------------------------------------------------

  # CD_method : "None" "CDQ17" "CDQ195" "CDQ165" "FlowMeter" "FlowMeter_CDQ"
  # eligibility_criteria : "All patients" "Symptomatic" "Eversmokers"

  BIA_simul <- function(CD_method, eligibility_criteria){

    CD_method <- "CDQ17" ; eligibility_criteria <- "All patients"

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
      input$values$diagnosis$p_case_detection[(CD_start_year+1):(CD_start_year+length(p_case_detection))] <- p_case_detection
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

    #set.seed(444)
    run(input = input$values)
    
    ag <- Cget_agent_events(1)[,121:130]
    for(i in 1:(settings$n_base_agents)){
      ag <- rbind(ag,Cget_agent_events(i+1)[,121:130])
    }
    
    
    max(ag$norm_count.1)
    
    max(ag$unif_count.1)
    
    max(ag$exp_count.1)
    
    max(ag$gamma_COPD_count.1[which(ag$gamma_COPD_count.1<50000)])
    
    max(ag$gamma_NCOPD_count.1[which(ag$gamma_NCOPD_count.1<50000)])
    
    
    
    # > summary(ag$norm_count.1)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 7.00   31.00   31.00   29.27   32.00   34.00 
    # > summary(ag$unif_count.1)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 12.00   71.00   75.00   70.05   77.00  638.00 
    # > summary(ag$exp_count.1)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 3.0   101.0   113.0   110.8   126.0  1329.0 
    
    
    output <- Cget_output()
    output_ex <- Cget_output_ex()

    terminate_session()

    
    
    
    
    
    # number alive
    base <- data.frame(year=1:(time_horizon),
                       criteria=eligibility_criteria,
                       CD_method=CD_method,
                       alive=rowSums(output_ex$n_alive_by_ctime_sex),
                       sex=output_ex$n_alive_by_ctime_sex[,2], # is this male or female???
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
                                 CD_true_pos=output_ex$n_case_detection_by_ctime[,3],
                                 CD_false_pos=output_ex$n_case_detection_by_ctime[,2])

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
                     input$values$cost$cost_case_detection+input$values$cost$cost_outpatient_diagnosis,
                     input$values$cost$cost_case_detection+input$values$cost$cost_outpatient_diagnosis)
    caseds_cost_yr <- rowSums(t(c(caseds_cost)*t(caseds)))

    medyrs <- output_ex$medication_time_by_ctime_class[,1:4]
    medyrs_combos <- data.frame(SABA=medyrs[,1],
                                LAMA=medyrs[,3]-medyrs[,2],
                                LAMA_LABA=medyrs[,2]-medyrs[,4],
                                ICS_LAMA_LABA=medyrs[,4])
    meds_cost <- input$values$medication$medication_costs[c(2,5,7,15)]
    meds_cost_yr <- rowSums(t(c(meds_cost)*t(medyrs_combos)))

    costs_subcs <- data.frame(cost_total=costs,
                              cost_case_detection=caseds_cost_yr,
                              cost_treat=meds_cost_yr,
                              cost_hosp=hosps_cost_yr,
                              cost_maint=copd_cost_yr)

    results <- cbind(base,diags,healthOC,case_detection,costs_subcs)
    
    
    
    # overall results
    overall <- data.frame(n_agents=output_ex$n_agents_post_CD,
                          n_diagnosed=output_ex$n_diagnosed_true_post_CD,
                          n_CD_eligible=output_ex$n_case_detection_eligible)
    


    return(list(overall=overall,results=results))

    }
  }



  # * BIA table - compare two scenarios -------------------------------------
  
  bia_table <- function(table,base_scenario,alt_scenario){
    
    base_table <- table %>%
      filter(scenario==base_scenario) %>%
      select(year,scenario,criteria,CD_method,cost_total:cost_maint)
    
    alt_table <- table %>%
      filter(scenario==alt_scenario) %>%
      select(year,scenario,criteria,CD_method,cost_total:cost_maint)
    
    bia_table <- data.frame(year=0:5,
                            scenario=alt_scenario,
                            criteria=alt_table$criteria[1],
                            CD_method=alt_table$CD_method[1],
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
    
    return(table)
    
  }





  # * BIA per COPD patient table - compare two scenarios --------------------
  
  bia_table_pp <- function(table,base_scenario,alt_scenario){
    
    base_table <- table %>%
      filter(scenario==base_scenario) %>%
      select(year,scenario,criteria,CD_method,cost_total_pp:cost_maint_pp)
    
    alt_table <- table %>%
      filter(scenario==alt_scenario) %>%
      select(year,scenario,criteria,CD_method,cost_total_pp:cost_maint_pp)
    
    bia_table <- data.frame(year=0:5,
                            scenario=alt_scenario,
                            criteria=alt_table$criteria[1],
                            CD_method=alt_table$CD_method[1],
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









# Run simulations ---------------------------------------------------------


  # * S1: All patients ------------------------------------------------------
  
  s10 <- BIA_simul(CD_method = "None", eligibility_criteria = "All patients")
  saveRDS(s10,"s10_3.Rda")
  s1a <- BIA_simul(CD_method = "CDQ17", eligibility_criteria = "All patients")
  saveRDS(s1a,"s1a_3.Rda")
  s1b <- BIA_simul(CD_method = "FlowMeter", eligibility_criteria = "All patients")
  saveRDS(s1b,"s1b_3.Rda")
  s1c <- BIA_simul(CD_method = "FlowMeter_CDQ", eligibility_criteria = "All patients")
  saveRDS(s1c,"s1c_3.Rda")
  
  
  # * S2: Symptomatic patients ----------------------------------------------
  
  #s20 <- BIA_simul(CD_method = "None", eligibility_criteria = "Symptomatic")
  s2a <- BIA_simul(CD_method = "FlowMeter", eligibility_criteria = "Symptomatic")
  saveRDS(s2a,"s2a_3.Rda")
  
  
  # * S3: Eversmokers + age >= 50 -------------------------------------------
  
  #s30 <- BIA_simul(CD_method = "None", eligibility_criteria = "Eversmokers")
  s3a <- BIA_simul(CD_method = "CDQ195", eligibility_criteria = "Eversmokers")
  saveRDS(s3a,"s3a_3.Rda")
  s3b <- BIA_simul(CD_method = "CDQ165", eligibility_criteria = "Eversmokers")
  saveRDS(s3b,"s3b_3.Rda")
  s3c <- BIA_simul(CD_method = "FlowMeter", eligibility_criteria = "Eversmokers")
  saveRDS(s3c,"s3c_3.Rda")
  s3d <- BIA_simul(CD_method = "FlowMeter_CDQ", eligibility_criteria = "Eversmokers")
  saveRDS(s3d,"s3d_3.Rda")






# Read in data ------------------------------------------------------------

s10 <- readRDS("s10_3.Rda")
s1a <- readRDS("s1a_3.Rda")
s1b <- readRDS("s1b_3.Rda")
s1c <- readRDS("s1c_3.Rda")
s2a <- readRDS("s2a_3.Rda")
s3a <- readRDS("s3a_3.Rda")
s3b <- readRDS("s3b_3.Rda")
s3c <- readRDS("s3c_3.Rda")
s3d <- readRDS("s3d_3.Rda")









# Results -----------------------------------------------------------------



# * prepare data ----------------------------------------------------------

# overall output
s10o <- s10$overall %>% mutate(scenario="s10") %>% relocate(scenario)
s1ao <- s1a$overall %>% mutate(scenario="s1a") %>% relocate(scenario)
s1bo <- s1b$overall %>% mutate(scenario="s1b") %>% relocate(scenario)
s1co <- s1c$overall %>% mutate(scenario="s1c") %>% relocate(scenario)
s2ao <- s2a$overall %>% mutate(scenario="s2a") %>% relocate(scenario)
s3ao <- s3a$overall %>% mutate(scenario="s3a") %>% relocate(scenario)
s3bo <- s3b$overall %>% mutate(scenario="s3b") %>% relocate(scenario)
s3co <- s3c$overall %>% mutate(scenario="s3c") %>% relocate(scenario)
s3do <- s3d$overall %>% mutate(scenario="s3d") %>% relocate(scenario)

all_overall <- bind_rows(s10o,s1ao,s1bo,s1co,s2ao,s3ao,s3bo,s3co,s3do) # combine



# set baseline year and case detection years
baseline_yr <- CD_start_year+1-1
CD_yrs <- (CD_start_year+1):(CD_start_year+length(p_case_detection))

# year-by-year results
s10 <- s10$results %>% mutate(scenario="s10", .after=year)
s1a <- s1a$results %>% mutate(scenario="s1a", .after=year)
s1b <- s1b$results %>% mutate(scenario="s1b", .after=year)
s1c <- s1c$results %>% mutate(scenario="s1c", .after=year)

#s20 <- s20 %>% mutate(scenario="s20", .after=year)
s2a <- s2a$results %>% mutate(scenario="s2a", .after=year)

#s30 <- s30 %>% mutate(scenario="s30", .after=year)
s3a <- s3a$results %>% mutate(scenario="s3a", .after=year)
s3b <- s3b$results %>% mutate(scenario="s3b", .after=year)
s3c <- s3c$results %>% mutate(scenario="s3c", .after=year)
s3d <- s3d$results %>% mutate(scenario="s3d", .after=year)


# combine all year-by-year results
all_results <- bind_rows(s10,s1a,s1b,s1c,s2a,s3a,s3b,s3c,s3d) %>% 
  filter(year %in% c(baseline_yr,CD_yrs)) %>% # remove unnecessary data
  mutate(year=year-baseline_yr) %>% # renumber so baseline year = 0 then 1:5 case detection years
  mutate(copd_prev=copd/alive*100, .after=copd) %>%
  mutate(diags_true_prev=diags_true/alive*100, .after=diags_true) %>%
  mutate(diags_false_prev=diags_false/(alive-copd)*100, .after=diags_false) %>%
  mutate(diags_total_prev=diags_total/alive*100, .after=diags_total) %>%
  mutate(cost_total_pp=cost_total/alive) %>%
  mutate(cost_case_detection_pp=cost_case_detection/alive) %>%
  mutate(cost_treat_pp=cost_treat/alive) %>%
  mutate(cost_hosp_pp=cost_hosp/alive) %>%
  mutate(cost_maint_pp=cost_maint/alive)









# * compare baseline years ------------------------------------------------

baseline_compar_raw <- all_results %>%
  filter(year==0) %>% 
  select(-year)

baseline_compar_percent <- baseline_compar_raw %>%
  mutate(across(alive:cost_maint_pp, function(x) (x - mean(x))/mean(x)*100))


baseline_avg <- baseline_compar_raw %>% 
  summarize(alive=mean(alive),
            sex=mean(sex),
            smoke=mean(smoke),
            copd=mean(copd),
            copd_prev=mean(copd_prev),
            copd_sev=mean(copd_sev),
            diags_true=mean(diags_true),
            diags_true_prev=mean(diags_true_prev),
            diags_false=mean(diags_false),
            diags_false_prev=mean(diags_false_prev),
            diags_total=mean(diags_total),
            diags_total_prev=mean(diags_total_prev),
            exacs=mean(exacs),
            exacs_mild=mean(exacs_mild),
            exacs_mod=mean(exacs_mod),
            exacs_sev=mean(exacs_sev),
            exacs_vsev=mean(exacs_vsev),
            deaths=mean(deaths),
            case_detection=0,
            CD_true_pos=0,
            CD_false_pos=0,
            cost_total=mean(cost_total),
            cost_case_detection=0,
            cost_treat=mean(cost_treat),
            cost_hosp=mean(cost_hosp),
            cost_maint=mean(cost_maint),
            cost_total_pp=mean(cost_total_pp),
            cost_case_detection_pp=0,
            cost_treat_pp=mean(cost_treat_pp),
            cost_hosp_pp=mean(cost_hosp_pp),
            cost_maint_pp=mean(cost_maint_pp)
            )


all_results[which(all_results$year==0),5:35] <- baseline_avg
  



# * overall results -------------------------------------------------------

# number eligible / tested / true positives / false positives
overall <- all_overall %>% 
  inner_join(
    all_results %>% 
      filter(year>0) %>% 
      group_by(scenario) %>% 
      summarize(case_detection=sum(case_detection),
                true_positive=sum(CD_true_pos),
                false_positive=sum(CD_false_pos),
                exacs=sum(exacs)),
    by="scenario"
  ) %>% 
  mutate(percent_pop_CD=case_detection/n_agents*100, .after=case_detection) %>% 
  mutate(percent_eligible_CD=case_detection/n_CD_eligible*100, .after=percent_pop_CD) %>%   
  mutate(percent_eligible=n_CD_eligible/n_agents*100, .after=n_CD_eligible) %>% 
  mutate(true_positive_rate=true_positive/case_detection*100, .after=true_positive) %>% 
  mutate(false_positive_rate=false_positive/case_detection*100, .after=false_positive)






# * budget impact ---------------------------------------------------------


# run budget impact function
bia_s10_s1a <- bia_table(all_results,"s10","s1a")
bia_s10_s1b <- bia_table(all_results,"s10","s1b")
bia_s10_s1c <- bia_table(all_results,"s10","s1c")
bia_s10_s2a <- bia_table(all_results,"s10","s2a")
bia_s10_s3a <- bia_table(all_results,"s10","s3a")
bia_s10_s3b <- bia_table(all_results,"s10","s3b")
bia_s10_s3c <- bia_table(all_results,"s10","s3c")
bia_s10_s3d <- bia_table(all_results,"s10","s3d")


bia_all_long <- bind_rows(bia_s10_s1a,
                          bia_s10_s1b,
                          bia_s10_s1c,
                          bia_s10_s2a,
                          bia_s10_s3a,
                          bia_s10_s3b,
                          bia_s10_s3c,
                          bia_s10_s3d) %>% 
  distinct() %>% 
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="CAD") %>% 
  mutate(year=fct_recode(year,"2021"="year_0","2022"="year_1","2023"="year_2","2024"="year_3","2025"="year_4","2026"="year_5"))





# cost plots
p1 <- ggplot(bia_all_long %>% filter(cost_group=="cost_total"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Total cost (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p2 <- ggplot(bia_all_long %>% filter(cost_group=="cost_case_detection"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Case detection cost (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p3 <- ggplot(bia_all_long %>% filter(cost_group=="cost_treat"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Medication cost (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p4 <- ggplot(bia_all_long %>% filter(cost_group=="cost_hosp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Hospitalization cost (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
grid.arrange(p1,p2,p3,p4,nrow=2)





# bia plots
p5 <- ggplot(bia_all_long %>% filter(cost_group=="bia_total" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Total budget impact (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p6 <- ggplot(bia_all_long %>% filter(cost_group=="bia_case_detection" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Case detection budget impact (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p7 <- ggplot(bia_all_long %>% filter(cost_group=="bia_treat" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Medication budget impact (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p8 <- ggplot(bia_all_long %>% filter(cost_group=="bia_hosp" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Hospitalization budget impact (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
grid.arrange(p5,p6,p7,p8,nrow=2)




# total cost and BI over time horizon
bia_total <- bia_all_long %>% 
  filter(year!="2021") %>% 
  group_by(scenario,cost_group) %>% 
  summarize(total_CAD=sum(CAD)) %>% 
  pivot_wider(names_from=scenario,values_from=total_CAD) 


cost_increase_percent <- bia_total %>% 
  filter(cost_group=="bia_total") %>% 
  select(-s10,-cost_group) %>% 
  pivot_longer(cols=s1a:s3d,names_to="scenario",values_to="bia") %>% 
  mutate(cost_increase=-bia/bia_total$s10[which(bia_total$cost_group=="cost_total")]*100) %>% 
  mutate(percent_of_HC_budget=-bia/308000000000*100)
  






# * budget impact per COPD patient ------------------------------------------


# run budget impact function
bia_pp_s10_s1a <- bia_table_pp(all_results,"s10","s1a")
bia_pp_s10_s1b <- bia_table_pp(all_results,"s10","s1b")
bia_pp_s10_s1c <- bia_table_pp(all_results,"s10","s1c")
bia_pp_s10_s2a <- bia_table_pp(all_results,"s10","s2a")
bia_pp_s10_s3a <- bia_table_pp(all_results,"s10","s3a")
bia_pp_s10_s3b <- bia_table_pp(all_results,"s10","s3b")
bia_pp_s10_s3c <- bia_table_pp(all_results,"s10","s3c")
bia_pp_s10_s3d <- bia_table_pp(all_results,"s10","s3d")


bia_pp_all_long <- bind_rows(bia_pp_s10_s1a,
                             bia_pp_s10_s1b,
                             bia_pp_s10_s1c,
                             bia_pp_s10_s2a,
                             bia_pp_s10_s3a,
                             bia_pp_s10_s3b,
                             bia_pp_s10_s3c,
                             bia_pp_s10_s3d) %>% 
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="CAD") %>% 
  mutate(year=fct_recode(year,"2021"="year_0","2022"="year_1","2023"="year_2","2024"="year_3","2025"="year_4","2026"="year_5"))





# cost plots
p1 <- ggplot(bia_pp_all_long %>% filter(cost_group=="cost_total_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Total cost per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p2 <- ggplot(bia_pp_all_long %>% filter(cost_group=="cost_case_detection_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Case detection cost per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p3 <- ggplot(bia_pp_all_long %>% filter(cost_group=="cost_treat_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Medication cost per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p4 <- ggplot(bia_pp_all_long %>% filter(cost_group=="cost_hosp_pp"),
             aes(x=year,y=CAD,group=scenario,col=scenario)) +
  ylab("Hospitalization cost per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
grid.arrange(p1,p2,p3,p4,nrow=2)





# bia plots
p5 <- ggplot(bia_pp_all_long %>% filter(cost_group=="bia_total_pp" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Total budget impact per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p6 <- ggplot(bia_pp_all_long %>% filter(cost_group=="bia_case_detection_pp" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Case detection budget impact per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p7 <- ggplot(bia_pp_all_long %>% filter(cost_group=="bia_treat_pp" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Medication budget impact per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
p8 <- ggplot(bia_pp_all_long %>% filter(cost_group=="bia_hosp_pp" & scenario!="s10"),
             aes(x=year,y=CAD*(-1),group=scenario,col=scenario)) +
  ylab("Hospitalization budget impact per patient (CAD)") + xlab("Year") +
  geom_point() +
  geom_line() +
  theme_bw()
grid.arrange(p5,p6,p7,p8,nrow=2)



# total cost and BI over time horizon
bia_pp_total <- bia_pp_all_long %>% 
  filter(year!="2021") %>% 
  group_by(scenario,cost_group) %>% 
  summarize(total_CAD=sum(CAD)) %>% 
  pivot_wider(names_from=scenario,values_from=total_CAD)






# * health outcomes results -----------------------------------------------


all_results_ho <- all_results %>%
  select(year,scenario,criteria,CD_method,copd_sev,diags_true,diags_false,exacs,exacs_mild,exacs_mod,exacs_sev,exacs_vsev,deaths) %>%
  pivot_longer(cols=copd_sev:deaths,names_to="ho_group",values_to="count") %>%
  pivot_wider(names_from=year,names_prefix="year_",values_from=count)

all_results_ho_long <- all_results_ho %>%
  pivot_longer(cols=year_0:year_5,names_to="year",values_to="count") %>% 
  mutate(year=fct_recode(year,"2021"="year_0","2022"="year_1","2023"="year_2","2024"="year_3","2025"="year_4","2026"="year_5"))

p1 <- ggplot(all_results_ho_long %>% filter(ho_group=="copd_sev"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Severe COPD") +
  geom_point() +
  geom_line() +
  theme_bw()
p2 <- ggplot(all_results_ho_long %>% filter(ho_group=="diags_true"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("True COPD diagnoses") +
  geom_point() +
  geom_line() +
  theme_bw()
p3 <- ggplot(all_results_ho_long %>% filter(ho_group=="diags_false"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Overdiagnosis") +
  geom_point() +
  geom_line() +
  theme_bw()
p4 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Total exacerbations)") +
  geom_point() +
  geom_line() +
  theme_bw()
p5 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_mild"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Mild exacerbations") +
  geom_point() +
  geom_line() +
  theme_bw()
p6 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_mod"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Moderate exacerbations") +
  geom_point() +
  geom_line() +
  theme_bw()
p7 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_sev"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Severe exacerbations") +
  geom_point() +
  geom_line() +
  theme_bw()
p8 <- ggplot(all_results_ho_long %>% filter(ho_group=="exacs_vsev"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("very severe exacerbations") +
  geom_point() +
  geom_line() +
  theme_bw()
p9 <- ggplot(all_results_ho_long %>% filter(ho_group=="deaths"),
             aes(x=year,y=count,group=scenario,col=scenario)) +
  ylab("Deaths") +
  geom_point() +
  geom_line() +
  theme_bw()
grid.arrange(p2,p3,nrow=2)
grid.arrange(p5,p6,p7,p4,nrow=2)














