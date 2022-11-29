library(tidyverse)
library(here)
library(ggplot2)
library(dplyr)


##########################################

# build data set


# icerdata <- data.frame(Scenario=rep(c("NoCD","S1a","S1b","S1c","S2a","S3a","S3b","S3c","S3d"),each=2),
#                        Testing_Interval=rep(c("3 years","5 years"),9),
#                        CostpAgent=c(2151,2151,2438,2356,2363,2296,2386,2313,2286,2246,2234,2207,2292,2250,2256,2224,2263,2227),
#                        QALYpAgent=c(12.546,12.546,12.560,12.556,12.554,12.552,12.551,12.550,12.553,12.551,12.548,12.548,12.553,12.552,12.550,12.549,12.549,12.548)
#                        )
# 
# icerdata$IncrementalCosts <- icerdata$CostpAgent - 2151
# icerdata$IncrementalQALYs <- icerdata$QALYpAgent - 12.546
# icerdata$ICER <- icerdata$IncrementalCosts/icerdata$IncrementalQALYs





##########################################


icerdata <- read_csv("CEPlane.csv")

icerdata$Scenario <- ifelse(icerdata$Scenario=="S1NoCDAvg" | icerdata$Scenario=="S2NoCD" | icerdata$Scenario=="S3NoCD",
                            "NoCD", icerdata$Scenario)

icerdata <- icerdata %>% filter(Testing_Interval=="5 years")

### Cost-effectiveness  plane ###
p <- ggplot(data=icerdata, aes(x=IncrementalQALYs, y=IncrementalCosts)) + 
  geom_point(aes(shape=Testing_Interval)) +
  geom_text(aes(label=Scenario),vjust=-0.5, size=3) + 
  #scale_y_continuous(expand=expand_scale(mult=c(0.07,0), add=c(0,100))) +
  #scale_x_continuous(expand=expand_scale(mult=c(0.07,0)), limits=c(0,0.016)) +
  geom_abline(intercept=0, slope=25000, colour="grey", linetype="dashed") +
  theme_bw() + labs(x="Incremental QALYs", y="Incremental Costs", col="Testing Interval",
                    shape="Testing Interval") +
  theme(legend.position="none")



##############################
### Assessing dominance ####
### First round ###
# Most cost-effective scenario- S3b-5 yrs NEW REF
icerdata %>% 
  arrange(ICER) %>% 
  slice(1)

### Second round ###
# S3b-5 years new REF, Reassess S1a, S1b, S2a-3 yrs, S3b-3yrs
icerdata.red <- icerdata %>% 
  filter(Scenario=="S1a"|Scenario=="S1b"|Scenario=="S3b"|
           (Scenario=="S2a"&Testing_Interval=="3 years")) %>%
  arrange(CostpAgent) 

icerdata.red <- icerdata.red %>% 
  mutate(IncrementalCosts2 = CostpAgent - as.numeric(icerdata.red[1,"CostpAgent"]),
         IncrementalQALYs2 = QALYpAgent - as.numeric(icerdata.red[1,"QALYpAgent"]),
         ICER2= IncrementalCosts2/IncrementalQALYs2)

icerdata.red2 <- filter(icerdata.red, ICER2<25000) # S3b-3 yrs, S1a, S2a-3yrs 

icerdata.red2 %>% arrange(ICER2) %>% 
  slice(1)  %>% 
  select(Scenario, Testing_Interval, IncrementalCosts2, IncrementalQALYs2, ICER2)

# S1a-3 yrs lowest ICER, nothing left!
#### S1a at 3 years is the most cost-effective with an ICER of $21,108/QALY
##############################

### Adding lines to CEPlane ###
ceplane <- data.frame(icerdata %>%  
                        select(Scenario, Testing_Interval, IncrementalCosts, IncrementalQALYs)) 

p <- p + geom_segment(data=ceplane, aes(x=0, xend=ceplane[ceplane$Scenario=="S3b" & ceplane$Testing_Interval=="5 years",4], 
                                        y=0, yend=ceplane[ceplane$Scenario=="S3b" & ceplane$Testing_Interval=="5 years",3]), 
                      size=0.01, col="blue") #+
  geom_segment(data=ceplane, aes(x=ceplane[ceplane$Scenario=="S3b" & ceplane$Testing_Interval=="5 years",4],
                                 xend=ceplane[ceplane$Scenario=="S1a" & ceplane$Testing_Interval=="3 years",4],
                                 y=ceplane[ceplane$Scenario=="S3b" & ceplane$Testing_Interval=="5 years",3],
                                 yend=ceplane[ceplane$Scenario=="S1a" & ceplane$Testing_Interval=="3 years",3]),
               size=0.01, col="blue")

#### Highlight winning scenario  ######
ceplane_dom <- filter(ceplane, Scenario=="S3b", Testing_Interval=="5 years")

p <- p + geom_point(data=ceplane_dom, aes(x=IncrementalQALYs, y=IncrementalCosts), col="red", shape=16) +
  geom_text(data=ceplane_dom, aes(label=Scenario),vjust=-0.5, size=3, col="red")

p + annotate("text", label="WTP $20,000 \n S3b most CE - ICER=$16,251", 
             x=0.0125, y=30, size=3)

####

ggsave("WTP20.pdf",width=10,height=6)
#ggsave("Figures/CEPlaneFigure.tiff", width=6, height=6, device="tiff", dpi=300)  


