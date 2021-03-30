# Plot Guangzhou Cases
rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggpubr)
library(stringr)
options(stringsAsFactors = F)
casefile <- "./data/Guangzhou_Cases.csv"
casedata <- read.csv(casefile, header=T, sep=",")
casedatamelt <- melt(casedata, id.vars=c("DateID","Date"), 
                     measure.vars=c("New_cases", "Imported_cases"),
                     variable.name="Cases",value.name="Numbers")

# figure
casemeltfigure <- ggplot(casedatamelt) + 
                geom_bar( aes(DateID, Numbers, fill=Cases), 
                         stat = "identity" , position='identity',
                         alpha=0.5, width=1 ) +
                scale_fill_manual(values = c("black" , "red")) +
                scale_x_discrete(limits = seq(1,length(casedata$DateID),14),
                                 labels = casedata$Date[seq(1,length(casedata$DateID),14) ]) + 
                labs(x= "", y="Daily new cases") +
                theme( axis.line = element_line(colour='black'),
                       panel.background = element_rect(fill = "transparent",colour = "black"),
                       legend.key.size = unit(0.2,'cm'))

cumulativefigure <- ggplot(casedata) +
                    geom_line(aes(x=DateID, y=Cumulatively_Import_case), color="red" ) + 
                    geom_line(aes(x=DateID, y=Cumulatively_Local_Imported_associated_NonS), color="black" ) +
  scale_x_discrete(limits = seq(1,length(casedata$DateID),14),
                   labels = casedata$Date[seq(1,length(casedata$DateID),14) ]) + 
  labs(x= "", y="Daily new cases") +
  theme( axis.line = element_line(colour='black'),
         panel.background = element_rect(fill = "transparent",colour = "black"))

mergeplot <- ggarrange(casemeltfigure, cumulativefigure, nrow=2, align= "h")
ggsave(filename = str_replace( casefile, ".csv", ".pdf"), width=12, height=8, units="cm", scale=2) 
