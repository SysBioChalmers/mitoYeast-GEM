# Create plots related to protein import, Fe/S cluster biosynthesis and cofactor requirements
#install.packages("extrafont")
# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(extrafont)

setwd('~/Documents/MATLAB/mitoYeast-GEM/ecMitoYeast/Results/')

# load data
filename <- 'proteinImport_w_carriers_comparison.csv'
data <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

# Convert g/g total protein to %
data <- data %>%
  rowwise() %>%
  mutate(model_wo_carriers = g.g.D...0.1.h.1.without.carriers*100,
         model_w_carriers = g.g.D...0.1.h.1.with.carriers*100,
         proteomics = g.g.Proteomics*100)

# Load file with info on protein import complexes
protImportInfo <- read.delim('../ComplementaryData/proteinImport/importMachineryComponents.tsv',
                             sep = '\t',dec = '.',header = TRUE)[,c(1,2,3,5)]
# Merge TIM23 and PAM into one complex for counting abundance
protImportInfo[which(protImportInfo$Component == 'PAM'),4] <- 'TIM23'

# Merge data from model with and import complex info
combinedData <- merge(protImportInfo,data,by.x = 'UniProtID',by.y = 'UniprotID')
combinedData[which(combinedData$Component == 'Disulfide relay'),4] <- 'MIA'

# Plot comparison between model prediction using estimated kcats and curated kcats to proteomis data at D = 0.1 h-1
dataToPlot <- combinedData[,c(4,9,10,11)]
dataToPlot <- aggregate(cbind(dataToPlot$model_wo_carriers,
                              dataToPlot$model_w_carriers,
                              dataToPlot$proteomics), 
                              by = list(Component = dataToPlot$Component), FUN = sum)
colnames(dataToPlot)[2:4] <- c('Model w/o metabolite carriers',
                                     'Model w/ metabolite carriers',
                                     'Proteomics')
# Melt data
dataToPlot <- melt(dataToPlot,id=c('Component'))
# Create a grouped bar plot
fig3_suppl <- ggplot(dataToPlot,aes(fill=variable,y=value,x=Component)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        axis.text.x = element_text(angle = 90,size = 7,face = 'bold',color = 'black', family = 'Arial'),
        axis.ticks = element_line(size = 1,color = 'black'),
        axis.text.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.position = c(0.4,0.8),
        panel.border = element_rect(color = 'black',size = 1,fill = NA),
        panel.background = element_blank()) +
  scale_fill_brewer(palette = 'Blues') +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.5)) +
  ylab('Percent of total protein mass (%)')
fig3_suppl
ggarrange(fig3_suppl,labels = 'B')

# Plot only TIM22 and TOM
dataToPlot <- dataToPlot[c(6,8),]
# melt data
dataToPlot <- melt(dataToPlot,id=c('Component'))
# Create a grouped bar plot
figS4 <- ggplot(dataToPlot,aes(fill=variable,y=value,x=Component)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        axis.text.x = element_text(angle = 90,size = 7,face = 'bold',color = 'black', family = 'Arial'),
        axis.ticks = element_line(size = 1,color = 'black'),
        axis.text.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.position = c(0.25,0.8),
        panel.border = element_rect(color = 'black',size = 1,fill = NA),
        panel.background = element_blank()) +
  scale_fill_brewer(palette = 'Blues') +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25)) +
  ylab('Percent of total protein mass (%)')
figS4
ggsave('~/Documents/Manuscripts/Modeling paper/Submission related/Submission/iScience/Revision/FigS4.pdf',
       plot=figS4,device = 'pdf', dpi = 400, height = 8,width = 8, units = 'cm')

# Plot only model w/ carriers vs proteomics at D = 0.1 h-1
dataToPlot <- combinedData[,c(4,10,11)]
dataToPlot <- aggregate(cbind(dataToPlot$model_w_carriers,
                              dataToPlot$proteomics), 
                        by = list(Component = dataToPlot$Component), FUN = sum)
colnames(dataToPlot)[2:3] <- c('Model','Proteomics')

# Melt data
dataToPlot <- melt(dataToPlot,id=c('Component'))
# Create a grouped bar plot
model_w_carriers_vs_proteomics <- ggplot(dataToPlot,aes(fill=variable,y=value,x=Component)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        axis.text.x = element_text(angle = 90,size = 7,face = 'bold',color = 'black', family = 'Arial'),
        axis.ticks = element_line(size = 1,color = 'black'),
        axis.text.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.position = c(0.4,0.8),
        panel.border = element_rect(color = 'black',size = 1,fill = NA),
        panel.background = element_blank()) +
  scale_fill_brewer(palette = 'Blues') +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.5)) +
  ylab('Percent of total protein mass (%)')
model_w_carriers_vs_proteomics

# Plot over range of dilution rates
# Load data
filename <- 'proteinImport_w_carriers.csv'
data <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

# Convert g/g total protein to %
data <- data %>%
  rowwise() %>%
  mutate(model_D_01 = g.g.D...0.1.h.1*100,
         proteomics_D_01 = g.g.Proteomics*100,
         model_D_02 = g.g.D...0.2.h.1*100,
         model_D_028 = g.g.D...0.28.h.1*100,
         model_D_035 = g.g.D...0.35.h.1*100,
         model_0.39 = g.g.D...0.4.h.1*100,
         model_batch = Batch*100)

# Merge data from model with and import complex info
combinedData <- merge(protImportInfo,data,by.x = 'UniProtID',by.y = 'UniprotID')
combinedData[which(combinedData$Component == 'Disulfide relay'),4] <- 'MIA'

# Extract data for plot for all dilution rates
dataToPlot <- combinedData[,c(4,13:19)]

# Compile data for each complex
dataToPlot <- aggregate(cbind(dataToPlot$model_D_01,dataToPlot$proteomics_D_01,dataToPlot$model_D_02,
                              dataToPlot$model_D_028,dataToPlot$model_D_035,dataToPlot$model_0.39,
                              dataToPlot$model_batch),
                        by = list(Component = dataToPlot$Component), FUN=sum)
colnames(dataToPlot)[2:8] <- c('Model D = 0.1 h-1','Proteomics D = 0.1 h-1','Model D = 0.2 h-1',
                               'Model D = 0.28 h-1','Model D = 0.35 h-1','Model D = 0.4 h-1',
                               'Model Batch')

# Melt data
dataToPlot <- melt(dataToPlot,id=c('Component'))

# Divide the data to plot TIM23 and TOM separate from other components
dataToPlot_TIM23andTOM <- dataToPlot[which(dataToPlot$Component == 'TIM23' | dataToPlot$Component == 'TOM' ),]
dataToPlot_remaining <- dataToPlot[-which(dataToPlot$Component == 'TIM23' | dataToPlot$Component == 'TOM' ),]

# Create a grouped bar plot for TIM23 and TOM
fig3_supplB <- ggplot(dataToPlot_TIM23andTOM,aes(fill=variable,y=value,x=Component)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        axis.ticks.x = element_line(size = 1,color = 'black'),
        axis.ticks.y = element_line(size = 1,color = 'black'),
        axis.text.x = element_text(angle = 0,size = 7,face = 'bold',color = 'black',family = 'Arial'),
        axis.text.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7,face = 'bold'),
        legend.position = c(0.7,0.8),
        panel.border = element_rect(color = 'black',size = 1,fill = NA),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 1,color = 'black')) +
  ylab('Percent of total protein mass (%)') +
  scale_fill_brewer(palette = 'Blues') +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.9))
fig3_supplB

# Create a bar plot for the remaining components
fig3_supplC <- ggplot(dataToPlot_remaining,aes(fill=variable,y=value,x=Component)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        axis.text.x = element_text(angle = 0,size = 7,face = 'bold',color = 'black',family = 'Arial'),
        axis.text.y = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.title = element_blank(),
        legend.text = element_text(size = 7,face = 'bold',color = 'black',family = 'Arial'),
        legend.position = c(0.7,0.8),
        panel.border = element_rect(color = 'black',size = 1,fill = NA),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 1,color = 'black')) +
  ylab('Percent of total protein mass (%)') +
  scale_fill_brewer(palette = 'Blues') +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.15001))
fig3_supplC

# Create a scatter plot comparing model predictions after adding metabolite carriers
dataForScatter <- combinedData[,c(4,6,7)]
colnames(dataForScatter)[2:3] <- c('Model','Proteomics')
dataForScatter <- aggregate(cbind(dataForScatter$Model,dataForScatter$Proteomics),
                            by = list(Component = dataForScatter$Component),FUN = sum)
colnames(dataForScatter)[2:3] <- c('Model','Proteomics')
# log transform values
dataForScatter <- dataForScatter %>%
  rowwise() %>% 
  mutate(Model = log10(Model),
         Proteomics = log10(Proteomics))
# Plot data
fig3_supplD <- ggscatter(dataForScatter, x = 'Proteomics', y = 'Model',
                    add = 'reg.line', conf.int = FALSE,
                    cor.coef = TRUE, cor.method = 'pearson',
                    label = 'Component',repel = TRUE,
                    xlab = 'in vivo log10(g protein/g total protein mass)',
                    ylab = 'in silico log10(g protein/g total protein mass)',
                    font.label = c(8,'black','bold'),label.rectangle = TRUE,
                    cor.coef.coord = c(-4,-2.5)) +
  theme(axis.text = element_text(size = 8,color = 'black',face = 'bold'))
fig3_supplD

# Combine plots
ggarrange(fig3_suppl,fig3_supplD,nrow = 1,labels = c('A','B'))
ggsave('~/Documents/Manuscripts/Modeling paper/Submission related/Submission/iScience/Revision/FigureS3.pdf',
       device = 'pdf',height = 7, width = 14,units = 'cm',dpi = 300)
