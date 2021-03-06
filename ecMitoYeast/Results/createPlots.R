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

### Create plots related to protein import ###
# load data
filename <- 'proteinImport.csv'
data <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

# Convert g/g total protein to %
data <- data %>%
  rowwise() %>%
  mutate(model_estimated_kcats = g.g.estimated.kcats*100,
         model_D_01 = g.g.D...0.1.h.1*100,
         proteomics_D_01 = g.g.Proteomics*100,
         model_D_02 = g.g.D...0.2.h.1*100,
         model_D_028 = g.g.D...0.28.h.1*100,
         model_D_035 = g.g.D...0.35.h.1*100,
         model_0.39 = g.g.D...0.4.h.1*100,
         model_batch = Batch*100)

# Load file with info on protein import complexes
setwd('~/Documents/MATLAB/createMitoYeastGEM/')
protImportInfo <- read.delim('../ComplementaryData/proteinImport/importMachineryComponents.tsv',
                             sep = '\t',dec = '.',header = TRUE)[,c(1,2,3,5)]
# Merge TIM23 and PAM into one complex for counting abundance
protImportInfo[which(protImportInfo$Component == 'PAM'),4] <- 'TIM23'

# Merge data from model with and import complex info
combinedData <- merge(protImportInfo,data,by.x = 'UniProtID',by.y = 'UniprotID')
combinedData[which(combinedData$Component == 'Disulfide relay'),4] <- 'MIA'

# Plot comparison between model prediction using estimated kcats and curated kcats to proteomis data at D = 0.1 h-1
dataToPlot_fig3B <- combinedData[,c(4,14,15,16)]
dataToPlot_fig3B <- aggregate(cbind(dataToPlot_fig3B$model_estimated_kcats,dataToPlot_fig3B$model_D_01,
                                    dataToPlot_fig3B$proteomics_D_01), 
                              by = list(Component = dataToPlot_fig3B$Component), FUN = sum)
colnames(dataToPlot_fig3B)[2:4] <- c('Model w/ estimated kcats','Model w/ curated kcats',
                                     'Proteomics')
# Melt data
dataToPlot_fig3B <- melt(dataToPlot_fig3B,id=c('Component'))
# Create a grouped bar plot
fig3B <- ggplot(dataToPlot_fig3B,aes(fill=variable,y=value,x=Component)) +
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
fig3B
ggarrange(fig3B,labels = 'B')
ggsave("~/Documents/Manuscripts/Modeling paper/Figures_latest/fig3B.pdf",
       plot = fig3B,device = "pdf",dpi = 300,width = 7,height = 7,units = 'cm')
#pdf('~/Documents/Manuscripts/Modeling paper/Figures_latest/Figure3/fig3B.pdf')
# Extract data for plot for all dilution rates
dataToPlot <- combinedData[,c(4,15,16,17,18,19,20,21)]

# Compile data for each complex
dataToPlot <- aggregate(cbind(dataToPlot$model_D_01,dataToPlot$proteomics_D_01,dataToPlot$model_D_02,
                              dataToPlot$model_D_028,dataToPlot$model_D_035,dataToPlot$model_0.39,
                              dataToPlot$model_batch),
                        by = list(Component = dataToPlot$Component), FUN=sum)
colnames(dataToPlot)[2:8] <- c('Model D= 0.1 h-1','Proteomics D = 0.1 h-1','Model D = 0.2 h-1',
                               'Model D = 0.28 h-1','Model D = 0.35 h-1','Model D = 0.4 h-1',
                               'Model Batch')
# Melt data
dataToPlot <- melt(dataToPlot,id=c('Component'))

# Divide the data to plot TIM23 and TOM separate from other components
dataToPlot_TIM23andTOM <- dataToPlot[which(dataToPlot$Component == 'TIM23' | dataToPlot$Component == 'TOM' ),]
dataToPlot_remaining <- dataToPlot[-which(dataToPlot$Component == 'TIM23' | dataToPlot$Component == 'TOM' ),]

# create a grouped bar plot
ggplot(dataToPlot,aes(fill=variable,y=value,x=Component)) +
  geom_bar(position = 'dodge',stat = 'identity') +
  theme_bw(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = 'bold'),
        axis.text.x = element_text(angle = 30,size = 12,face = 'bold'),
        axis.text.y = element_text(size = 12,face = 'bold'),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        legend.position = c(0.3,0.8),
        panel.border = element_rect(color = 'black',size = 1,fill = NA),
        panel.background = element_blank()) +
  ylab('Percent of total protein mass (%)') +
  scale_fill_brewer(palette = 'Dark2') +
  ylim(0,0.65) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6))

# Create a grouped bar plot for TIM23 and TOM
fig3D <- ggplot(dataToPlot_TIM23andTOM,aes(fill=variable,y=value,x=Component)) +
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
fig3D
ggsave("",
       plot = fig3D,device = "pdf",dpi = 300,width = 7,height = 7,units = 'cm')

# Create a bar plot for the remaining components
fig3C <- ggplot(dataToPlot_remaining,aes(fill=variable,y=value,x=Component)) +
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
fig3C
ggsave("",
       plot = fig3C,device = "pdf",dpi = 300,width = 7,height = 7,units = 'cm')

# Combine subpanels into 1 graph
ggarrange(fig3C,fig3D,ncol = 2,labels = c('C','D'),common.legend = TRUE,legend='bottom')
ggsave("",
       device = 'pdf',dpi = 300,width = 17.4,height = 7,units = 'cm')
# Create a scatter plot comparing model predictions (estimated kcats) to proteomics data
dataForScatter <- combinedData[,c(4,6,8)]
dataForScatter <- aggregate(cbind(dataForScatter$g.g.estimated.kcats,dataForScatter$g.g.Proteomics),
                            by = list(Component = dataForScatter$Component), FUN = sum)
colnames(dataForScatter)[2:3] <- c('Model','Proteomics')
# log transform values
dataForScatter <- dataForScatter %>%
  rowwise() %>% 
  mutate(Model = log10(Model),
         Proteomics = log10(Proteomics))
# Plot data
fig3X1 <- ggscatter(dataForScatter, x = 'Proteomics', y = 'Model',
          add = 'reg.line', conf.int = FALSE,
          cor.coef = TRUE, cor.method = 'pearson',
          label = 'Component',repel = TRUE,
          xlab = 'in vivo log10(g protein/g total protein mass)',
          ylab = 'in silico log10(g protein/g total protein mass)',
          font.label = c(8,'black','bold'),label.rectangle = TRUE,
          cor.coef.coord = c(-4,-2.4)) +
        theme(axis.text = element_text(size = 8,color = 'black',face = 'bold'))
fig3X1
ggsave('',
       plot = fig3X1,device = "png",dpi = 300)

# Check correlation excluding TIM22, disulfide relay
dataForScatter <- dataForScatter[-c(2,6),]
# Plot data
fig3X2 <- ggscatter(dataForScatter, x = 'Proteomics', y = 'Model',
                    add = 'reg.line', conf.int = FALSE,
                    cor.coef = TRUE, cor.method = 'pearson',
                    label = 'Component',repel = TRUE,
                    xlab = 'in vivo log10(g protein/g total protein mass)',
                    ylab = 'in silico log10(g protein/g total protein mass)',
                    font.label = c(8,'black','bold'),label.rectangle = TRUE,
                    cor.coef.coord = c(-4,-2.4)) +
  theme(axis.text = element_text(size = 8,color = 'black',face = 'bold'))
fig3X2
ggsave('',
       plot = fig3X2,device = "png",dpi = 300)
ggarrange(fig3X1,fig3X2,ncol=2)
ggsave('',
       device = "pdf",dpi = 300,height = 10,width=17.4,units = 'cm')
# Create a scatter plot comparing model predictions (after kcat curation) to proteomics data
dataForScatter <- combinedData[,c(4,7,8)]
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
fig3Y1 <- ggscatter(dataForScatter, x = 'Proteomics', y = 'Model',
                   add = 'reg.line', conf.int = FALSE,
                   cor.coef = TRUE, cor.method = 'pearson',
                   label = 'Component',repel = TRUE,
                   xlab = 'in vivo log10(g protein/g total protein mass)',
                   ylab = 'in silico log10(g protein/g total protein mass)',
                   font.label = c(8,'black','bold'),label.rectangle = TRUE,
                   cor.coef.coord = c(-4,-2.5)) +
  theme(axis.text = element_text(size = 8,color = 'black',face = 'bold'))
fig3Y1

ggsave('',
       plot = fig3Y1,device = "png",dpi = 300)
# Check correlation excluding TIM22, disulfide relay
dataForScatter_modified <- dataForScatter[-c(2,6),]
fig3Y2 <- ggscatter(dataForScatter_modified, x = 'Proteomics', y = 'Model',
                   add = 'reg.line', conf.int = FALSE,
                   cor.coef = TRUE, cor.method = 'pearson',
                   label = 'Component',repel = TRUE,
                   xlab = 'in vivo log10(g protein/g total protein mass)',
                   ylab = 'in silico log10(g protein/g total protein mass)',
                   font.label = c(8,'black','bold'),label.rectangle = TRUE,
                   cor.coef.coord = c(-4,-2.5)) +
  theme(axis.text = element_text(size = 8,color = 'black',face = 'bold'))
fig3Y2
ggsave('',
       plot = fig3Y2,device = "png",dpi = 300)
ggarrange(fig3Y1,fig3Y2,ncol=2)
ggsave('',
       device = "pdf",dpi = 300,height = 10,width=17.4,units = 'cm')

### Create plots related to Fe/S clusters ###
# load data
filename <- 'abundancesFeSbiosynthesis.csv'
data <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

# Convert data from mmol/gDW to mg/gDW
data <- data %>%
  rowwise() %>%
  mutate(model_D_0.1 = mmol.gDW.D...0.1.h.1*MW_kDa*1000,
         proteomics_D_0.1 = mmol.gDW.Proteomics*MW_kDa*1000,
         model_D_0.2 = mmol.gDW.D...0.2.h.1*MW_kDa*1000,
         model_D_0.28 = mmol.gDW.D...0.28.h.1*MW_kDa*1000,
         model_D_0.35 = mmol.gDW.D...0.35.h.1*MW_kDa*1000,
         model_D_0.4 = mmol.gDW.D...0.4.h.1*MW_kDa*1000,
         model_batch = Batch*MW_kDa*1000)

# Select data for plotting (comparing model prediction at D = 0.1 h-1 to proteomics)
dataToPlot <- data[,c(2,12,13)]
colnames(dataToPlot)[2:3] <- c('Model','Proteomics')
# melt data
dataToPlot <- melt(dataToPlot,id=c('Gene'))

# create grouped bar plot (converting to mg/g protein)
fig4C <- ggplot(dataToPlot,aes(fill=variable,y=value/0.46,x=Gene)) +
      geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
            axis.text.x = element_text(angle = 90,size = 8,face = 'bold',color = 'black'),
            axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
            legend.title = element_blank(),
            legend.text = element_text(size = 8,face = 'bold'),
            legend.position = 'bottom',
            panel.border = element_rect(color = 'black',size = 1,fill = NA),
            panel.background = element_blank(),
            axis.ticks = element_line(color = 'black')) +
      ylab('Abundance (mg/g total protein)') +
      scale_fill_brewer(palette = 'Blues') +
      scale_y_continuous(expand = c(0,0),limits = c(0,1.15))
fig4C
ggsave('',
       plot = fig4C,device = "pdf",dpi = 300,height = 6,width = 8.5,units = 'cm')

# select data for plotting (all dilution rates)
dataToPlot <- data[,c(2,12,13,14,15,16,17,18)]
colnames(dataToPlot)[2:8] <- c('D = 0.1 h-1','Proteomics D = 0.1 h-1','D = 0.2 h-1',
                               'D = 0.28 h-1','D = 0.35 h-1','D = 0.4 h-1',
                               'Batch')
# melt data
dataToPlot <- melt(dataToPlot,id=c('Gene'))

# Separate Mge1 and Grx5 from the rest since they have a much higher abundance
dataToPlot_Mge1andGrx5 <- dataToPlot[which(dataToPlot$Gene == 'MGE1' | dataToPlot$Gene == 'GRX5'),]
dataToPlot_remaining <- dataToPlot[-which(dataToPlot$Gene == 'MGE1' | dataToPlot$Gene == 'GRX5'),]

# create grouped bar plot (convert mg/gDW to mg/g total protein)
ggplot(dataToPlot,aes(fill=variable,y=value/0.46,x=Gene)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = 'bold',color = 'black'),
        axis.text.x = element_text(angle = 30,size = 12,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 12,face = 'bold',color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        legend.position = c(0.3,0.8),
        panel.border = element_rect(color = 'black',fill = NA),
        panel.background = element_blank(),
        axis.ticks = element_line(color = 'black')) +
  ylab('Abundance (mg/g total protein)') +
  scale_fill_brewer(palette = 'Dark2') +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.175))

# create grouped bar plot for Mge1 and Grx5
fig4G <- ggplot(dataToPlot_Mge1andGrx5,aes(fill=variable,y=value/0.46,x=Gene)) +
              geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
              theme_bw() +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
                    axis.text.x = element_text(angle = 0,size = 8,face = 'bold',color = 'black'),
                    axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 8,face = 'bold',color = 'black'),
                    legend.position = c(0.3,0.75),
                    panel.border = element_rect(color = 'black',fill = NA),
                    panel.background = element_blank(),
                    axis.ticks = element_line(color = 'black')) +
              ylab('Abundance (mg/g protein)') +
              scale_fill_brewer(palette = 'Blues') +
              scale_y_continuous(expand = c(0,0),limits = c(0,1.15))
fig4G
ggsave('',
       plot = fig4G,device = "pdf",dpi = 300,width = 5.4,height = 7,units = 'cm')

# create grouped bar plot for remaining proteins
fig4F <- ggplot(dataToPlot_remaining,aes(fill=variable,y=value/0.46,x=Gene)) +
              geom_bar(position = position_dodge(0.85),stat = 'identity',color = 'black') +
              theme_bw() +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
                    axis.text.x = element_text(angle = 90,size = 8,face = 'bold',color = 'black'),
                    axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 8,face = 'bold'),
                    legend.position = 'none',
                    panel.border = element_rect(color = 'black',fill = NA),
                    panel.background = element_blank(),
                    axis.ticks = element_line(color = 'black')) +
              ylab('Abundance (mg/g protein)') +
              scale_fill_brewer(palette = 'Blues') +
              scale_y_continuous(expand = c(0,0),limits = c(0,0.115))
fig4F
ggsave('',
       plot = fig4F,device = "pdf",dpi = 300,height=7,width = 12,units = 'cm')
ggarrange(fig4F,fig4G,ncol = 2,labels = c('F','G'),common.legend = TRUE,legend='bottom')
ggsave('',
       device = "pdf",dpi = 300,height=7,width = 17.4,units = 'cm')
# Divide the data by stage of synthesis
dataByStage <- data[,c(10,12,13,14,15,16,17,18)]
colnames(dataByStage)[2:8] <- c('Model D= 0.1 h-1','Proteomics D = 0.1 h-1','Model D = 0.2 h-1',
                                'Model D = 0.28 h-1','Model D = 0.35 h-1','Model D = 0.4 h-1',
                                'Model Batch')
dataByStage <- aggregate(cbind(dataByStage$`Model D= 0.1 h-1`,dataByStage$`Proteomics D = 0.1 h-1`,
                               dataByStage$`Model D = 0.2 h-1`,dataByStage$`Model D = 0.28 h-1`,
                               dataByStage$`Model D = 0.35 h-1`,dataByStage$`Model D = 0.4 h-1`,
                               dataByStage$`Model Batch`),
             by = list(Stage = dataByStage$Stage), FUN=sum)
colnames(dataByStage)[2:8] <- c('Model D= 0.1 h-1','Proteomics D = 0.1 h-1','Model D = 0.2 h-1',
                                'Model D = 0.28 h-1','Model D = 0.35 h-1','Model D = 0.4 h-1',
                                'Model Batch')
dataByStage <- melt(dataByStage,id=c('Stage'))

# Create a grouped bar plot (convert g/gDW to g/g protein)
ggplot(dataByStage,aes(fill=variable,y=value/0.46,x=Stage)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 11,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 12,face = 'bold',color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        legend.position = c(0.75,0.8),
        panel.border = element_rect(color = 'black',fill = NA,size = 0.75),
        panel.background = element_blank(),
        axis.ticks = element_line(color = 'black')) +
  ylab('Abundance (g/g protein)') +
  scale_fill_brewer(palette = 'Dark2') +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.675))

# Load data for model proteins containing Fe/S clusters
filename <- 'abundancesFeSproteins.csv'
data <- read.csv(filename,header = TRUE,sep = ';',dec = ',')
# Convert data from mmol/gDW to mg/g protein
data <- data %>%
  rowwise() %>%
  mutate(Chemostat = ((mmol.gDW.D...0.1.h.1*MW_kDa)/0.46)*1000,
         Batch = ((mmol.gDW.batch*MW_kDa)/0.46)*1000)
  
dataForScatter <- data[,c(2,10,11)]

# Plot scatter plot
fig4X <- ggplot(dataForScatter,aes(x=Chemostat,y=Batch)) +
              geom_point() +
              geom_abline(slope = 1,intercept = 0, size = 1) +
              geom_label_repel(aes(label = Protein),point.padding = 0.5,segment.color = 'black') +
  #geom_text(aes(label = Protein),hjust = 0,vjust = 0)
              xlab('mg/g protein growth rate 0.1 h-1') +
              ylab('mg/g protein max growth rate') +
              theme_pubr() +
              theme(axis.text = element_text(size = 12,face = 'bold',color = 'black'),
                    axis.title = element_text(size = 12,face = 'bold',color = 'black')) +
              scale_x_continuous(limits = c(0,5)) +
              scale_y_continuous(limits = c(0,3.5))
fig4X
ggsave('',
       plot = fig4X,device = "png",dpi = 300)
# Divide data by higher and lower abundance
dataForScatter1 <- dataForScatter[which(dataForScatter$Batch > 0.13),]
  
# Plot scatter plot
fig4X <- ggplot(dataForScatter1,aes(x=Chemostat,y=Batch)) +
              geom_point() +
              geom_abline(slope = 1,intercept = 0, size = 1) +
              geom_label_repel(aes(label = Protein),point.padding = 0.5,segment.color = 'black') +
  #geom_text(aes(label = Protein),hjust = 0,vjust = 0)
              xlab('in silico abundance (mg/g protein) at growth rate 0.1 h-1') +
              ylab('in silico abundance (mg/g protein) max growth rate') +
              theme_pubr() +
              theme(axis.text = element_text(size = 8,face = 'bold',color = 'black'),
                    axis.title = element_text(size = 8,face = 'bold',color = 'black')) +
              scale_x_continuous(limits = c(0,5)) +
              scale_y_continuous(limits = c(0,3.2))
fig4X
ggsave('',
       plot = fig4X,device = "png",dpi = 300)
# Collect data for proteins with lower abundance
dataForScatter2 <- dataForScatter[which(dataForScatter$Batch < 0.13),]
# Create another scatter plot
fig4Y <- ggplot(dataForScatter2,aes(x=Chemostat,y=Batch)) +
              geom_point() +
              geom_abline(slope = 1,intercept = 0, size = 1) +
              geom_label_repel(aes(label = Protein),point.padding = 0.5,segment.color = 'black') +
              #geom_text(aes(label = Protein),hjust = 0,vjust = 0)
              xlab('in silico abundance (mg/g protein) at growth rate 0.1 h-1') +
              ylab('in silico abundance (mg/g protein) at max growth rate') +
              theme_pubr() +
              theme(axis.text = element_text(size = 8,face = 'bold',color = 'black'),
                    axis.title = element_text(size = 8,face = 'bold',color = 'black')) +
              scale_x_continuous(limits = c(0,0.171)) +
              scale_y_continuous(limits = c(0,0.122))# +
              #xlim(0,5.5e-5) +
              #ylim(0,5.5e-5)
fig4Y
ggsave('',
       plot = fig4Y,device = "png",dpi = 300)
ggarrange(fig4X,fig4Y,ncol=2)
ggsave('',
       device = "pdf",dpi = 300,height = 10,width = 17.4,units = 'cm')

### Create plot for Fe/S cluster requirements
# Load data
filename <- 'FeSrequirement.csv'
data <- read.csv(filename,header = TRUE,sep = ';',dec = ',')
colnames(data)[2:7] <- c('D = 0.1 h-1','D = 0.2 h-1','D = 0.28 h-1',
                         'D = 0.35 h-1','D = 0.4 h-1','Batch')
data <- melt(data,id=c('Cofactor'))
# Create grouped bar plot
fig4D <- ggplot(data,aes(fill=variable,y=value,x=Cofactor)) +
              geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
              theme_bw() +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
                    axis.text.x = element_text(angle = 30,size = 8,face = 'bold',color = 'black'),
                    axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 8,face = 'bold'),
                    legend.position = c(0.3,0.8),
                    panel.border = element_rect(color = 'black',size = 1,fill = NA),
                    panel.background = element_blank()) +
              scale_y_continuous(labels = scientific, expand = c(0,0),limits=c(0,7.25e-05)) +
              scale_fill_brewer(palette = 'Blues') +
              ylab('Abundance (mmol/gDW)')
fig4D
ggsave('',
       plot = fig4D,device = "png",dpi = 300)

### Create plots for Cofactor (other than Fe/S clusters) requirement
# Load data
filename <- 'cofactorRequirements.csv'
data <- read.csv(filename,header = TRUE,sep = ';',dec = ',')
colnames(data)[2:7] <- c('D = 0.1 h-1','D = 0.2 h-1','D = 0.28 h-1',
                         'D = 0.35 h-1','D = 0.4 h-1','Batch')
data_siroheme <- data[which(data$cofactor == 'Siroheme'),]
data_siroheme <- melt(data_siroheme,id=c('cofactor'))
data <- data[-which(data$cofactor == 'Siroheme'),]
data <- melt(data,id=c('cofactor'))

# Create grouped bar plot
fig4A <- ggplot(data,aes(fill=variable,y=value,x=cofactor)) +
              geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
              theme_bw() +
              theme(axis.title.x = element_blank(),
                    axis.title.y = element_text(size = 8,face = 'bold'),
                    axis.text.x = element_text(angle = 30,size = 8,face = 'bold',color = 'black'),
                    axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 8,face = 'bold'),
                    legend.position = c(0.25,0.8),
                    panel.border = element_rect(color = 'black',size = 1,fill = NA),
                    panel.background = element_blank()) +
              scale_y_continuous(labels = scientific,expand = c(0,0),limits = c(0,1.35e-04)) +
              scale_fill_brewer(palette = 'Blues') +
              ylab('Abundance (mmol/gDW)')
fig4A
ggsave('',
       plot = fig4A,device = "pdf",dpi = 300,height = 5,width = 5.8,units = 'cm')
# Create a separate plot for siroheme
fig4A2 <- ggplot(data_siroheme,aes(fill=variable,y=value,x=cofactor)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = 'bold'),
        axis.text.x = element_text(angle = 0,size = 8,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 8,face = 'bold'),
        legend.position = c(0.2,0.8),
        panel.border = element_rect(color = 'black',size = 1,fill = NA),
        panel.background = element_blank()) +
  scale_y_continuous(labels = scientific,expand = c(0,0),limits = c(0,5e-08)) +
  scale_fill_brewer(palette = 'Blues') +
  ylab('Abundance (mmol/gDW)')
fig4A2
ggsave('',
       plot = fig4A2,device = "png",dpi = 300)
ggarrange(fig4A,fig4A2,fig4D,ncol = 3,labels = c('A','B','C'),common.legend = TRUE,legend='bottom')
ggsave("",
              device = "pdf",dpi = 300,width = 17.4,height = 5.4,units = 'cm')
##### Plot P/O and Yxs at increasing dilution rate
# load data
filename <- 'POandBiomassYield.txt'
data <- read.delim(filename,header = TRUE)
colnames(data) <- c('growthRate','P/O','Biomass yield (gDW/g glucose)')
data <- melt(data,id=c('growthRate'))

# Create a grouped bar plot
fig2E <- ggplot(data,aes(fill=variable,y=value,x=growthRate)) +
              geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
              theme_bw() +
              theme(axis.title.x = element_text(size = 8,face = 'bold',color = 'black'),
                    axis.title.y = element_text(size = 8,face = 'bold'),
                    axis.text.x = element_text(angle = 0,size = 8,face = 'bold',color = 'black'),
                    axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
                    legend.title = element_blank(),
                    legend.text = element_text(size = 8,face = 'bold'),
                    legend.position = 'top',
                    panel.border = element_rect(color = 'black',size = 1,fill = NA),
                    panel.background = element_blank()) +
              scale_y_continuous(expand = c(0,0),limits = c(0,1.1)) +
              scale_fill_brewer(palette = 'Blues') +
              ylab('P/O, Biomass yield (gDW/g glucose)') +
              xlab('Growth rate (h-1)')
fig2E
ggsave('',
       plot = fig2E,device = "pdf",dpi = 300,width = 6,height = 6,units = 'cm')


