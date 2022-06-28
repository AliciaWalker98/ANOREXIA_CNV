

args=commandArgs(trailingOnly = T)
cohort = args[1]
refvec=args[2]
pcvec=args[3]
popstr=args[4]


# examples for testing
# cohort="MND_B6_GSResults_20210517"
# refvec="/scratch/60days/uqtlin5/RTP_Pipeline_test/Temp//MND_B6_GSResults_20210517/PCA/1000G_phase3_20130502_combined_snpsonly.05.common_pca3"
# pcvec="/scratch/60days/uqtlin5/RTP_Pipeline_test/Temp//MND_B6_GSResults_20210517/PCA/cleaned_flipped_MND_B6_GSResults_20210517.05.common_pca3.proj.eigenvec"
# popstr="/afm01/UQ/Q3046/Reference//1000GP_Phase3.sample"



library(ggplot2)
library(ggpubr)
str = read.table(popstr, header=T)
ref = read.table(paste0(refvec,".eigenvec" ))
pc = read.table(pcvec)

ref$str = str[match(ref$V2, str$ID),"GROUP"]
pc$str <- "Samples"
colnames(ref) = c("FID", "IID", "PC1", "PC2", "PC3", "str")
colnames(pc) = c("FID", "IID", "PC1", "PC2", "PC3", "str")

pc$str = as.character(pc$str)
pc$FID = as.character(pc$FID)
pc$IID = as.character(pc$IID)
ref$str = as.character(ref$str)
ref$FID = as.character(ref$FID)
ref$IID = as.character(ref$IID)

data = rbind((ref), (pc))



data$str = factor(data$str, levels = c("AFR", "AMR", "EAS", "EUR", "SAS",  "Samples"))

cbPalette <- c("#999999",    
				"#F0E442", 	
				"#E69F00",   
				"#009E73",  
				"#FC4E07", 
				"#0072B2"
				)	




#####################################################

## ancestry parameters

### EAS
pc1_SD_a =   sd(data[data$str=="EAS","PC1"])
pc2_SD_a =   sd(data[data$str=="EAS","PC2"])   
pc3_SD_a =   sd(data[data$str=="EAS","PC3"])   
pc1_mn_a = mean(data[data$str=="EAS","PC1"])  
pc2_mn_a = mean(data[data$str=="EAS","PC2"])     
pc3_mn_a = mean(data[data$str=="EAS","PC3"])     


### EUR
pc1_SD_e =   sd(data[data$str=="EUR","PC1"])
pc2_SD_e =   sd(data[data$str=="EUR","PC2"])   
pc3_SD_e =   sd(data[data$str=="EUR","PC3"])   
pc1_mn_e = mean(data[data$str=="EUR","PC1"])  
pc2_mn_e = mean(data[data$str=="EUR","PC2"])     
pc3_mn_e = mean(data[data$str=="EUR","PC3"])     


### AFR
pc1_SD_f =   sd(data[data$str=="AFR","PC1"])
pc2_SD_f =   sd(data[data$str=="AFR","PC2"])   
pc3_SD_f =   sd(data[data$str=="AFR","PC3"])   
pc1_mn_f = mean(data[data$str=="AFR","PC1"])  
pc2_mn_f = mean(data[data$str=="AFR","PC2"])     
pc3_mn_f = mean(data[data$str=="AFR","PC3"])     


### SAS
pc1_SD_s =   sd(data[data$str=="SAS","PC1"])
pc2_SD_s =   sd(data[data$str=="SAS","PC2"])   
pc3_SD_s =   sd(data[data$str=="SAS","PC3"])   
pc1_mn_s = mean(data[data$str=="SAS","PC1"])  
pc2_mn_s = mean(data[data$str=="SAS","PC2"])     
pc3_mn_s = mean(data[data$str=="SAS","PC3"])     


### AMR
pc1_SD_m =   sd(data[data$str=="AMR","PC1"])
pc2_SD_m =   sd(data[data$str=="AMR","PC2"])   
pc3_SD_m =   sd(data[data$str=="AMR","PC3"])   
pc1_mn_m = mean(data[data$str=="AMR","PC1"])  
pc2_mn_m = mean(data[data$str=="AMR","PC2"])     
pc3_mn_m = mean(data[data$str=="AMR","PC3"])     


## set boundaries
n = 6

x1_asn = pc1_mn_a - 8 * pc1_SD_a
x2_asn = pc1_mn_a + 8 * pc1_SD_a
y1_asn = pc2_mn_a - 8 * pc2_SD_a
y2_asn = pc2_mn_a + 8 * pc2_SD_a
z1_asn = pc3_mn_a - n * pc3_SD_a
z2_asn = pc3_mn_a + n * pc3_SD_a


x1_eur = pc1_mn_e - n * pc1_SD_e
x2_eur = pc1_mn_e + n * pc1_SD_e
y1_eur = pc2_mn_e - n * pc2_SD_e
y2_eur = pc2_mn_e + n * pc2_SD_e
z1_eur = pc3_mn_e - n * pc3_SD_e
z2_eur = pc3_mn_e + n * pc3_SD_e


x1_afr = pc1_mn_f - n * pc1_SD_f
x2_afr = pc1_mn_f + n * pc1_SD_f
y1_afr = pc2_mn_f - n * pc2_SD_f
y2_afr = pc2_mn_f + n * pc2_SD_f
z1_afr = pc3_mn_f - n * pc3_SD_f
z2_afr = pc3_mn_f + n * pc3_SD_f


x1_sas = pc1_mn_s - n * pc1_SD_s
x2_sas = pc1_mn_s + n * pc1_SD_s
y1_sas = pc2_mn_s - 4 * pc2_SD_s
y2_sas = pc2_mn_s + 4 * pc2_SD_s
z1_sas = pc3_mn_s - 4 * pc3_SD_s
z2_sas = pc3_mn_s + 4 * pc3_SD_s

x1_amr = pc1_mn_m - 2 * pc1_SD_m
x2_amr = pc1_mn_m + 3 * pc1_SD_m
y1_amr = pc2_mn_m - 2 * pc2_SD_m
y2_amr = pc2_mn_m + 2 * pc2_SD_m
z1_amr = pc3_mn_m - 2 * pc3_SD_m
z2_amr = pc3_mn_m + 3 * pc3_SD_m



#####################################################

## figure option 1

fig <- ggplot(data=data, aes(x=PC1, y=PC2, color=str)) + 
  geom_point(size=0.8)  + 
  scale_color_manual(values=cbPalette) +
  ylim(-0.04, 0.04) +
  xlim(-0.02, 0.06) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
	panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())




framed.fig =   fig + 
  geom_vline(xintercept = x1_eur, color='#009E73', size=0.3) +
  geom_vline(xintercept = x2_eur, color='#009E73', size=0.3) +
  geom_hline(yintercept = y1_eur, color='#009E73', size=0.3) +
  geom_hline(yintercept = y2_eur, color='#009E73', size=0.3) +
  geom_vline(xintercept = x1_asn, color='#E69F00', size=0.3) +
  geom_vline(xintercept = x2_asn, color='#E69F00', size=0.3) +
  geom_hline(yintercept = y1_asn, color='#E69F00', size=0.3) +
  geom_hline(yintercept = y2_asn, color='#E69F00', size=0.3) +
  geom_vline(xintercept = x1_afr, color='#999999', size=0.3) +
  geom_vline(xintercept = x2_afr, color='#999999', size=0.3) +
  geom_hline(yintercept = y1_afr, color='#999999', size=0.3) +
  geom_hline(yintercept = y2_afr, color='#999999', size=0.3) +
  geom_vline(xintercept = x1_sas, color='#FC4E07', size=0.3) +
  geom_vline(xintercept = x2_sas, color='#FC4E07', size=0.3) +
  geom_hline(yintercept = y1_sas, color='#FC4E07', size=0.3) +
  geom_hline(yintercept = y2_sas, color='#FC4E07', size=0.3) +
  geom_vline(xintercept = x1_amr, color='#F0E442', size=0.3) +
  geom_vline(xintercept = x2_amr, color='#F0E442', size=0.3) +
  geom_hline(yintercept = y1_amr, color='#F0E442', size=0.3) +
  geom_hline(yintercept = y2_amr, color='#F0E442', size=0.3)



fig2 <- ggplot(data=data, aes(x=PC1, y=PC3, color=str)) + 
  geom_point(size=0.8)  + 
  scale_color_manual(values=cbPalette)+
  ylim(-0.06, 0.06) +
  xlim(-0.02, 0.04) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


framed.fig2 =   fig2 + 
  geom_vline(xintercept = x1_eur, color='#009E73', size=0.3) +
  geom_vline(xintercept = x2_eur, color='#009E73', size=0.3) +
  geom_hline(yintercept = z1_eur, color='#009E73', size=0.3) +
  geom_hline(yintercept = z2_eur, color='#009E73', size=0.3) +
  geom_vline(xintercept = x1_asn, color='#E69F00', size=0.3) +
  geom_vline(xintercept = x2_asn, color='#E69F00', size=0.3) +
  geom_hline(yintercept = z1_asn, color='#E69F00', size=0.3) +
  geom_hline(yintercept = z2_asn, color='#E69F00', size=0.3) +
  geom_vline(xintercept = x1_afr, color='#999999', size=0.3) +
  geom_vline(xintercept = x2_afr, color='#999999', size=0.3) +
  geom_hline(yintercept = z1_afr, color='#999999', size=0.3) +
  geom_hline(yintercept = z2_afr, color='#999999', size=0.3) +
  geom_vline(xintercept = x1_sas, color='#FC4E07', size=0.3) +
  geom_vline(xintercept = x2_sas, color='#FC4E07', size=0.3) +
  geom_hline(yintercept = z1_sas, color='#FC4E07', size=0.3) +
  geom_hline(yintercept = z2_sas, color='#FC4E07', size=0.3) +
  geom_vline(xintercept = x1_amr, color='#F0E442', size=0.3) +
  geom_vline(xintercept = x2_amr, color='#F0E442', size=0.3) +
  geom_hline(yintercept = z1_amr, color='#F0E442', size=0.3) +
  geom_hline(yintercept = z2_amr, color='#F0E442', size=0.3)


framed.fig.c = ggarrange(framed.fig, framed.fig2, ncol = 2, nrow =  1)


#ggsave(framed.fig.c, file=paste0("PCA/",cohort, "_PC_plot_meanSD_framed.pdf"), width = 12, height = 6)



#ggsave(framed.fig.c, file = paste0("PCA/", cohort, "_PC_plot_meanSD_framed.png"), width = 12, height = 6)
png(paste0("PCA/", cohort, "_PC_plot_meanSD_framed.png"),res=100, width=12, height=6, units='in', bg="white",type=c("cairo"))
framed.fig.c
dev.off()
# fig.c = ggarrange(fig, fig2, ncol = 2, nrow =  1)
# ggsave(fig.c, file=paste0("../Data/GRP_B8_GSResults_20211028/Imputation/PCA/",cohort, "_PC_plot_noframe.png"), width = 12, height = 6)

#####################################################3


### Assign population 
pc$population = "other"

pc[which(pc$PC1>x1_amr & pc$PC1<x2_amr  & pc$PC2>y1_amr & pc$PC2<y2_amr  & pc$PC3>z1_amr & pc$PC3<z2_amr), "population"] = "AMR"
pc[which(pc$PC1>x1_sas & pc$PC1<x2_sas  & pc$PC2>y1_sas & pc$PC2<y2_sas  & pc$PC3>z1_sas & pc$PC3<z2_sas), "population"] = "SAS"
pc[which(pc$PC1>x1_asn & pc$PC1<x2_asn  & pc$PC2>y1_asn & pc$PC2<y2_asn  & pc$PC3>z1_asn & pc$PC3<z2_asn), "population"] = "EAS"
pc[which(pc$PC1>x1_eur & pc$PC1<x2_eur  & pc$PC2>y1_eur & pc$PC2<y2_eur  & pc$PC3>z1_eur & pc$PC3<z2_eur), "population"] = "EUR"
pc[which(pc$PC1>x1_afr & pc$PC1<x2_afr  & pc$PC2>y1_afr & pc$PC2<y2_afr  & pc$PC3>z1_afr & pc$PC3<z2_afr), "population"] = "AFR"

write.table(pc, file=paste0(cohort , "_ancestry_assigned_with_meanSD_method.txt"), row.names=F, col.names=T, quote=F, sep="\t")






#####################################################3





#### figure option 2

if(nrow(pc[which(pc$population == "EUR"),])>0) {
  fig = fig + 
    geom_vline(xintercept = x1_eur, color='#009E73', size=0.3) +
    geom_vline(xintercept = x2_eur, color='#009E73', size=0.3) +
    geom_hline(yintercept = y1_eur, color='#009E73', size=0.3) +
    geom_hline(yintercept = y2_eur, color='#009E73', size=0.3) 
}

if(nrow(pc[which(pc$population == "EAS"),])>0) {
  fig = fig + 
    geom_vline(xintercept = x1_asn, color='#E69F00', size=0.3) +
    geom_vline(xintercept = x2_asn, color='#E69F00', size=0.3) +
    geom_hline(yintercept = y1_asn, color='#E69F00', size=0.3) +
    geom_hline(yintercept = y2_asn, color='#E69F00', size=0.3) 
}

if(nrow(pc[which(pc$population == "AFR"),])>0) {
  fig = fig + 
    geom_vline(xintercept = x1_afr, color='#999999', size=0.3) +
    geom_vline(xintercept = x2_afr, color='#999999', size=0.3) +
    geom_hline(yintercept = y1_afr, color='#999999', size=0.3) +
    geom_hline(yintercept = y2_afr, color='#999999', size=0.3)
}

if(nrow(pc[which(pc$population == "SAS"),])>0) {
  fig2 = fig2 + 
    geom_vline(xintercept = x1_sas, color='#FC4E07', size=0.3) +
    geom_vline(xintercept = x2_sas, color='#FC4E07', size=0.3) +
    geom_hline(yintercept = z1_sas, color='#FC4E07', size=0.3) +
    geom_hline(yintercept = z2_sas, color='#FC4E07', size=0.3)
}

if(nrow(pc[which(pc$population == "AMR"),])>0) {
  fig2 = fig2 + 
    geom_vline(xintercept = x1_amr, color='#F0E442', size=0.3) +
    geom_vline(xintercept = x2_amr, color='#F0E442', size=0.3) +
    geom_hline(yintercept = z1_amr, color='#F0E442', size=0.3) +
    geom_hline(yintercept = z2_amr, color='#F0E442', size=0.3)
}

fig.c = ggarrange(fig, fig2, ncol = 2, nrow =  1)


#ggsave(fig.c, file=paste0("PCA/",cohort, "_PC_plot_meanSD.pdf"), width = 12, height = 6)

ggsave(fig.c, file = paste0( cohort, "_PC_plot_meanSD.png"), width = 12, height = 6)


##########################################




### side by side fig

colnames(ref)[6] = "population"
ref$str="1000G"
ref = ref[,c(1:5, 7, 6) ]
data.to.plot = rbind(pc, ref)
colnames(data.to.plot)[6:7] = c("sample", "POP") 



cbPalette.for.all = data.frame(
  ancestry = c("AFR", "AMR", "EAS", "EUR", "SAS", "other"), 
  color = c("#999999",    
            "#F0E442", 	
            "#E69F00",   
            "#009E73",  
            "#FC4E07",
            '#0072B2')	)

cbPalette = cbPalette.for.all[cbPalette.for.all$ancestry %in% unique(pc$population),]
cbPalette.1000g = cbPalette.for.all[c(1:5),]

data.to.plot$POP = factor(data.to.plot$POP, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "other"))

sep.fig1 = ggplot(data=data.to.plot) + 
  geom_point(data = subset(data.to.plot, sample == "Samples"),  aes(x=PC1, y=PC2, color=POP), size = 1)  + 
  scale_color_manual(values=c(cbPalette$color))+
  ylim(-0.06, 0.06) +
  xlim(-0.02, 0.04) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) + guides(alpha = FALSE)



sep.fig2 = ggplot(data=data.to.plot) + 
  geom_point(data = subset(data.to.plot, sample == "Samples"),  aes(x=PC1, y=PC3, color=POP ), size = 1)  + 
  scale_color_manual(values=c(cbPalette$color))+
  ylim(-0.06, 0.06) +
  xlim(-0.02, 0.04) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) + guides(alpha = FALSE)

sep.fig3 = ggplot(data=data.to.plot) + 
  geom_point(data = subset(data.to.plot, sample == "Samples"),  aes(x=PC2, y=PC3, color=POP), size = 1)  + 
  scale_color_manual(values=c(cbPalette$color))+
  ylim(-0.06, 0.06) +
  xlim(-0.04, 0.04) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) + guides(alpha = FALSE)


sep.fig4 =  ggplot(data=data.to.plot) + 
  geom_point(data = subset(data.to.plot, sample == "1000G"),  aes(x=PC1, y=PC2, color=POP), shape = 3, size = 0.5)  + 
  scale_color_manual(values=c(cbPalette.for.all$color))+
  ylim(-0.06, 0.06) +
  xlim(-0.02, 0.04) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) 



sep.fig5 = ggplot(data=data.to.plot) + 
  geom_point(data = subset(data.to.plot, sample == "1000G"),  aes(x=PC1, y=PC3, color=POP), shape = 3 , size = 0.5)  + 
  scale_color_manual(values=c(cbPalette.for.all$color))+
  ylim(-0.06, 0.06) +
  xlim(-0.02, 0.04) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) 


sep.fig6 = ggplot(data=data.to.plot) + 
  geom_point(data = subset(data.to.plot, sample == "1000G"),  aes(x=PC2, y=PC3, color=POP), shape = 3 , size = 0.5)  + 
  scale_color_manual(values=c(cbPalette.for.all$color))+
  ylim(-0.06, 0.06) +
  xlim(-0.04, 0.04) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) 



sep.fig = ggarrange(sep.fig1, sep.fig2, sep.fig3, sep.fig4, sep.fig5, sep.fig6,  ncol = 3, nrow = 2)

ggsave(sep.fig, filename = paste0( cohort, "_PC_plot_meanSD_compared_1000G.png"), width = 18, height = 10)









