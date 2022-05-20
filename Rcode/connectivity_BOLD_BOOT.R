library(readxl)
library(HDInterval)
library(parallel)
library(doParallel)

library(boot)
library(bayesboot)
library(ggplot2)
library(plyr)
library(lsmeans)
library(raster)
library(lattice)
library(magrittr)
library(car)
library(RColorBrewer)
library(factoextra)
library(FactoMineR)
library(clue)
library(NbClust)
library(ggrepel)
library(circlize)
library(stringr)
library(lme4)
library(merTools)
library(splines)
library(lsmeans)


########
#### Functions
########


boot_Time <- function(d){

	set.seed(123)
	model <- loess(R ~ Time, degree=2, span=1, d)
	newdat <- data.frame("Time"=seq(33,85,1))
	newdat$pred <- predict(model, newdat)
	return(c(newdat$pred)) 
}



#################
### Load data ###
#################

R_FUS <- lapply(list.files(path = "D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\All_and_scripts\\dFC_results\\v2\\FUS", pattern=glob2rx("R_*.txt"), full.names=T), read.table, header = T, sep = "", dec = ",")
names(R_FUS) <- gsub('.txt', "", list.files(path = "D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\All_and_scripts\\dFC_results\\v2\\FUS", pattern=glob2rx("R_*.txt")))

Data <- read_excel("D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\All_and_scripts\\rsfMRI_BOLD.xlsx",1)

Data <- Data[!(Data$Suj == "7u1" & Data$Run == "Run2" ),]
Data <- Data[!(Data$Suj == "HS1_s2018" & Data$Run == "Run2" ),]

###########
### Run ###
###########

res2 = data.frame(ROI1=numeric(), ROI2=numeric(), N_diff=numeric(), model_glm_p=numeric(), model_glm_coef=numeric(), stringsAsFactors=FALSE)

cl=makeCluster(detectCores()-2)    
registerDoParallel(cl) 

for (x in 1:140){

	for (j in 1:140){

		if (x >= j) {next}

		## Build data

		df_fus = data.frame("Time"=numeric(), "R"=numeric(), "ID"=character(), stringsAsFactors=FALSE)
		for(z in 1:length(R_FUS)){ df_fus[nrow(df_fus) +1,] = data.frame(as.numeric(gsub("m", '', gsub("e", '', gsub("_", '', str_sub(names(R_FUS)[[z]], -5))))),
										R_FUS[[z]][x,j],										
										gsub(as.numeric(gsub("m", '', gsub("e", '', gsub("_", '', str_sub(names(R_FUS)[[z]], -5))))), '', names(R_FUS)[[z]]))
		}

		df_cont = subset(Data, Group == "Control")
		df_cont$Suj2 = paste(df_cont$Suj, "_", df_cont$Run, sep="")
		df_cont = df_cont[c(x+4, j+4, ncol(df_cont))]
		colnames(df_cont)[c(1,2)] = c("ROI1", "ROI2")
		df_cont <- ddply(df_cont, .(Suj2), summarize, R = cor.test(ROI1, ROI2)$estimate)

		

		## Boot

		out <- bayesboot(df_fus, boot_Time, R=500, .parallel=TRUE)

		out_clean <- data.frame("Time"=seq(33,85,1), "mean"=colMeans(out, na.rm=T))
		for (z in 1:ncol(out)){out_clean$HDI_L[z] <- hdi(out[,z], credMass=0.95)[[1]]}
		for (z in 1:ncol(out)){out_clean$HDI_H[z] <- hdi(out[,z], credMass=0.95)[[2]]}

		out_clean$sig = ifelse(out_clean$HDI_L > hdi(df_cont$R)[[2]], 1, 
					ifelse(out_clean$HDI_H < hdi(df_cont$R)[[1]], 1, 0))

		model_glm = glm(sig ~ Time, out_clean, family="binomial")

		res2[nrow(res2) +1,] = data.frame(x, j, sum(out_clean$sig), Anova(model_glm)$`Pr(>Chisq)`, coef(model_glm)[[2]])

		write.table(res2, file = "D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\All_and_scripts\\Boot.txt", row.names=FALSE, sep="\t", dec=",")

	#	plot(df_fus$Time, df_fus$R)
	#	lines(out_clean$Time, out_clean$mean, col="blue")
	#	lines(out_clean$Time, out_clean$HDI_H, col="blue")
	#	lines(out_clean$Time, out_clean$HDI_L, col="blue")

	#	out_clean$pred = predict(model_glm, type="response")
	#	plot(out_clean$Time, out_clean$pred, type="o")
	#	test = subset(out_clean, pred > 0.5)

		cat("\014")
		cat(paste("ROI-1 = ", x, "\nROI-2 = ", j, sep=""))

	}
}

stopCluster(cl)


##################################
### Select significant results ###
##################################

total = read.table("D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\All_and_scripts\\Boot.txt", header = T, sep="\t", dec = ",")

total2 = subset(total, model_glm_p <= 0.05 & model_glm_coef < 0)

total2$FWE_glm <- p.adjust(total2$model_glm_p, method ="bonferroni", n = 9660)


#################
## Whole-Brain ##

total = read.table("D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\All_and_scripts\\New.txt", header = T, sep="\t", dec = ",")

total$FWE_t <- p.adjust(total$model_t, method ="bonferroni", n = (140*140))
total$FWE_glm <- p.adjust(total$model_p, method ="bonferroni", n = (140*140))

total2 = subset(total, FWE_t <= 0.05 & FWE_glm <= 0.05 & model_p_coef < 0)

## OR

total = read.table("D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\All_and_scripts\\New.txt", header = T, sep="\t", dec = ",")

total = subset(total, ROI1 != 37) ; total = subset(total, ROI2 != 37)
total = subset(total, ROI1 != 40) ; total = subset(total, ROI2 != 40)
total = subset(total, ROI1 != 46) ; total = subset(total, ROI2 != 46)
total = subset(total, ROI1 != 65) ; total = subset(total, ROI2 != 65)
total = subset(total, ROI1 != 107) ; total = subset(total, ROI2 != 107)
total = subset(total, ROI1 != 110) ; total = subset(total, ROI2 != 110)
total = subset(total, ROI1 != 116) ; total = subset(total, ROI2 != 116)
total = subset(total, ROI1 != 135) ; total = subset(total, ROI2 != 135)

DTI_files <- lapply(list.files(path = "D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\", pattern=glob2rx("Streamlines_*.txt"), full.names=T), read.table, header = F, sep = "", dec = ",")
names(DTI_files) <- gsub('.txt', "", list.files(path = "D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\", pattern=glob2rx("Streamlines_*.txt")))

DTI = data.frame(matrix(0, nrow = 140, ncol = 140))

for(x in 1:ncol(DTI)){
	for(j in 1:nrow(DTI)){
		DTI[x,j] = mean(c(DTI_files[[1]][x,j], DTI_files[[2]][x,j], DTI_files[[3]][x,j], DTI_files[[4]][x,j], DTI_files[[5]][x,j])) } }

total$tract = 0  ; total$tract_p = 0

for(m in 1:nrow(total)){
	total$tract[m] = DTI[total$ROI1[m], total$ROI2[m]]
	total$tract_p[m] = try(t.test(c(DTI_files[[1]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[2]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[3]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[4]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[5]][total$ROI1[m], total$ROI2[m]]))$p.value, TRUE)
}




total2 = subset(total, tract >= 50 & model_p_coef < 0)

total2$FWE_t <- p.adjust(total2$model_t, method ="bonferroni", n = length(total2$model_t))
total2$FWE_glm <- p.adjust(total2$model_p, method ="bonferroni", n = length(total2$model_p))

total3 = subset(total2, FWE_t <= 0.05 & FWE_glm <= 0.05)

## OR




## Count how many tract & if > 0

DTI_files <- lapply(list.files(path = "D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\", pattern=glob2rx("Streamlines_*.txt"), full.names=T), read.table, header = F, sep = "", dec = ",")
names(DTI_files) <- gsub('.txt', "", list.files(path = "D:\\Projets\\Projet - FUS\\Macaque Paul\\Matrices\\", pattern=glob2rx("Streamlines_*.txt")))

DTI = data.frame(matrix(0, nrow = 140, ncol = 140))

for(x in 1:ncol(DTI)){
	for(j in 1:nrow(DTI)){
		DTI[x,j] = mean(c(DTI_files[[1]][x,j], DTI_files[[2]][x,j], DTI_files[[3]][x,j], DTI_files[[4]][x,j], DTI_files[[5]][x,j])) } }

total$tract = 0  ; total$tract_p = 0

for(m in 1:nrow(total)){
	total$tract[m] = DTI[total$ROI1[m], total$ROI2[m]]
	total$tract_p[m] = try(t.test(c(DTI_files[[1]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[2]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[3]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[4]][total$ROI1[m], total$ROI2[m]],
		DTI_files[[5]][total$ROI1[m], total$ROI2[m]]))$p.value, TRUE)
}

total2 = subset(total, tract >= 50 & tendency == "Decreasing diff")

total2$FWE_model <- p.adjust(total2$model_p, method ="bonferroni", n = length(total2$model_p))
total2$FWE_max_diff <- p.adjust(total2$max_diff, method ="bonferroni", n = length(total2$max_diff))

total3 = subset(total2, FWE_model <= 0.05)
total4 = subset(total2, FWE_max_diff <= 0.05)


