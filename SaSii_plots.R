#define data directory
setwd("./results")


if (!require(tidyr)) install.packages('tidyr')
if (!require(ggplot2)) install.packages('ggplot2')
data_summary <- function(data, varname, groupnames){
  if (!require(plyr)) install.packages('plyr')
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#Change file extension
files <- list.files(pattern="*.sso")
newfiles <- gsub(".sso$", ".txt", files)
file.rename(files, newfiles)


#Fig 1. Fraction of common alleles (freq > 5%) detected
####
#import and adjust data frame
df <- read.table("1-5percent_rate.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(x=n, y=var,group=n)) + 
  geom_boxplot()+
  geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
  scale_y_continuous(name="Fraction of common alleles detected")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("1-5percent_rate.png", width= 6, heigh=4, units="in", dpi=300)

#Fig . Fraction of common alleles (freq > 5%) detected
####
#import and adjust data frame
df <- read.table("2-all_percent_rate.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(x=n, y=var,group=n)) + 
  geom_boxplot()+
  scale_y_continuous(name="Fraction of all alleles detected")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("2-all_percent_rate.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 3. Fraction of common alleles (freq > 5%) detected
####
#import and adjust data frame
df <- read.table("3-rare_percent_rate.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(x=n, y=var,group=n)) + 
  geom_boxplot()+
  scale_y_continuous(name="Fraction of rare alleles detected")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("3-rare_percent_rate.png", width= 6, heigh=4, units="in", dpi=300)


#Fig 4. Impact of sample size on the accuracy of mean sample allele frequency

#import and adjust data frame
df <- read.table("4-freq_dif.txt", skip=3)
colnames(df) <- c("n", "locus", "diff")
df <- separate(data = df, col = locus, into = c("locus", "allele"))

#organize data for the plot
df2 <- data_summary(df, varname="diff", groupnames=c("n", "locus"))

#plot the graphic and save as vector
p <- ggplot(df, aes(x=n, y=diff,group=n)) + 
  geom_boxplot()+
  geom_hline(yintercept=0.01, linetype="dashed", color = "red")+
  scale_y_continuous(name="Mean difference from real allele frequency")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("4-freq_dif.png", width= 6, heigh=4, units="in", dpi=300)


#Fig 5. Impact of sample size on the precision of sample allele frequencies (most common and most rare alleles)
df <- read.table("5-freq_impact.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd", "min", "max")
df <- separate(data = df, col = locus, into = c("freq", "allele"), sep=8)


#plot the graphic and save as vector
p <- ggplot(df, aes(n, var, group=allele, color=allele))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var), fill="white")+
  scale_y_continuous(name="Allele frquency")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("5-freq_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 6. Impact of sample size on the accuracy and precision of expected heterosigosity
#import and adjust data frame
df <- read.table("6-He_impact.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd", "min", "max")

loci <- nlevels(factor(df$locus))
#plot the graphic and save as vector
p <- ggplot(df, aes(n, var, group=locus))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var), fill="white")+
  facet_wrap(~locus, nrow=round(loci/3))+
  scale_y_continuous(name="Expected heterozigosity")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("6-He_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 7. Impact of sample size on the accuracy and precision of multiloci mean expected heterozigosity
df <- read.table("7-meanHe_impact.txt")
colnames(df) <- c("n", "var", "sd", "min", "max")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var),fill="white" )+
  scale_y_continuous(name="Mean expected heterozigosity")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("7-meanHe_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 8. Impact of sample size on the accuracy and precision of observed heterosigosity
#import and adjust data frame
df <- read.table("8-Ho_impact.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd", "min", "max")

loci <- nlevels(factor(df$locus))
#plot the graphic and save as vector
p <- ggplot(df, aes(n, var, group=locus))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var),fill="white")+
  facet_wrap(~locus, nrow=round(loci/3))+
  scale_y_continuous(name="Observed heterozigosity")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("8-Ho_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 9. Imapct of sample size on the accuracy and precision of multiloci mean observed heterozigosity
####
#import and adjust data frame
df <- read.table("9-meanHo_impact.txt")
colnames(df) <- c("n", "var", "sd", "min", "max")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var),fill= "white")+
  scale_y_continuous(name="Mean observed heterozigosity")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("9-meanHo_impact.png", width= 6, heigh=4, units="in", dpi=300)


#Fig 10. Impact of sample size on the accuracy of mean Ho
#import and adjust data frame
df <- read.table("10-Ho_diff.txt", skip=3)
colnames(df) <- c("n", "locus", "diff")

#plot the graphic and save as vector
p <- ggplot(df, aes(x=n, y=diff,group=n)) + 
  geom_boxplot()+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  scale_y_continuous(name="Mean difference from real Ho")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("10-Ho_diff.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 11. Impact of sample size on the accuracy of mean He
#import and adjust data frame
df <- read.table("11-He_diff.txt", skip=3)
colnames(df) <- c("n", "locus", "diff")

#plot the graphic and save as vector
p <- ggplot(df, aes(x=n, y=diff,group=n)) + 
  geom_boxplot()+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  scale_y_continuous(name="Mean difference from real He")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("11-He_diff.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 12. Impact pf sample size on mean Fst between samples and the true population
df <- read.table("12-Fst.txt")
colnames(df) <- c("n", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  geom_errorbar(aes(ymin=var-sd, ymax=var+sd), width=.2) +
  scale_y_continuous(name="Mean pairwise Fst")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("12-Fst.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 13. Impact pf sample size on mean Nei's genetic distance between samples and the true population
df <- read.table("13-Nei.txt")
colnames(df) <- c("n", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_errorbar(aes(ymin=var-sd, ymax=var+sd), width=.2) +
  scale_y_continuous(name="Mean pairwise Nei's genetic distance")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("13-Nei.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 14. Impact pf sample size on mean Roger's genetic distance between samples and the true population
df <- read.table("14-Roger.txt")
colnames(df) <- c("n", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_errorbar(aes(ymin=var-sd, ymax=var+sd), width=.2) +
  scale_y_continuous(name="Mean pairwise Roger's genetic distance")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("14-Roger.png", width= 6, heigh=4, units="in", dpi=300)

#Change file extension
files <- list.files(pattern="*.txt")
newfiles <- gsub(".txt", ".sso", files)
file.rename(files, newfiles)
