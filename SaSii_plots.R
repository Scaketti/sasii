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
# Q? O que é ML?
####
#import and adjust data frame
df <- read.table("1-5percent_rate.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var, group=locus))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_errorbar(aes(ymin=var-sd, ymax=var+sd), width=.2) +
  scale_y_continuous(name="Fraction of common alleles detected")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

p <- ggplot(df, aes(x=n, y=var,group=n)) + 
  geom_boxplot()+
  scale_y_continuous(name="Fraction of common alleles detected")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("1-5percent_rate.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 2. Impact of sample size on the accuracy of mean sample allele frequency

#import and adjust data frame
df <- read.table("2-freq_dif.txt", skip=3)
colnames(df) <- c("n", "locus", "diff")
df <- separate(data = df, col = locus, into = c("locus", "allele"))

#organize data for the plot
df2 <- data_summary(df, varname="diff", groupnames=c("n", "locus"))

#plot the graphic and save as vector
p <- ggplot(df2, aes(n, diff, group=locus))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_errorbar(aes(ymin=diff-sd, ymax=diff+sd), width=.2) +
  scale_y_continuous(name="Mean difference from real allele frequency")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

p <- ggplot(df, aes(x=n, y=diff,group=n)) + 
  geom_boxplot()+
  scale_y_continuous(name="Mean difference from real allele frequency")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("2-freq_dif.png", width= 6, heigh=4, units="in", dpi=300)


#Fig 3. Impact of sample size on the precision of sample allele frequencies (most common and most rare alleles)
df <- read.table("3-freq_impact.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd", "min", "max")
df <- separate(data = df, col = locus, into = c("freq", "allele"), sep=8)


#plot the graphic and save as vector
p <- ggplot(df, aes(n, var, group=allele, color=allele))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var), fill="white")+
  scale_y_continuous(name="Allele frquency")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("3-freq_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 4. Impact of sample size on the accuracy and precision of expected heterosigosity
#import and adjust data frame
df <- read.table("4-He_impact.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd", "min", "max")

loci <- nlevels(factor(df$locus))
#plot the graphic and save as vector
p <- ggplot(df, aes(n, var, group=locus))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var), fill="white")+
  facet_wrap(~locus, nrow=loci/3)+
  scale_y_continuous(name="Expected heterozigosity")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("4-He_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 5. Impact of sample size on the accuracy and precision of multiloci mean expected heterozigosity
df <- read.table("5-meanHe_impact.txt")
colnames(df) <- c("n", "var", "sd", "min", "max")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var),fill="white" )+
  scale_y_continuous(name="Mean expected heterozigosity")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("5-meanHe_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 6. Impact pf sample size on mean Fst between samples and the true population
df <- read.table("6-Fst.txt")
colnames(df) <- c("n", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_errorbar(aes(ymin=var-sd, ymax=var+sd), width=.2) +
  scale_y_continuous(name="Mean pairwise Fst")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("6-Fst.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 6. Impact pf sample size on mean Nei's genetic distance between samples and the true population
df <- read.table("7-Nei.txt")
colnames(df) <- c("n", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_errorbar(aes(ymin=var-sd, ymax=var+sd), width=.2) +
  scale_y_continuous(name="Mean pairwise Nei's genetic distance")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("7-Nei.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 8. Impact pf sample size on mean Roger's genetic distance between samples and the true population
df <- read.table("8-Roger.txt")
colnames(df) <- c("n", "var", "sd")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, var))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  geom_errorbar(aes(ymin=var-sd, ymax=var+sd), width=.2) +
  scale_y_continuous(name="Mean pairwise Roger's genetic distance")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("8-Roger.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 9. Imapct of sample size on the accuracy and precision of multiloci mean observed heterozigosity
####
# Trocar He por Ho, rever col names (sd, min, max)
#####
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

ggsave("9-Ho_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 10. Impact of sample size on the accuracy and precision of observed heterosigosity
#import and adjust data frame
df <- read.table("10-Ho_impact.txt", skip=3)
colnames(df) <- c("n", "locus", "var", "sd", "min", "max")

loci <- nlevels(factor(df$locus))
#plot the graphic and save as vector
p <- ggplot(df, aes(n, var, group=locus))+
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_crossbar(aes(ymin = var-sd, ymax = var+sd, x = n, y = var),fill="white")+
  facet_wrap(~locus, nrow=loci/3)+
  scale_y_continuous(name="Observed heterozigosity")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("10-Ho_impact.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 11. Impact of sample size on the accuracy of mean Ho
#import and adjust data frame
df <- read.table("11-Ho_diff.txt", skip=3)
colnames(df) <- c("n", "locus", "diff")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, diff, group=locus))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  scale_y_continuous(name="Mean allele frequency")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

p <- ggplot(df, aes(x=n, y=diff,group=n)) + 
  geom_boxplot()+
  scale_y_continuous(name="Mean difference from real Ho")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("11-Ho_diff.png", width= 6, heigh=4, units="in", dpi=300)

#Fig 12. Impact of sample size on the accuracy of mean He
#import and adjust data frame
df <- read.table("12-He_diff.txt", skip=3)
colnames(df) <- c("n", "locus", "diff")

#plot the graphic and save as vector
p <- ggplot(df, aes(n, diff, group=locus))+
  geom_point(stat="identity")+
  geom_line(stat="identity")+  
  scale_y_continuous(name="Mean difference from real He")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

p <- ggplot(df, aes(x=n, y=diff,group=n)) + 
  geom_boxplot()+
  scale_y_continuous(name="Mean difference from real He")+
  scale_x_continuous(breaks = seq(0, tail(df$n, n=1), by = 10))+
  theme_classic()

ggsave("12-He_diff.png", width= 6, heigh=4, units="in", dpi=300)

#Change file extension
files <- list.files(pattern="*.txt")
newfiles <- gsub(".txt", ".sso", files)
file.rename(files, newfiles)