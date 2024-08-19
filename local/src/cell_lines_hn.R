hnc <- read.xlsx("/mnt/cold1/snaketree/prj/hn/local/share/data/linee_hnc_eugy.xlsx", sheet = 1)
colnames(hnc) <- c("drug", 0, 6, 7, 8, 9, 10, 11)
hnc[c(1:3),"drug"] <- "BMS-833923 + Cetuximab 20ng/µL"
hnc[c(4:6),"drug"] <- "BMS-833923"
hnc <- hnc %>%
  pivot_longer(cols = -drug, names_to = "variable", values_to = "value")

summary_stats <- hnc %>%
  group_by(drug, variable) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value)
  )

#summary_stats$variable <- gsub("μ.M", "", summary_stats$variable)
summary_stats$variable <- as.numeric(summary_stats$variable)
summary_stats$drug <- as.factor(summary_stats$drug)
summary_stats <- summary_stats[order(summary_stats$variable),]

custom_breaks <- unique(summary_stats$variable)

ggplot(summary_stats, aes(x = variable, y = mean_value, color = drug)) +
  geom_point(size = 3) + 
  geom_line() +  
  theme_minimal()+  
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2)+
  scale_x_continuous(breaks = custom_breaks)+scale_color_manual(values = c("forestgreen", "magenta4"))+
  xlab("BMS-833923 µM")+ylab("Cell number")+ggtitle("HSC-2 Cell Lines")

hnc <- read.xlsx("/mnt/cold1/snaketree/prj/hn/local/share/data/linee_hnc_eugy.xlsx", sheet = 2)
colnames(hnc) <- c("drug", 0, 6, 7, 8, 9, 10, 11)
hnc[c(1:3),"drug"] <- "BMS-833923 + Cetuximab 20ng/µL"
hnc[c(4:6),"drug"] <- "BMS-833923"
hnc <- hnc %>%
  pivot_longer(cols = -drug, names_to = "variable", values_to = "value")

summary_stats <- hnc %>%
  group_by(drug, variable) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value)
  )

#summary_stats$variable <- gsub("μ.M", "", summary_stats$variable)
summary_stats$variable <- as.numeric(summary_stats$variable)
summary_stats$drug <- as.factor(summary_stats$drug)
summary_stats <- summary_stats[order(summary_stats$variable),]

custom_breaks <- unique(summary_stats$variable)

ggplot(summary_stats, aes(x = variable, y = mean_value, color = drug)) +
  geom_point(size = 3) + 
  geom_line() +  
  theme_minimal()+  
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2)+
  scale_x_continuous(breaks = custom_breaks)+scale_color_manual(values = c("forestgreen", "magenta4"))+
  xlab("BMS-833923 µM")+ylab("Cell number")+ggtitle("CAL-27 Cell Lines")


hnc <- read.xlsx("/mnt/cold1/snaketree/prj/hn/local/share/data/hnc2.xlsx", sheet = 1)
colnames(hnc) <- c("drug", 0, 6, 7, 8, 9, 10, 11)
hnc[c(4:6),"drug"] <- "BMS-833923 + Cetuximab 20ng/µL"
hnc[c(1:3),"drug"] <- "BMS-833923"
hnc <- hnc %>%
  pivot_longer(cols = -drug, names_to = "variable", values_to = "value")

summary_stats <- hnc %>%
  group_by(drug, variable) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value)
  )

#summary_stats$variable <- gsub("μ.M", "", summary_stats$variable)
summary_stats$variable <- as.numeric(summary_stats$variable)
summary_stats$drug <- as.factor(summary_stats$drug)
summary_stats <- summary_stats[order(summary_stats$variable),]

custom_breaks <- unique(summary_stats$variable)

ggplot(summary_stats, aes(x = variable, y = mean_value, color = drug)) +
  geom_point(size = 3) + 
  geom_line() +  
  theme_minimal()+  
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2)+
  scale_x_continuous(breaks = custom_breaks)+scale_color_manual(values = c("forestgreen", "magenta4"))+
  xlab("BMS-833923 µM")+ylab("Cell number")+ggtitle("HNC0246")

hnc <- read.xlsx("/mnt/cold1/snaketree/prj/hn/local/share/data/hnc2.xlsx", sheet = 2)
colnames(hnc) <- c("drug", 0, 6, 7, 8, 9, 10, 11)
hnc[c(4:6),"drug"] <- "BMS-833923 + Cetuximab 20ng/µL"
hnc[c(1:3),"drug"] <- "BMS-833923"
hnc <- hnc %>%
  pivot_longer(cols = -drug, names_to = "variable", values_to = "value")

summary_stats <- hnc %>%
  group_by(drug, variable) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value)
  )

#summary_stats$variable <- gsub("μ.M", "", summary_stats$variable)
summary_stats$variable <- as.numeric(summary_stats$variable)
summary_stats$drug <- as.factor(summary_stats$drug)
summary_stats <- summary_stats[order(summary_stats$variable),]

custom_breaks <- unique(summary_stats$variable)

ggplot(summary_stats, aes(x = variable, y = mean_value, color = drug)) +
  geom_point(size = 3) + 
  geom_line() +  
  theme_minimal()+  
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2)+
  scale_x_continuous(breaks = custom_breaks)+scale_color_manual(values = c("forestgreen", "magenta4"))+
  xlab("BMS-833923 µM")+ylab("Cell number")+ggtitle("HNC0032")

hnc <- read.xlsx("/mnt/cold1/snaketree/prj/hn/local/share/data/hnc2.xlsx", sheet = 3)
colnames(hnc) <- c("drug", 0, 6, 7, 8, 9, 10, 11)
hnc[c(4:6),"drug"] <- "BMS-833923 + Cetuximab 20ng/µL"
hnc[c(1:3),"drug"] <- "BMS-833923"
hnc <- hnc %>%
  pivot_longer(cols = -drug, names_to = "variable", values_to = "value")

summary_stats <- hnc %>%
  group_by(drug, variable) %>%
  summarise(
    mean_value = mean(value),
    sd_value = sd(value)
  )

#summary_stats$variable <- gsub("μ.M", "", summary_stats$variable)
summary_stats$variable <- as.numeric(summary_stats$variable)
summary_stats$drug <- as.factor(summary_stats$drug)
summary_stats <- summary_stats[order(summary_stats$variable),]

custom_breaks <- unique(summary_stats$variable)

ggplot(summary_stats, aes(x = variable, y = mean_value, color = drug)) +
  geom_point(size = 3) + 
  geom_line() +  
  theme_minimal()+  
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2)+
  scale_x_continuous(breaks = custom_breaks)+scale_color_manual(values = c("forestgreen", "magenta4"))+
  xlab("BMS-833923 µM")+ylab("Cell number")+ggtitle("HNC0044")


