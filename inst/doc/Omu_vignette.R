## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE, warning=FALSE, message=FALSE---------------------------------
library(omu)
library(knitr)
load("../data/c57_nos2KO_mouse_countDF.rda")
load("../data/c57_nos2KO_mouse_metadata.rda")
df_trunc <- c57_nos2KO_mouse_countDF[,1:3]
kable(df_trunc[1:4,])

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
kable(c57_nos2KO_mouse_metadata[1:4,])

## ----echo=FALSE,message=FALSE, warning=FALSE----------------------------------
data = c57_nos2KO_mouse_countDF
data = assign_hierarchy(count_data = data, keep_unknowns = TRUE, identifier = "KEGG")
c57_nos2KO_mouse_countDF = c57_nos2KO_mouse_countDF[,1:2]
DF <- assign_hierarchy(count_data = c57_nos2KO_mouse_countDF, keep_unknowns = TRUE, identifier = "KEGG")
kable(DF[1:4,])

## ----echo=FALSE, warning=FALSE------------------------------------------------
library(ggfortify)
load("../data/c57_nos2KO_mouse_countDF.rda")
c57_nos2KO_mouse_countDF_log <- c57_nos2KO_mouse_countDF
c57_nos2KO_mouse_countDF_log <- log(c57_nos2KO_mouse_countDF_log[,3:31])
c57_nos2KO_mouse_countDF_log <- cbind(c57_nos2KO_mouse_countDF[,1:2], c57_nos2KO_mouse_countDF_log)
PCA <- PCA_plot(count_data = c57_nos2KO_mouse_countDF_log, metadata = c57_nos2KO_mouse_metadata, variable = "Grouped", color = "Grouped", response_variable = "Metabolite") + theme_bw() + theme(panel.grid = element_blank())
plot(PCA)

## ----echo=FALSE,message=FALSE, warning=FALSE----------------------------------
DF_stats <- omu_summary(count_data = data, metadata = c57_nos2KO_mouse_metadata, numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment", log_transform = TRUE, p_adjust = "BH", test_type = "welch")
DF_stats_trunc = DF_stats[,c(1,6,7,8,9,10)]
kable(DF_stats_trunc[1:3,])

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
c57_nos2KO_mouse_metadata$Grouped <-factor(paste0(c57_nos2KO_mouse_metadata$Background,
c57_nos2KO_mouse_metadata$Treatment))
kable(c57_nos2KO_mouse_metadata[1:4,])

## ----echo=FALSE,message=FALSE, warning=FALSE----------------------------------
DF_stats_grouped <- omu_summary(count_data = data, metadata = c57_nos2KO_mouse_metadata, numerator = "WTStrep", denominator = "WTMock", response_variable = "Metabolite", Factor = "Grouped", log_transform = TRUE, p_adjust = "BH", test_type = "welch")
DF_stats_grouped_trunc = DF_stats_grouped[,c(1,6,7,8,9,10)]
kable(DF_stats_grouped_trunc[1:3,])

## ----echo=FALSE,message=FALSE, warning=FALSE----------------------------------
DF_stats_sub <- subset(DF_stats, Class=="Organic acids")
DF_stats_sub <- DF_stats_sub[which(DF_stats_sub[,"padj"] <= 0.05),]
DF_s_trunc <- DF_stats_sub[c(1:5),]
DF_s_trunc <- DF_s_trunc[,c(1,9,12,13,14,43)]
kable(DF_s_trunc)

## ----echo=FALSE,message=FALSE, warning=FALSE,cache = FALSE, results = 'hide'----
DF_s_trunc_g <- KEGG_gather(DF_stats_sub)

kable(head(DF_s_trunc_g, n = 5))

## ----echo=FALSE, warning=FALSE------------------------------------------------
DF_stats_counts <- count_fold_changes(count_data = DF_stats, column = "Class", sig_threshold = 0.05, keep_unknowns = FALSE)
kable(DF_stats_counts)

## ----echo=FALSE, fig.keep='all', results='hide'-------------------------------
library(ggplot2)
Class_Bar_Plot <- plot_bar(fc_data = DF_stats_counts, fill = c("dodgerblue2", "firebrick2"), outline_color = c("black", "black"), size = c(1,1)) + labs(x = "Class") + theme(panel.grid = element_blank())
plot(Class_Bar_Plot)

## ----echo=FALSE---------------------------------------------------------------
DF_ra <- ra_table(fc_data = DF_stats_counts, variable = "Class")
kable(DF_ra)

## ----echo=FALSE, warning=FALSE------------------------------------------------
Pie_Chart <- pie_chart(ratio_data = DF_ra, variable = "Class", column = "Decrease", color = "black")
plot(Pie_Chart)

## ----echo=FALSE, warning=FALSE, results='hide', fig.keep='all'----------------
Volcano_Plot <- plot_volcano(count_data = DF_stats_grouped, size = 2, column = "Class", strpattern = c("Organic acids", "Carbohydrates"), fill = c("firebrick2","white","dodgerblue2"), color = c("black", "black", "black"), alpha = c(1,1,1), shape = c(21,21,21)) + theme_bw() + theme(panel.grid = element_blank())
plot(Volcano_Plot)

## ----echo=FALSE, warning=FALSE------------------------------------------------
library(ggfortify)
load("../data/c57_nos2KO_mouse_countDF.rda")
c57_nos2KO_mouse_countDF_log <- c57_nos2KO_mouse_countDF
c57_nos2KO_mouse_countDF_log <- log(c57_nos2KO_mouse_countDF_log[,3:31])
c57_nos2KO_mouse_countDF_log <- cbind(c57_nos2KO_mouse_countDF[,1:2], c57_nos2KO_mouse_countDF_log)
PCA <- PCA_plot(count_data = c57_nos2KO_mouse_countDF_log, metadata = c57_nos2KO_mouse_metadata, variable = "Grouped", color = "Grouped", response_variable = "Metabolite") + theme_bw() + theme(panel.grid = element_blank())
plot(PCA)

## ----echo=FALSE, warning=FALSE------------------------------------------------
DF <- assign_hierarchy(count_data = c57_nos2KO_mouse_countDF, keep_unknowns = TRUE, identifier = "KEGG")
heatmap <- plot_heatmap(count_data = DF, metadata = c57_nos2KO_mouse_metadata, Factor = "Treatment", response_variable = "Metabolite", log_transform = TRUE, high_color = "goldenrod2", low_color = "midnightblue", aggregate_by = "Class") + theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size = 6), axis.text.y = element_text(size = 6))

plot(heatmap)

## ----echo=FALSE, warning=FALSE------------------------------------------------
DF_c <- subset(DF, Class == "Carbohydrates")
heatmap <- plot_heatmap(count_data = DF_c, metadata = c57_nos2KO_mouse_metadata, Factor = "Treatment", response_variable = "Metabolite", log_transform = TRUE, high_color = "goldenrod2", low_color = "midnightblue") + theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size = 6), axis.text.y = element_text(size = 6))
plot(heatmap)

## ----echo=FALSE, warning=FALSE------------------------------------------------
DF_carbs <- subset(DF, Class == "Carbohydrates")
heatmap <- plot_heatmap(count_data = DF_carbs, metadata = c57_nos2KO_mouse_metadata, Factor = "Treatment", response_variable = "Metabolite", log_transform = TRUE, high_color = "goldenrod2", low_color = "midnightblue", aggregate_by = "Subclass_2") + theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size = 6), axis.text.y = element_text(size = 6))
plot(heatmap)

## ----echo=FALSE, warning=FALSE------------------------------------------------
DF_carbs_trunc <- DF_carbs[1:10,]
boxplot_carbs <- plot_boxplot(count_data = DF_carbs_trunc, metadata = c57_nos2KO_mouse_metadata, log_transform = TRUE, Factor = "Treatment", response_variable = "Metabolite", fill_list = c("darkgoldenrod1", "dodgerblue2"))
plot(boxplot_carbs)

## ----echo=FALSE, warning=FALSE------------------------------------------------
boxplot_carbs_sc2 <- plot_boxplot(count_data = DF_carbs, metadata = c57_nos2KO_mouse_metadata, log_transform = TRUE, Factor = "Treatment", response_variable = "Metabolite", fill_list = c("darkgoldenrod1", "dodgerblue2"), aggregate_by = "Subclass_2")
plot(boxplot_carbs_sc2)

## ----echo=FALSE, warning=FALSE------------------------------------------------
boxplot_class <- plot_boxplot(count_data = DF, metadata = c57_nos2KO_mouse_metadata, log_transform = TRUE, Factor = "Treatment", response_variable = "Metabolite", fill_list = c("darkgoldenrod1", "dodgerblue2"), aggregate_by = "Class")
plot(boxplot_class)

