setwd("/mnt/GEMMA/pca")
# Load necessary libraries
library(ggplot2)
library(data.table)
)

# Read the .sscore file (PCA ran on unrelated samples, and related samples were projected)
score_data <- fread("pca_all_projected.sscore")

colnames(score_data)[1]="FID"

# Select the columns for PCA components (PC1 to PC10)
pca_scores <- score_data[, .(FID, IID, PC1 = PC1_AVG, PC2 = PC2_AVG, PC3 = PC3_AVG, 
                             PC4 = PC4_AVG, PC5 = PC5_AVG, PC6 = PC6_AVG, 
                             PC7 = PC7_AVG, PC8 = PC8_AVG, PC9 = PC9_AVG, PC10 = PC10_AVG)]

pca_scores$country=substr(pca_scores$FID, 1,2)

cols=c("dodgerblue","tomato","gold1")

# ggplot to plot the PCA, change the PC:s and name to plot PC:s 1-14
pdf("pca.pdf",width=6, height=6)
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = country), size=3) +   # Color by family ID or any other variable
  #geom_text(aes(label = IID), check_overlap = TRUE, size = 2, hjust = 0.5, vjust = 0.5) +  # Add sample ID labels to detect outliers
  theme_minimal() +
  scale_color_manual(values=cols) +
  labs(title = "PCA Plot: Projected Related Samples")
dev.off()
