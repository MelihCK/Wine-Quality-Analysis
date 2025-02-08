#Melih Can KANMAZ
#Mustafa Ugur YALCIN
#Tevfik OGUZ


library(ggplot2)
library(tidyr)
library(GGally)
library(corrplot)
library(dplyr)

# 1 - Explanatory Data Analysis

wine = read.csv("Wine_Quality_Data.csv")

head(wine)
tail(wine)

summary(wine)

str(wine)

colSums(is.na(wine)) # No NA values

#Selecting continuous variables to visualize appropriately, if it is not numeric we use boxplot
cont_var = names(wine)[sapply(wine, is.numeric)]

par(mfrow = c(3, 4)) 

#Boxplot to see the relationship of variables with color
for (var in cont_var) {
  boxplot(wine[[var]] ~ wine$color, 
          main = paste("Color vs", var), 
          col = "orange",
          ylab = var,
          xlab="color")
}

#Plot for relationship of continuous var with quality
for (var in cont_var) {
  plot(wine[[var]], wine$quality, 
       main = paste(var, "vs Quality"), 
       xlab = var, ylab = "Quality", col = "blue", pch = 19)
}

#Plot that includes all the continuous variables' relations with each other
pairplot = ggpairs(wine[,cont_var])

#to save to the pc
ggsave("EDA_pair_plot.png", plot = pairplot, width = 40, height = 25, units = "cm", dpi = 300)

#Correlation plot
par(mfrow = c(1, 1))
corrplot(cor(wine[,cont_var]),method = "color", addCoef.col = "black",)

#Outlier removing according to the quartiles -1.5Q1<x<1.5Q3
clean_outliers = function(df) {
  for (col in cont_var) {
      IQR = IQR(df[[col]], na.rm = TRUE)
      
      lw = quantile(df[[col]], 0.25, na.rm = TRUE) - 1.5 * IQR
      up = quantile(df[[col]], 0.75, na.rm = TRUE) + 1.5 * IQR
      
      df = df[df[[col]] >= lw & df[[col]] <= up, ]
  }
  return(df)
}

wine_cleaned = clean_outliers(wine)
wine_cleaned = clean_outliers(wine_cleaned) # needed to run couple of times to get rid of all the outliers

#After cleaning, checking the outliers in the plot
par(mfrow = c(3, 4))
for (var in cont_var) {
  boxplot(wine_cleaned[[var]] ~ wine_cleaned$color, 
          main = paste("Color vs", var), 
          col = "orange",
          ylab = var,
          xlab="color")
}


# 2- Inference of Mean Vector
#necessary libraries 
library("ICSNP")
library(MVN)
library(bestNormalize) #Pelin Hoca's suggestion

#bestNormalize gave the best transformation as orderNorm
bn_ph = bestNormalize(wine_cleaned$pH)
bn_quality = bestNormalize(wine_cleaned$quality)

wine_cleaned$bn_ph = predict(bn_ph)
wine_cleaned$bn_quality = predict(bn_quality)

par(mfrow = c(1, 1))
plot(bn_ph, leg_loc = "bottomright")

#to have a close prediction on mu0 vector
mean_vector = colMeans(wine_cleaned[,cont_var])
mean_vector

y = wine_cleaned[,c("bn_ph","bn_quality")]
dim(y)

#Still after transformation the data was not normal
result=mvn(y,mvnTest = "mardia")
result$multivariateNormality

result=mvn(y,mvnTest = "mardia")
result$univariateNormality

#We tried many transformations but none of the make the data normal so after asking Vilda Hoca
#we continued as assuming the data is normal for the sake of the practice

#library(MASS)
#boxcox_result = boxcox(lm(alcohol ~ 1, data = wine_cleaned))  
#lambda = boxcox_result$x[which.max(boxcox_result$y)] 
#clean_data$bc_alcohol = (clean_data$alcohol^lambda - 1) / lambda  

#boxcox_result = boxcox(lm(density ~ 1, data = wine_cleaned))  
#lambda = boxcox_result$x[which.max(boxcox_result$y)]  
#clean_data$bc_density = (clean_data$density^lambda - 1) / lambda  

#y_clean = y[y > 0, ]  
#log_y = log(y_clean)  

#test = mvn(log_y, mvnTest = "mardia")
#test$multivariateNormality

#After trying multiple ways to make the response variable normal for Hotelling T test, the data was still not normal. For the sake of the practice we continue with non_normal data.

#violin graph
library (psych)
error.bars (y, ylab="Group Means", xlab=" Dependent Variables")

#mu0 
bn_mu0 = colMeans(wine_cleaned[,c("bn_ph","bn_quality")])
bn_mu0
dim(bn_mu0)
#one-at-a-time confidence interval for each variable 
ph = lm(bn_ph ~ 1, data = wine_cleaned)
confint(ph)

quality = lm(bn_quality ~ 1, data = wine_cleaned)
confint(quality)

mu0 = c(1,1)
#Hotelling T test
HotellingsT2(y,mu=mu0)

# Reject H0 and conclude that mean of the ph and quality differ from 1.

library(mvdalab)
MVcis(y)


#3.COMPARISONS OF SEVERAL MULTIVARIATE MEANS
subset_wine = wine_cleaned %>% select(pH, quality, color) %>% mutate(log_pH = log(pH), log_quality = log(quality))
subset_wine %>% head()

subset_wine %>% group_by(color) %>%  summarise(n = n(), 
                                               mean_quality = mean(quality), 
                                               sd_quality = sd(quality),
                                               mean_ph = mean(pH),
                                               sd_ph = sd(pH))

library(gridExtra)
p1 = ggplot(subset_wine, aes(x = color, y = quality, fill = color)) + geom_boxplot(outlier.shape = NA) + theme(legend.position="top")+
  labs(title = "Boxplot of quality by color" )
p2 = ggplot(subset_wine, aes(x = color, y = pH, fill = color)) + geom_boxplot(outlier.shape = NA) + theme(legend.position="top")+
  labs(title = "Boxplot of pH by color")
grid.arrange(p1, p2, ncol=2)

library(rstatix)
subset_wine %>% group_by(color) %>%  shapiro_test(log_pH,log_quality)


library(heplots)
boxM(Y = cbind(subset_wine$log_pH,subset_wine$log_quality), group = factor(subset_wine$color))

#since data is not normal we used Levene test
library(car)
leveneTest(log_pH ~ color, data = subset_wine)

m1 = manova(cbind(log_pH,log_quality) ~ color, data = subset_wine)
summary(m1)

summary.aov(m1)



#4.PRINCIPAL COMPONENTS ANALYSIS and PRINCIPAL COMPONENTS REGRESSION

# Standardization for the PCA as first step
scaled_wine = scale(wine_cleaned[, cont_var])

#correlation and covariance matrix as a second step
cor(scaled_wine)

cov(scaled_wine)
#both cor and cov tables looks same as expected 

# Exclude the 'quality' since it is response variable
scaled_wine1=scaled_wine[, colnames(scaled_wine) != "quality"]

# Perform PCA
pca1 = prcomp(scaled_wine1)
summary(pca1)

# PCA features
pca1$rotation
pca1$x
pca1$sdev


library(factoextra)

#Scree plot to represent the proportion of variance explained by each PC
fviz_eig(pca1,addlabels=TRUE) 

#Since the first 7 PC's are better predictors than others we continue with them
pca=pca1$x[,1:7]

head(pca)

# Correlation analysis
res1 = cor(pca, method="pearson")
corrplot::corrplot(res1, method= "color", order = "hclust")

cor(scaled_wine1,pca)


fviz_pca_var(pca1, axes = c(1, 7)) +
  ggtitle("PCA 1 and PCA 7")

fviz_pca_var(pca1,axes = c(3, 5)) +
  ggtitle("PCA 3 and PCA 5")


# Highlight variable contributions with color gradients blue, red, green
fviz_pca_var(pca1, col.var = "contrib")+ scale_color_gradient2( low="red", mid="green",
                                                                high="blue", midpoint=96, space = "Lab")
#top 5 contributing variables to the PCs
fviz_pca_var(pca1, select.var = list(contrib =5))

#since data has numereous instances the following two plots are not that helpful
fviz_pca_ind(pca1, col.ind = "#FFA721")

fviz_contrib(pca1, choice = "ind", axes = 1:2) + coord_flip()

# PCA with grouping by 'pH' and adding confidence ellipses
fviz_pca_ind(pca1, label="none", habillage=wine_cleaned$pH,
             addEllipses=TRUE, ellipse.level=0.95)


# Principal Components Regression (PCR)
ols.data = data.frame(quality=scaled_wine[,12],pca)

lmodel = lm(quality ~ ., data = ols.data)
summary(lmodel)

#mean((ols.data$quality - predict(lmodel))^2) 

library(psych)
library(party)

num_data <- char2numeric(wine) #Categorical to Numerical

# FACTOR ANALYSIS AND FACTOR ROTATION
# We prefer Varimax and Quartimax as a Rotation Method since the estimated weights for the factor scores are probably unmatched as the other methods.

pa <- fa(r = num_data,
         nfactors = 3,
         rotate = "quartimax",
         fm="pa",
         residuals=TRUE) #Principal Axis Factor Analysis

ml <- fa(r = num_data,
         nfactors = 3,
         rotate = "varimax",
         fm="ml",
         residuals = TRUE) #Maximum Likelihood

gls <- fa(r = num_data,
         nfactors = 3,
         rotate = "varimax",
         fm="gls",
         residuals = TRUE) #Generalized Least Square Method

ols <- fa(r = num_data,
          nfactors = 3,
          rotate = "varimax",
          fm="ols",
          residuals = TRUE) #Ordinary Least Square Method


#Since the Cumulative Variance and BIC values for ML method is larger than OLS and GLS Method, We will work on the ML Method for the analysis.
#As we can see on the output of ML method;

# ML1 is related to Fixed Acidity, Citric Acid, pH.
# We can correspond this factor loading as a Tendency to Acidity.

# ML2 is related to Residual Sugar, Density, Alcohol, Quality.
# We can correspond this factor loading as a Sweetness Profile.

# ML3 is related to Volatile Acidity, Chlorides, Free Sulfur Dioxide, Total Sulfur Dioxide, Sulphates, Color.
# We can correspond this factor loading as a Wine Stability.

#What is the best and worst feature that is explained in factor analysis? Since h^2 value of color of wine is equal to 0.93, it says 93% of the total variance in the variable is explained by the factors. As for worst,quality variable has a  0.89 value of u^2, this tells otherwise.


# DISCRIMINATION AND CLASSIFICATION

# We want the graphical presentation of all these variables and their relation with each other, so,

pairs.panels(num_data)
set.seed(467)
ind <- sample(2, nrow(num_data),
              replace = TRUE, prob = c(0.7,0.3)) #70% for training set, and the rest is for validation set.

train.data <- num_data[ind==1,]
test.data <- num_data[ind==2,]
myf <- data$quality~alcohol+pH+density+residual_sugar+citric_acid
wine_ctree <- ctree(myf,data = train.data)

table(predict(wine_ctree),train.data$quality)

# In the table, as for row, we see the actual data and as for column, we see the predicted values in the tree.
# For example, in the first row the actual value is 5.016, there are 59 items and only 40 of the predicted items are accurately identified.

plot(wine_ctree)

# The ctree function shows us the conditional inferences about the variables. As we can see the conditional tree, it can interpreted as you like.

# Clustering
# Since the dataset is too large for a readable dendrogram,
# we'll sample a subset of the data for visualization
# Used examples from ODTUCLASS
sample_size = 50
sample_indices = sample(1:nrow(wine_cleaned), sample_size)
dm_sample = dist(wine_cleaned[sample_indices, cont_var])

jpeg("single_linkage_dendrogram_sample.jpg", width = 5000, height = 4000, res = 600)
plot(hclust(dm_sample, method = "single"),
     main = "Single Linkage Dendrogram",
     xlab = "Sample Index",
     ylab = "Distance",
     sub = paste("Random sample of", sample_size, "observations")) 
dev.off()                                        

jpeg("complete_linkage_dendrogram.jpg", width = 5000, height = 4000, res = 600)
plot(hclust(dm_sample, method = "complete"),
     main = "Complete Linkage Dendrogram",
     xlab = "Sample Index",
     ylab = "Distance",
     sub = paste("Random sample of", sample_size, "observations"))
dev.off()

jpeg("average_linkage_dendrogram.jpg", width = 5000, height = 4000, res = 600)
plot(hclust(dm_sample, method = "average"),
     main = "Average Linkage Dendrogram",
     xlab = "Sample Index", 
     ylab = "Distance",
     sub = paste("Random sample of", sample_size, "observations"))
dev.off()

jpeg("ward_linkage_dendrogram.jpg", width = 5000, height = 4000, res = 600)
plot(hclust(dm_sample, method = "ward.D2"),
     main = "Ward Method Dendrogram",
     xlab = "Sample Index",
     ylab = "Distance",
     sub = paste("Random sample of", sample_size, "observations"))
dev.off()

# K-means Clustering

pca_data = data.frame(PC1=pca1$x[,1], PC2=pca1$x[,2])
wss = sapply(1:20, function(k) {
  sum(kmeans(pca_data, centers=k)$withinss)
})

#Elbow Plot
jpeg("kmeans_elbow_plot.jpg", width = 3000, height = 2000, res = 300)
plot(1:20, wss,
     type="b", 
     pch = 19,
     frame = FALSE,
     xlab="Number of Clusters",
     ylab="Within-Cluster Sum of Squares",
     main="Elbow Plot for Optimal k")
dev.off()

# Silhouette
library(cluster)
library(factoextra)

sil_score = sapply(2:20, function(k) {
  temp_clusters = kmeans(pca_data, centers=k)$cluster
  mean(silhouette(temp_clusters, 
                  dist(pca_data))[,3])
})

sil_df = data.frame(
  k = 2:20,
  score = sil_score
)

# silhouette scores
jpeg("kmeans_silhouette_plot.jpg", width = 3000, height = 2000, res = 300)
ggplot(sil_df, aes(x = k, y = score)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Silhouette Score for different cluster amount",
       x = "Number of clusters",
       y = "Silhouette Score") +
  scale_x_continuous(breaks = 2:20) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
dev.off()

# Kmeans
km = kmeans(pca_data, centers=3)
km$centers
km$cluster

# Cluster Scatter plot
cluster_viz_data = data.frame(
  PC1 = pca1$x[,1],
  PC2 = pca1$x[,2],
  Cluster = as.factor(km$cluster)
)

jpeg("kmeans_cluster_visualization3.jpg", width = 3000, height = 2000, res = 300)
ggplot(cluster_viz_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.5) +
  stat_ellipse(aes(fill = Cluster), geom = "polygon", alpha = 0.2) +
  theme_minimal() +
  labs(title = "K-means Clusters Visualization",
       x = "First Principal Component", 
       y = "Second Principal Component") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")
dev.off()


#CCA
# Split variables into two groups
# Group 1: Chemical properties
# Group 2: Sensory properties
chemical_vars = wine_cleaned[, c("fixed_acidity", "volatile_acidity", "citric_acid", 
                               "residual_sugar", "chlorides", "free_sulfur_dioxide", 
                               "total_sulfur_dioxide", "density", "pH", "sulphates")]

sensory_vars = wine_cleaned[, c("alcohol", "quality")]

library(CCA)
cca_result = cancor(chemical_vars, sensory_vars)
print("Canonical Correlations:")
print(cca_result$cor)

# Bar Plot of the correlations
canonical_correlations = data.frame(
  Dimension = 1:length(cca_result$cor),
  Correlation = cca_result$cor
)

jpeg("canonical_correlations.jpg", width = 3000, height = 2000, res = 300)
ggplot(canonical_correlations, aes(x = factor(Dimension), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Canonical Correlations by Dimension",
       x = "Canonical Dimension",
       y = "Correlation") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
dev.off()

# Calculating the loadings
library(CCP)
loadings = comput(chemical_vars, sensory_vars, cca_result)
loadings

# Wilks Test
n = nrow(chemical_vars)
p = ncol(chemical_vars)
q = ncol(sensory_vars)
p.asym(cca_result$cor, n, p, q, tstat = "Wilks")

