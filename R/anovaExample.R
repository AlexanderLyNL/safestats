# 1. Create a dummy dataset
plantData <- data.frame(
  treatment = factor(rep(c("control", "fertiliserA", "fertiliserB"), each = 10)),
  height = c(
    c(12, 14, 11, 13, 12, 15, 13, 11, 14, 12), # Control group
    c(18, 20, 19, 22, 21, 20, 17, 19, 21, 23), # Fertiliser A group
    c(15, 17, 16, 15, 18, 14, 16, 17, 15, 16)  # Fertiliser B group
  )
)

# 2. View a quick summary of the data
aggregate(height ~ treatment, data = plantData, mean)

# 3. Fit the one-way ANOVA model
anovaModel <- aov(height ~ treatment, data = plantData)


# 4. View the ANOVA summary table
summary(anovaModel)

# 5. Run the post-hoc Tukey HSD test (since the p-value will be significant)
TukeyHSD(anovaModel)
