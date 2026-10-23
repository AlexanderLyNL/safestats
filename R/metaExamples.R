# Install and load metafor
install.packages("metafor")
library(metafor)

# 1. Prepare sample data (e.g., 5 studies with effect size yi and variance vi)
data <- data.frame(
  study = c("Study 1", "Study 2", "Study 3", "Study 4", "Study 5"),
  yi = c(0.20, 0.50, -0.10, 0.30, 0.40),
  vi = c(0.04, 0.06, 0.05, 0.03, 0.04)
)

# 2. Fit a random-effects meta-analysis model
res <- rma(yi, vi, data = data, method = "REML")
print(res)

# 3. Generate a forest plot
forest(res, slab = data$study)


# Install and load meta
install.packages("meta")
library(meta)

# Example dataset with treatment and control groups
m.cont <- metacont(
  n.e = c(15, 20, 25),
  mean.e = c(12.5, 14.0, 11.2),
  sd.e = c(2.1, 2.5, 1.8),
  n.c = c(15, 20, 25),
  mean.c = c(10.1, 11.5, 10.9),
  sd.c = c(2.0, 2.4, 1.9),
  studlab = c("Trial A", "Trial B", "Trial C"),
  sm = "MD",          # Mean Difference
  random = TRUE       # Random-effects model
)

# Print summary and create forest plot
summary(m.cont)
forest(m.cont)
