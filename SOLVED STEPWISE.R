# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(broom)  # For tidy summaries
library(car)    # For VIF calculation
library(MASS)   # For stepAIC

# Set working directory
setwd("C:/Users/rober/OneDrive/Documents/Research MDC")

# Read data
SENIC <- read_csv("SENIC.csv")

# Display column names to verify the presence of 'id'
cat("Columns before removal:\n")
print(names(SENIC))

# Remove the 'id' variable using base R
columns_to_remove <- c("id")
SENIC <- SENIC[, !(names(SENIC) %in% columns_to_remove)]

# Convert region and medSchool to factors
SENIC$region <- as.factor(SENIC$region)
SENIC$medSchool <- as.factor(SENIC$medSchool)

# Display column names to verify removal
cat("Columns after removal:\n")
print(names(SENIC))

# Define the response variable
response_var <- "infectionRisk"

# Create formula using all other variables as predictors
predictors <- setdiff(names(SENIC), response_var)
formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = "+")))

# Fit the initial linear model (full model)
multiplelinear_model <- lm(formula, data = SENIC)

# Perform stepwise regression
# The Akaike information criterion (AIC)
stepwise_model <- stepAIC(multiplelinear_model, direction = "both", trace = 0) # Both means forward selection and backward elimination

# Get model summary
model_summary <- summary(stepwise_model)

# Create tidy summary of coefficients with confidence intervals
tidy_model <- tidy(stepwise_model, conf.int = TRUE)

# Prepare formatted regression formula
intercept <- round(coef(stepwise_model)[1], 4)
terms <- paste(round(coef(stepwise_model)[-1], 4), names(coef(stepwise_model)[-1]), sep = "*", collapse = " + ")
regression_formula <- paste(response_var, "=", intercept, "+", terms)

# Calculate VIF for each predictor
vif_values <- vif(stepwise_model)

# Debugging output to check VIF values
cat("VIF values:\n")
print(vif_values)
cat("Length of vif_values:\n")
print(length(vif_values))
cat("Names of vif_values:\n")
print(names(vif_values))

# Manually add predictor names if they are missing
if (is.null(names(vif_values))) {
  names(vif_values) <- rownames(vif_values)
}

# Check if VIF values are correctly calculated
if (length(vif_values) > 0) {
  # Ensure that vif_values is a named vector
  vif_df <- data.frame(variable = rownames(vif_values), VIF = vif_values[,1])
  
  # Plot VIF values
  vif_plot <- ggplot(vif_df, aes(x = reorder(variable, VIF), y = VIF)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +
    labs(title = "VIF Values of Predictors",
         x = "Predictor Variables",
         y = "VIF")
} else {
  cat("No VIF values to display. The stepwise model may have excluded all predictors.\n")
  vif_plot <- NULL
}

# Create predictions and residuals
estimated <- predict(stepwise_model)
residuals <- residuals(stepwise_model)
predictions_df <- data.frame(estimated = estimated, residuals = residuals)
predictions_df[[response_var]] <- SENIC[[response_var]]

# Plot estimated vs actual values
estimated_vs_actual_plot <- ggplot(predictions_df, aes(x = estimated, y = .data[[response_var]])) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Estimated vs Actual Values",
       x = "Estimated Values",
       y = "Actual Values")

# Plot residuals vs estimated values
residuals_vs_estimated_plot <- ggplot(predictions_df, aes(x = estimated, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Estimated Values",
       x = "Estimated Values",
       y = "Residuals")

# Plot residuals vs actual values
residuals_vs_actual_plot <- ggplot(predictions_df, aes(x = .data[[response_var]], y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Actual Values",
       x = "Actual Values",
       y = "Residuals")

# Print all results at once at the end
cat("Stepwise Model Summary:\n")
print(model_summary)

cat("\nTidy Stepwise Model Coefficients with Confidence Intervals:\n")
print(tidy_model)

cat("\nStepwise Regression Formula:\n", regression_formula, "\n")

cat("\nMultiple R-squared: ", model_summary$r.squared, "\n")
cat("Adjusted R-squared: ", model_summary$adj.r.squared, "\n")

f_value <- model_summary$fstatistic[1]
f_df1 <- model_summary$fstatistic[2]
f_df2 <- model_summary$fstatistic[3]
f_pvalue <- pf(f_value, f_df1, f_df2, lower.tail = FALSE)
cat("F-statistic: ", f_value, " on ", f_df1, " and ", f_df2, " DF, p-value: ", f_pvalue, "\n")

# Print the VIF values if available
if (!is.null(vif_plot)) {
  cat("\nVIF Values:\n")
  print(vif_df)
  print(vif_plot)
}

# Print the other plots
print(estimated_vs_actual_plot)
print(residuals_vs_estimated_plot)
print(residuals_vs_actual_plot)

