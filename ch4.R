# Required packages
library(tidymodels)
library(readxl)
library(tidyverse)
library(here)
library(performance)

dir.create("data",
           showWarnings = FALSE)

dir.create("outputs",
           showWarnings = FALSE)

dir.create("scripts",
           showWarnings = FALSE)

# read the data into R document 
ch4 <- read_xlsx(here("data", "DNA methylation data.xlsm"), sheet = 1)

# Explore the first few rows
head(ch4)

# Clean data 
skim(ch4)

# Check for any duplications 
ch4 |>
  distinct()

#Remove any missing data from the dataset
ch4 <- ch4 |>
  drop_na()


linear_model <- linear_reg() |> 
  set_engine("lm") 

age_recipe <- recipe(Age ~ ., data = ch4) |> 
  step_center(all_predictors())  # Centering the predictors

workflow_model <- workflow() |> 
  add_model(linear_model) |> 
  add_recipe(age_recipe)

library(tidymodels)  # make sure this is at the top

# Example assuming you're predicting 'Age' using all other variables
age_recipe <- recipe(Age ~ ., data = ch4) %>%
  step_center(all_numeric_predictors()) %>%
  step_scale(all_numeric_predictors())

# Example model: linear regression
linear_model <- linear_reg() %>%
  set_engine("lm")

# Build workflow
age_workflow <- workflow() %>%
  add_model(linear_model) %>%
  add_recipe(age_recipe)

# Fit the model
fit_model <- fit(age_workflow, data = ch4)

fit_model_summary <- tidy(fit_result)
fit_model_summary

# Get predictions on the training data

predictions_lm <- augment(fit_model, new_data = ch4)

# Plot observed vs. predicted values
ggplot(data = ch4, aes(x = Age, y = predictions_lm$.pred)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Observed vs Predicted Bat Age")

mse <- mean((ch4$Age - predictions_lm$.pred)^2)

# A cleaner function for calculating MSE
mse_impl <- function(model, data, predictor) {
  augment(model, new_data = data) |> 
    mutate(squared_error = (.pred - {{predictor}})^2) |> 
    summarise(mse = mean(squared_error)) |> 
    pull(mse)
}

mse <-  mse_impl(fit_model, ch4, Age)

rmse <- sqrt(mse)
rmse

# Get more comprehensive model statistics
glance(fit_model)

fit_model |> 
  extract_fit_engine() |> 
  check_model()

# Fit a simple model outside tidymodels
lm_model <- lm(Age ~ ., data = ch4)
performance::check_model(lm_model)

# Split the data into 80% training and 20% testing
set.seed(123)  # For reproducibility
split <- initial_split(ch4, prop = 0.8)

# Extract training and testing sets
train_data <- training(split)
test_data <- testing(split)


glimpse(train_data)
glimpse(test_data)

# Fit the model on the training data
lm_fit <- fit(workflow_model, data = train_data)

age_recipe <- recipe(Age ~ ., data = ch4) %>%
  step_center(all_numeric_predictors()) %>%
  step_scale(all_numeric_predictors())

model <- linear_reg() %>%
  set_engine("lm")

workflow_model <- workflow() %>%
  add_model(model) %>%
  add_recipe(age_recipe)

lm_fit <- fit(workflow_model, data = train_data)

test_preds <- predict(lm_fit, new_data = test_data)

test_metrics <- metrics(test_data, truth = "Age", estimate = test_preds$.pred)


cat("RMSE is:", sqrt(mse_impl(lm_fit, test_data, Age)))

####


lm_predictions <- predict(lm_fit, new_data = test_data)

results <- tibble(
  truth = test_data$Age,  # Actual values (truth)
  estimate = lm_predictions$.pred  # Predicted values (estimate)
)

library(yardstick)
rsq_result <- rsq(results, truth = truth, 
                  estimate = estimate)

# Print the R-squared result
rsq_result

folds <- vfold_cv(train_data, v = 10)

# Fit models using cross-validation
cv_results <- fit_resamples(workflow_model, 
                            resamples = folds)

# Collect and visualize metrics
cv_metrics <- collect_metrics(cv_results)

ggplot(cv_metrics, aes(x = .metric, y = mean, color = .metric)) +
  geom_boxplot() +
  labs(title = "Cross-Validation Performance")

fold_numbers <- seq(2, 10, 1)

cv_results_by_fold <- map_dfr(fold_numbers, function(k_value) {
  # Create cross-validation with k folds
  k_folds <- vfold_cv(train_data, v = k_value)
  
  # Fit models using cross-validation
  cv_results <- fit_resamples(workflow_model, resamples = k_folds)
  
  # Collect metrics and add k value
  cv_metrics <- collect_metrics(cv_results)
  cv_metrics$.k <- k_value
  
  return(cv_metrics)
})

ggplot(cv_results_by_fold, aes(x = .k, y = mean, color = .metric)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ .metric, scales = "free_y") +
  labs(title = "Performance by Number of Cross-Validation Folds",
       x = "Number of Folds (k)",
       y = "Mean Performance")

