##Darby's Extreme Gradient Boosting code using XGBoost, 9/27/24

library(doParallel)
library(tidymodels)
library(xgboost)
library(lubridate)
library(vip)
library(tidyverse)
library(yardstick)
library(caret)


#Load data 

library(tidyverse)
library(dplyr)

flux_rates=read_csv("flux_rates_meta-analysis.csv")

# Specify the order of x-axis categories
desired_order2 <- c("Methanogenesis", "Aerobic methane oxidation", "Anaerobic methane oxidation")

# Convert "Pathway" column to factor with desired order
flux_rates$Pathway <- factor(flux_rates$Pathway, levels = desired_order2)

#create filtered dataframe with all the environmental variables I care about.
flux_rates_filtered <- flux_rates[, c("Pathway", "Mean_rate", "Wetland_type", "Inc_temp", "Inc_length_days","Depth_representative", "Salinity_ppt", "Salinity_cat", "dom_veg", "method","subpathway")]

#keep NAs in categorical values

# remove the continuous salinity column
flux_rates_filtered <- flux_rates_filtered %>%
  dplyr::select(-Salinity_cat)

#hot encoding of categorical values to integers (MLMs generally perform better when hot encoding categorical predictors...cite)
flux_rates_encoded <- flux_rates_filtered %>% mutate( Wetland_type = as.integer(factor(Wetland_type)), dom_veg = as.integer(factor(dom_veg)), method = as.integer(factor(method)), subpathway = as.integer(factor(subpathway)) ) 

#re-organize datasets for modeling of each pathway separately 
MGEN_filtered=flux_rates_encoded[flux_rates_encoded$Pathway == "Methanogenesis", ]
MGEN_filtered$Pathway <- NULL
names(MGEN_filtered) <- c("Mean_rate", "Wetland_type", "Incubation_temperature", "Incubation_length", "Depth", "Salinity", "Dominant_vegetation", "Method", "Subpathway")

MOX_filtered=flux_rates_encoded[flux_rates_encoded$Pathway == "Aerobic methane oxidation", ]
MOX_filtered$Pathway <- NULL
names(MOX_filtered) <- c("Mean_rate", "Wetland_type", "Incubation_temperature", "Incubation_length", "Depth", "Salinity", "Dominant_vegetation", "Method", "Subpathway")
MOX_filtered$Subpathway <- NULL #no subpathways exist for MOx

AOM_filtered=flux_rates_encoded[flux_rates_encoded$Pathway == "Anaerobic methane oxidation", ]
AOM_filtered$Pathway <- NULL
names(AOM_filtered) <- c("Mean_rate", "Wetland_type", "Incubation_temperature", "Incubation_length", "Depth", "Salinity", "Dominant_vegetation", "Method", "Subpathway")


##METHANOGENESIS

# Set variables for model
data <- MGEN_filtered
response <- data$Mean_rate
predictors <- data %>% select(-Mean_rate)


#set up xgboost model
xgb_spec1 <- boost_tree(
  trees = tune(),
  tree_depth = tune(),            
  min_n = tune(),
  loss_reduction = tune (),
  sample_size = tune(),
  mtry = tune (),
  learn_rate = tune(),
  stop_iter = tune()
) %>%
  set_engine("xgboost", validation = 0.2, nthread = ncores) %>%
  set_mode("regression")

# Define the workflow
workflow1 <- workflow() %>%
  add_formula(Mean_rate ~ .) %>%
  add_model(xgb_spec1)

# Define cross-validation folds
folds1 <- vfold_cv(data, v = 5, repeats = 2, strata = NULL)

params <- parameters(xgb_spec1) %>%
  finalize(data)

# Register parallel processing
doParallel::registerDoParallel()
ncores <- parallel::detectCores()


# Perform Bayesian optimization
Xgb_res_bayes <- tune_bayes(
  workflow1,
  resample = folds1,
  param_info = params,  # Use finalized tuning parameters
  initial = 10,          # Adjust as needed
  iter = 30,          # Adjust as needed
  control = control_bayes(save_pred = TRUE, save_workflow = TRUE)
)

final_xgb_bayes$
  
##Now exploring results

##View the data.frame that my model created
Xgb_res_bayes %>%
  collect_metrics() %>%
  View()

Xgb_res_bayes %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, mtry : sample_size) %>%
  pivot_longer (mtry: sample_size, names_to = "parameter", values_to = "value") %>%
  ggplot(aes(value, mean, color = parameter)) + geom_point(show.legend = FALSE) + facet_wrap(~parameter, scales = "free_x")


##Show the best predictors
show_best(Xgb_res, "rsq")
show_best(Xgb_res, "rsq") %>% View()
show_best(Xgb_res_bayes)
##select the model that did the best depending on the hyperparameters
best_hyperparameters <- select_best(Xgb_res_bayes)
?select_best()
Best_rsq <- select_best(Xgb_res, "rsq")
Best_rmse <- select_best(Xgb_res_bayes)
View(Best_rmse)

##finalizing my workflow based off of the best rsq model for all of the hyperparameters.
final_xgb <- finalize_workflow(workflow1, Best_rmse)
final_xgb

final_xgb.mod <- finalize_model(xgb_spec1, best_hyperparameters)

##Of these geom points on the plot from this output is the most important in getting the
final_xgb %>%
  fit(data = data) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point")

final_xgb_bayes <- workflow1 %>%
  finalize_workflow(
    select_best(Xgb_res_bayes)
  ) %>%
  fit(data)

final_xgb$fit$actions
?xgb.train()





final_fit_with_metrics <- final_xgb %>%
  fit_resamples(
    resamples = folds1,
    control = control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )



# Collect RMSE/R2 from the model. output table for model validation.
metrics <- final_fit_with_metrics %>%
  collect_metrics()


# Separate training and validation metrics
training_metrics <- metrics %>%
  filter(.metric == "rmse" & .estimator == "training")

validation_metrics <- metrics %>%
  filter(.metric == "rmse" & .estimator == "validation")



# Combine both training and validation metrics for plotting
combined_metrics <- bind_rows(
  training_metrics %>% mutate(type = "Training"),
  validation_metrics %>% mutate(type = "Validation")
)

# Plot the learning curves for both training and validation
ggplot(combined_metrics, aes(x = .iter, y = mean, color = type)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Learning Curve for XGBoost Model",
    x = "Boosting Iterations (Trees)",
    y = "RMSE",
    color = "Error Type"
  ) +
  theme_minimal()




# Prepare data as DMatrix ###ONLY IF THE ABOVE ISNT WORKING 
dtrain <- xgb.DMatrix(data = as.matrix(data %>% select(-Mean_rate)), label = data$Mean_rate)

# Set parameters
params <- list(
  objective = "reg:squarederror",
  eta = best_hyperparameters$learn_rate,
  max_depth = best_hyperparameters$tree_depth,
  subsample = best_hyperparameters$sample_size,
  colsample_bytree = best_hyperparameters$mtry,
  min_child_weight = best_hyperparameters$min_n
)

# Run cross-validation with watchlist tracking
cv <- xgb.cv(
  params = params,
  data = train,
  nrounds = best_params$trees,
  nfold = 10,  # 10-fold CV
  verbose = TRUE,
  early_stopping_rounds = 10,  # Stop if no improvement
  showsd = TRUE,  # Show standard deviation
  metrics = "rmse",  # Track RMSE
  print_every_n = 10  # Print results every 10 iterations
)


# Extract training and validation error
train_rmse <- cv$evaluation_log$train_rmse_mean
validation_rmse <- cv$evaluation_log$test_rmse_mean
iterations <- seq_along(train_rmse)

# Plot both training and validation learning curves
plot(iterations, train_rmse, type = "l", col = "blue", xlab = "Iterations", ylab = "RMSE", lwd = 2)
lines(iterations, validation_rmse, type = "l", col = "red", lwd = 2)
legend("topright", legend = c("Training RMSE", "Validation RMSE"), col = c("blue", "red"), lwd = 2)






## Different grid search 

set.seed(123)

##Setting up my model specification, early stop** here.
xgb_spec <- boost_tree(trees = tune(),
                       tree_depth = tune(),
                       min_n = tune(),
                       loss_reduction = tune (),
                       sample_size = tune(),
                       mtry = tune (),
                       learn_rate = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("regression")


##Now specifying my space - filling parameter grids
Xgb_grid <- grid_space_filling(trees(), tree_depth(), min_n(), loss_reduction(), sample_size = sample_prop(), finalize(mtry(), data), learn_rate(), size = 50)


##Now training the model with the parameters that I want
workflow <- workflow() %>%
  add_formula(Mean_rate~.) %>%
  add_model(xgb_spec)

##Now we need the data to tune on:
set.seed(123)

folds <-vfold_cv(data, strata = NULL)

## different fold, may be more robust. Give this a try.
set.seed(123)
folds <- vfold_cv(data, v = 5, repeats = 2, strata = NULL)


# Now tuning
doParallel::registerDoParallel()
set.seed(123)
Xgb_res <- tune_grid(workflow, resample = folds, grid = Xgb_grid, control = control_grid (save_pred = TRUE))

train_x = data.matrix(data[, -1])
train_y = data.matrix(data[,1])
test_x = data.matrix(data[, -1])
test_y = data.matrix(data[,1])
length(train_y)
length(train_x)
dtrain = xgb.DMatrix(data = train_x, label = train_y)
dtest = xgb.DMatrix(data = test_x, label = test_y)

Best_rmse <- select_best(Xgb_res)
final_xgb <- finalize_workflow(workflow, Best_rmse)
final_xgb$fit$actions

final_xgb %>%
  fit(data = data) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point")
# Find the best hyperparameters here


plot(final_xgb$fit$fit$fit$evaluation_log)
final_res$.metrics

watchlist = list(train = data, test = data)

xgb_m2 <- xgb.train(data = data, objective = "reg:squarederror", eta = 0.0125159363061706, max_depth = 10, gamma = 1.97332956554775e-07, colsample_bytree = 1, colsample_bynode = 0.4, min_child_weight = 4, subsample = 0.881933937387774, nthread = 1, watchlist = watchlist, nrounds = 2000, early_stopping_rounds = 200)


##Then do cross validations on xgboost model xgb_m2 and plot all of those.



##Now exploring results for the xgb_m2


xgb_ms.df <- as.data.frame(xgb_m2$evaluation_log)

head(xgb_ms.df)


ggplot(xgb_ms.df, aes(x = iter)) +
  geom_line(aes(y=train_rmse, color = "train_rmse")) +
  geom_line(aes(y=test_rmse, color = "test_rmse")) + ylab("Training and Testing") + ggtitle('Learning Curve') + scale_color_manual(values = c("tan", "salmon")) +
  theme_minimal() +
  theme(panel.grid = element_blank())




final_fit_with_metrics <- final_xgb %>%
  fit_resamples(
    resamples = folds,
    control = control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )

range(data$Mean_rate)
sd(data$Mean_rate)

# Collect metrics from the model. RMSE/R2 for model validation table.
metrics <- final_fit_with_metrics %>%
  collect_metrics()
