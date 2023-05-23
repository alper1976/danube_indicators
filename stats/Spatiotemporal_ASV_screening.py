# ASV screening as bioindicators by spatiotemporal modeling
# Author: Laurent Fontaine

# Set page width display
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))

# Load packages
from os import listdir
from os.path import isfile, join
import pandas as pd
from pandas import DataFrame
import numpy as np
import re
import sys
from random import sample
from xgboost import XGBRegressor
from xgboost import XGBClassifier
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from datetime import datetime
from datetime import timedelta
import plotly.express as px
import plotly.graph_objs as go

import warnings
warnings.filterwarnings("ignore")

# Scale variables into z-scores by mean and standard deviation for a given spatial window of n sites
def scale(df, n):
    df_scaled = np.where(df.rolling(n).std()!=0, (df - df.rolling(n).mean())/df.rolling(n).std(), df - df.rolling(n).mean())
    return df_scaled

# Format the response variable into a vector of next site values
def response_var(df):
    df_copy = df.copy()
    var = ''.join(df.columns)
    df_copy[var + '_next_site'] = df_copy.iloc[:,[0]]
    df_copy.iloc[0:df.shape[0]-1:1,[1]] = df_copy.iloc[1:df.shape[0]:1,[0]]
    df_copy.drop([var], inplace=True, axis=1)
    return(df_copy)

# Format explanatory variables into their appropriate lagged and scaled versions
def var_scaling(df, lag, n_sites_rolling_avg):
    var = ''.join(df.columns)
    df_copy = df.copy()
    if n_sites_rolling_avg==1:
        df_copy[var + '_scaled'] = df.copy()
    else:
        df_copy[var + '_scaled'] = scale(df, n_sites_rolling_avg).copy()
    for i in range(1, lag+1, 1):
        df_copy[var + "_scaled_lag_" + str(i)] = df_copy.iloc[:,[1]]
        df_copy.iloc[i:df.shape[0]:1,[i+1]] = df_copy.iloc[0:df.shape[0]-i:1,1]
    df_copy.drop(var, inplace=True, axis=1)
    return(df_copy)

# model_grid_search() is the main function for assessing each combination of ASVs in corresponding lag and scaling to find useful ones
def model_grid_search(df, df_all, response_variable, response_variable_best_settings, highest_n_first_rows_to_remove, depth=8, rate=0.075, weight=3):
    # Divide data into train and test sets
    df_X_index = []
    train_index = []
    test_index = []
    set_vec = []
    for i in np.unique(df['Side']):
        start = (df.Side.values == i).argmax() + highest_n_first_rows_to_remove
        s = pd.Series(df['Side']==i)
        stop = s.where(s).last_valid_index()
        df_index = list(range(start, stop + 1))
        train_length = int((len(df_index)-1)*train_share)
        test_length = len(df_index)- 1 - train_length
        train_index = train_index + list(range(start, start + train_length))
        test_index = test_index + list(range(start + train_length, start + train_length + test_length))

    df_Y = pd.DataFrame(df_all[response_variable + '_next_site'])
        
    X_train = df[:].iloc[train_index, :]
    X_train['Side'] = encoded_name.fit_transform(X_train['Side'].astype(str))
    X_train.drop(['Site'], inplace=True, axis=1)
    y_train = df_Y[:].iloc[train_index, :]
    X_test = df[:].iloc[test_index, :]
    X_test['Side'] = encoded_name.fit_transform(X_test['Side'].astype(str))
    X_test.drop(['Site'], inplace=True, axis=1)
    y_test = df_Y[:].iloc[test_index, :]

    # XGboost model
    model_seed = 100
    n_estimators = df_all.shape[1]    # Number of boosted trees to fit. default = 100
    max_depth = depth                 # Maximum tree depth for base learners. default = 3
    learning_rate = rate              # Boosting learning rate (xgb’s “eta”). default = 0.1
    min_child_weight = weight         # Minimum sum of instance weight(hessian) needed in a child. default = 1
    subsample = 1                     # Subsample ratio of the training instance. default = 1
    colsample_bytree = 1              # Subsample ratio of columns when constructing each tree. default = 1
    colsample_bylevel = 1             # Subsample ratio of columns for each split, in each level. default = 1
    gamma = 1.1

    # For classification, use 'multi:softprob' objective; for regression, use 'reg:squarederror'
    model = XGBClassifier(objective ='multi:softprob', 
                         booster='dart',
                         rate_drop=0.1,
                         skip_drop=0.5,
                         seed=model_seed,
                         n_estimators=n_estimators,
                         max_depth=max_depth,
                         learning_rate=learning_rate,
                         min_child_weight=min_child_weight,
                         subsample=subsample,
                         colsample_bytree=colsample_bytree,
                         colsample_bylevel=colsample_bylevel,
                         gamma=gamma)

    # Train the model
    # For classification, use 'mlogloss' eval metric; for regression, use 'rmse'
    model.fit(X_train, y_train, eval_set=[(X_train, y_train), (X_test, y_test)], #Training the model
              eval_metric='mlogloss', early_stopping_rounds=5, verbose=False)

    #Get the best model based on logloss or mean squared error
    XGboost_best_model_index = model.best_iteration #Index of the iteration yielding the lowest test error
    XGboost_best_iteration = model.get_booster().best_ntree_limit #Selection of the model by index of best iteration
    MSE_per_epoch = model.evals_result()  #Pulling vectors of prediction error for training and test sets for each iteration
    
    # Make predictions for the test set
    pred = model.predict(X_test, ntree_limit=XGboost_best_iteration)

    # Assemble true vs predicted values for the test set and return 
    true_vs_pred = pd.DataFrame(df_all[:].iloc[test_index, :][response_variable], columns=[response_variable])    
    true_vs_pred['pred_'+response_variable] = pred
    
    # Assess accuracy of the model; model_r_squared variable can be set to accuracy, r_squared, rmse 
    model_r_squared = accuracy_score(true_vs_pred[response_variable], true_vs_pred['pred_' + response_variable])    
    return pred, pred, round(model_r_squared, 5)

### Outside loop ###
# Loading data
file_dir = 'path/Data/' # Set path to files 'j3s_L.df.export.tsv', 'j3s_M.df.export.tsv', 'j3s_R.df.export.tsv'
files = [f for f in listdir(file_dir) if isfile(join(file_dir, f))]
files_dict = {}

for file in files:
    symbol_str = file.strip(".df.export.tsv")
    symbol_df = pd.read_csv(file_dir + file, sep='\t')
    files_dict[symbol_str] = symbol_df

files_list = list(files_dict.keys())
files_dict_array = np.array([files_dict[x] for x in files_list], dtype=object)
df_all = pd.concat(files_dict_array, axis=0)
df_all = pd.concat([df_all.reset_index(drop=True), response_var(df_all.loc[:,['MHS_WQ']]).reset_index(drop=True)], axis=1)

response_variable = 'MHS_WQ' # Set response variable: 'MHS_WQ' for ecological status; 'Chl_a' for chlorophyll a
train_share = 0.85 # How much of the dataset is allocated for training the model; the rest is available for testing it.

# Encode categorical variables
df_all_columns = [x for x in list(df_all.columns)]
df_all_columns_mask = [x.split('_')[0]=='MHS' for x in df_all_columns]
df_all_response_variable_columns = list(df_all.iloc[:, df_all_columns_mask].columns)
encoded_name = preprocessing.LabelEncoder()
for categorical_variable in df_all_response_variable_columns:
    df_all[categorical_variable] = encoded_name.fit_transform(df_all[categorical_variable].astype(str))

# Use 'MHS_WQ_next_site' for classification objective; 'Chl_a_next_site' for regression objective
response_variables = ['MHS_WQ_next_site']
metadata_variables = ["pH", "Conductivity", "Water_Temperature", "TSS_CHECK_comparability_JDS2vs3",
"DOC", "Ptot", "Carbon", "Ntot", 'Ptot', "Log1_E.coli", "Log1_Total_Coliforms", "Log1_AllBac", 'impounded',
"Site_ID", 'Site', 'Side', 'River_km_Distance_to_mouth', 'Chl_a', 'MHS_WQ', 'Water_temperature']

# ASVs are the screening variables; metadata variables are not used for building the design matrix
screening_variables = [x for x in list(df_all.columns) if x not in metadata_variables + response_variables]

# Set number models to generate
for random_start in range(10000):
    best_settings_dict = {} # Store best lag and scaling settings for each variable and update with each iteration
    highest_n_first_rows_to_remove = {} # Store values in a dictionary in order to update for each pass on a given variable
    
    #Set random model hyperparameters prior to optimization after design matrix optimization
    best_max_depth = sample(list(range(6, 9)), 1)[0]
    best_learning_rate = sample(list(np.round(np.arange(0.05, 1.005, 0.005),3)), 1)[0]
    best_min_child_weight = sample(list(range(1, 8)), 1)[0]
    print('max_depth, learning_rate, min_child_weight:', best_max_depth, best_learning_rate, best_min_child_weight)
    
    df = df_all.loc[:,['Site', 'Side', 'River_km_Distance_to_mouth']].copy()
    global_model_optimal = False
    design_matrix_optimal = False
    max_r_squared_reached = False
    iteration = 0
    global_start_time = datetime.now() ### Execution time start

    if random_start == 0: # Run only the first time when the most relevant ASVs should appear by decreasing order.
        variable_list = screening_variables.copy()
    else:
        variable_list = sample(screening_variables, k=len(screening_variables))
    print(variable_list)
    
    # Run variable screening process for ASVs and response variable
    while global_model_optimal == False or design_matrix_optimal==False:
        print('Iteration', iteration + 1)
        print('Design matrix optimization')
        start_time = datetime.now() ### Execution time start
        design_matrix_optimal = True

        # Optimize for response variable
        if response_variable in list(highest_n_first_rows_to_remove.keys()):
            # Remove from the completed design matrix all existing columns for the response variable to check if its lag and scaling are optimal.
            df_columns = [x for x in list(df.columns)]
            df_columns_mask = [x.split('_scaled')[0]==response_variable for x in df_columns]
            df_drop_columns = list(df.iloc[:, df_columns_mask].columns)
            df.drop(df_drop_columns, inplace=True, axis=1)

        df_temp = df.copy()

        response_variable_settings_list = []
        response_variable_r_squared_list = []

        for lag in range(0,1): # Set number of sites of lag
            for n_sites_rolling_avg in range(1,2):  # Set number of sites of scaling. No scaling for qualitative response variable.
                highest_n_first_rows_to_remove[response_variable] = lag + n_sites_rolling_avg - 1
                response_variable_best_settings = (lag, n_sites_rolling_avg)
                model_output = model_grid_search(df_temp, df_all, response_variable, response_variable_best_settings, max(list(highest_n_first_rows_to_remove.values())), depth=best_max_depth, rate=best_learning_rate, weight=best_min_child_weight)
                model_r_squared = model_output[2]
                
                # Store variable settings and model performance
                response_variable_settings_list.append((lag, n_sites_rolling_avg))
                response_variable_r_squared_list.append(model_r_squared)

        # Select the best model by performance and store corresponding best variable settings
        response_variable_best_settings = response_variable_settings_list[response_variable_r_squared_list.index(max(response_variable_r_squared_list))]
        response_variable_best_r_squared = max(response_variable_r_squared_list)
        highest_n_first_rows_to_remove[response_variable] = response_variable_best_settings[0] + response_variable_best_settings[1] - 1  # n_first_rows_to_remove = lag + n_sites_rolling_avg - 1

        # Only move on with model hyperparameter tuning after an iteration did not yield changes to the design matrix.
        if iteration > 0:
            if best_settings_dict[response_variable] != (response_variable_best_settings[0], response_variable_best_settings[1]):
                design_matrix_optimal = False
                global_model_optimal = False
        best_settings_dict[response_variable] = (response_variable_best_settings[0], response_variable_best_settings[1])

        print("Response variable:", response_variable)
        print("Best settings:", response_variable_best_settings)
        print("Best accuracy:", response_variable_best_r_squared, '\n')

        global_best_r_squared = response_variable_best_r_squared

        # Optimize for explanatory variables
        for i in variable_list:
            variable_settings_list = []
            variable_r_squared_list = []

            # Remove from the completed design matrix all existing columns for a feature to check if it is still useful
            removed_variable_model_r_squared = 0
            df_columns = [x for x in list(df.columns)]
            df_columns_mask = [x.split('_scaled')[0]==i for x in df_columns]
            df_drop_columns = list(df.iloc[:, df_columns_mask].columns)
            df.drop(df_drop_columns, inplace=True, axis=1)
            df_temp.drop(df_drop_columns, inplace=True, axis=1)
            if i in [x.split('_scaled')[0] for x in df_columns]:
                previous_iteration_feature_settings = best_settings_dict[i]
                previous_iteration_feature_n_first_rows_to_remove = highest_n_first_rows_to_remove[i]
                del highest_n_first_rows_to_remove[i]
                del best_settings_dict[i]

                # Compute accuracy/r-squared value for model after removing feature
                removed_variable_model_r_squared = model_grid_search(df, df_all, response_variable, response_variable_best_settings, max(list(highest_n_first_rows_to_remove.values())), depth=best_max_depth, rate=best_learning_rate, weight=best_min_child_weight)[2]
            
            # Test combinations of lag and scaling for a given ASV
            for lag in range(0,5):
                for n_sites_rolling_avg in range(1,15):
                    highest_n_first_rows_to_remove[i] = lag + n_sites_rolling_avg - 1

                    df_temp = df.copy()            
                    df_temp = pd.concat([df.reset_index(drop=True), var_scaling(df_all.loc[:,[i]], lag=lag, n_sites_rolling_avg=n_sites_rolling_avg).reset_index(drop=True)], axis=1)

                    model_output = model_grid_search(df_temp, df_all, response_variable, response_variable_best_settings, max(list(highest_n_first_rows_to_remove.values())), depth=best_max_depth, rate=best_learning_rate, weight=best_min_child_weight)

                    del highest_n_first_rows_to_remove[i]

                    model_r_squared = model_output[2]

                    variable_settings_list.append((lag, n_sites_rolling_avg))
                    variable_r_squared_list.append(model_r_squared)

            variable_best_settings = variable_settings_list[variable_r_squared_list.index(max(variable_r_squared_list))]
            variable_best_r_squared = max(variable_r_squared_list)

            # If the r-squared value is better from removing or altering an existing ASV, remove and change it, otherwise keep as is.
            if variable_best_r_squared > global_best_r_squared + 0.01:
                print(i, 'updated/added to design matrix.')
                df = pd.concat([df.reset_index(drop=True), var_scaling(df_all.loc[:,[i]], lag=variable_best_settings[0], n_sites_rolling_avg=variable_best_settings[1]).reset_index(drop=True)], axis=1)
                global_best_r_squared = variable_best_r_squared
                highest_n_first_rows_to_remove[i] = variable_best_settings[0] + variable_best_settings[1] - 1  # n_first_rows_to_remove = lag + n_sites_rolling_avg - 1
                best_settings_dict[i] = (variable_best_settings[0], variable_best_settings[1])

                # Only move on with model hyperparameter tuning after an iteration did not yield changes to the design matrix.
                if i not in [x.split('_scaled')[0] for x in df_columns]:
                    design_matrix_optimal = False
                    global_model_optimal = False
                elif best_settings_dict[i] != previous_iteration_feature_settings:
                    design_matrix_optimal = False
                    global_model_optimal = False

                print("Variable:", i)
                print("Best settings:", variable_best_settings)
                print("Best accuracy:", variable_best_r_squared, '\n')

            # Keep ASV if useful and its lag and scaling are unchanged
            elif iteration > 0 and removed_variable_model_r_squared < global_best_r_squared and i in [x.split('_scaled')[0] for x in df_columns]: # 
                print('Keeping variable settings from last iteration for ', i)
                df = pd.concat([df.reset_index(drop=True), var_scaling(df_all.loc[:,[i]], lag=previous_iteration_feature_settings[0], n_sites_rolling_avg=previous_iteration_feature_settings[1]).reset_index(drop=True)], axis=1)
                best_settings_dict[i] = previous_iteration_feature_settings
                highest_n_first_rows_to_remove[i] = previous_iteration_feature_n_first_rows_to_remove
            # Remove ASV if no longer useful
            elif iteration > 0 and removed_variable_model_r_squared > global_best_r_squared:
                print(i, 'removed from design matrix.')
                design_matrix_optimal = False
                global_model_optimal = False

            # End the screening process when a perfectly accurate model has been achieved
            if global_best_r_squared == 1.0:
                global_model_optimal = True
                if max_r_squared_reached == False:
                    max_r_squared_reached = True
                    variable_list = list(np.unique([x.split('_scaled')[0] for x in df_columns if x not in metadata_variables]))
                    break

        # Save time by skipping further ASV screening for poorly performing starting design matrices and model hyperparameters
        if iteration==0 and global_best_r_squared < 0.6:
            break
        if iteration > 0 and global_best_r_squared < 0.8:
            break
    
        # Proceed to model hyperparameter optimization when the design matrix cannot be improved further
        if design_matrix_optimal == True and global_model_optimal == False:
            print('Model optimization')
            max_depth_list = []
            learning_rate_list = []
            min_child_weight_list = []
            model_r_squared_list = []
            for depth in range(6,9):
                for rate in np.round(np.arange(0.05, 0.605, 0.005),3):
                    for weight in range(1,3):
                        max_depth_list.append(depth)
                        learning_rate_list.append(rate)
                        min_child_weight_list.append(weight)
                        model_r_squared_list.append(model_grid_search(df, df_all, response_variable, response_variable_best_settings, max(list(highest_n_first_rows_to_remove.values())), depth, rate, weight)[2])
            best_model_r_squared_index = model_r_squared_list.index(max(model_r_squared_list))
            best_optimization_max_depth = max_depth_list[best_model_r_squared_index]
            best_optimization_learning_rate = learning_rate_list[best_model_r_squared_index]
            best_optimization_min_child_weight = min_child_weight_list[best_model_r_squared_index]
            if best_optimization_max_depth==best_max_depth and best_optimization_learning_rate==best_learning_rate and best_optimization_min_child_weight==best_min_child_weight:
                global_model_optimal = True
            else:
                best_max_depth = best_optimization_max_depth
                best_learning_rate = best_optimization_learning_rate
                best_min_child_weight = best_optimization_min_child_weight
                design_matrix_optimal = False
                global_model_optimal = True

        finish_time = datetime.now()
        print('Best settings for iteration', iteration + 1)
        for feature in [x for x in list(best_settings_dict.items())]:
            print(feature[0], feature[1])
        print('Best max depth', best_max_depth)
        print('Best learning rate', best_learning_rate)
        print('Best min child weight', best_min_child_weight)
        print('\n', 'Best accuracy for iteration', iteration + 1, global_best_r_squared, '\n')
        print('Iteration', iteration + 1, 'took', (finish_time - start_time).total_seconds(), 'seconds.', '\n')
        iteration = iteration + 1

    global_finish_time = datetime.now()
    print('Best settings')
    for feature in [x for x in list(best_settings_dict.items())]:
        print(feature[0], feature[1])
    print('\n', 'Best accuracy', global_best_r_squared, '\n')
    print('Entire optimization took', (global_finish_time - global_start_time).total_seconds(), 'seconds.')
    
    # Write results to file
    f = open('Danube_WQ_predictors.txt', 'a') # When running the script for ecological status
	#f = open('Danube_chl_a_predictors.txt', 'a') # When running the script for chlorophyll a
    f.write('\n')
    f.write('\t'.join([str(x) for x in best_settings_dict.keys() if x not in [response_variable]]))
    f.close()

    f = open('Danube_WQ_dictionaries.txt', 'a') # When running the script for ecological status
	#f = open('Danube_chl_a_dictionaries.txt', 'a') # When running the script for chlorophyll a
    f.write('\n')
    f.write('\t'.join([str(global_best_r_squared), str(best_settings_dict), str([best_max_depth, best_learning_rate, best_min_child_weight])]))
    f.close()
