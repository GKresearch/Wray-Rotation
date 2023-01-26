#!/usr/bin/env python
# coding: utf-8

# In[23]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Documentation:
Final Updated Artificial Neural Network script for a rotation project in the Wray Lab
author=GKennedy 
Date=2022 December 12th

Notes:
    *This was developed using online resources for information
        -https://thinkingneuron.com/using-artificial-neural-networks-for-regression-in-python/
    **Uses Tensorflow and Keras for the artificial neural network. Originally this script was designed
    To be included in the full machine learning harness, but due to environment set-up issues, they are 
    separate.
"""

import tensorflow
import keras
import pandas as pd
import numpy as np
#Importing sklearn modules
from sklearn.metrics import mean_squared_error,confusion_matrix, precision_score, recall_score, auc,roc_curve
from sklearn import ensemble, linear_model, neighbors, svm, tree, neural_network
from sklearn.linear_model import Ridge
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn import svm,model_selection, tree, linear_model, neighbors, naive_bayes, ensemble, discriminant_analysis, gaussian_process
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import make_column_transformer
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error


# In[2]:


from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.losses import MeanSquaredLogarithmicError
from tensorflow.keras.optimizers import Adam
import matplotlib.pyplot as plt


# In[3]:


#Read in data
df2=pd.read_csv('C:/Users/GK/Documents/WrayLab/raw_data/OCRs_inTADs_unique.csv')
df2.isnull().sum()


# In[59]:


df5=pd.read_csv('C:/Users/GK/Documents/WrayLab/raw_data/pipeline_data/x0011_Gene-Element_TAD_Domain_Check_Confirmed_Promoter-Enhancer_2022-12-14.csv')
df5.isnull().sum()


# In[60]:


#Assign df to be tested
test=df5#.iloc[1:5000,::]


# In[61]:


# here we'll start by using wgCERES_score_nosig as the response vector, and loading in the important factors
x = test[["DHS_prop_repeat", 
                    "DHS_prop_GC", "DHS_length", "n_SNV_Zhou_per_bp", 
                    "distanceToTSS", "zeta.human", "zeta.chimp", "PP_con", "PP_acc", 
                    "PhastCons",
                    "chromHMM_cat_longest", 
                    "annotation", "PhyloP_primates_score"]]
y = test["wgCERES_score_nosig"]#.values
x = pd.get_dummies(x, columns = ['chromHMM_cat_longest','annotation'])


# In[62]:


#Split the data between test and train
x_train, x_test, y_train, y_test=train_test_split(x,y,test_size=0.2, random_state=0)


# In[50]:


#Scale datasets (later functions are using unscaled data so if you want scaled, change as appropriate with the scaled dataset below)
##Note: If you scale remember that you may need to unscale at the end for appropriate comparisons
def scale_datasets(x_train, x_test):
  standard_scaler = StandardScaler()
  x_train_scaled = pd.DataFrame(
      standard_scaler.fit_transform(x_train),
      columns=x_train.columns
  )
  x_test_scaled = pd.DataFrame(
      standard_scaler.transform(x_test),
      columns = x_test.columns
  )
  return x_train_scaled, x_test_scaled
x_train_scaled, x_test_scaled = scale_datasets(x_train, x_test)


# In[63]:


#First Test, with 3 units and sequential model
hidden_units1 = 160
hidden_units2 = 480
hidden_units3 = 256
learning_rate = 0.01
# Creating model using the Sequential in tensorflow
def build_model_using_sequential():
  model = Sequential([
    Dense(hidden_units1, kernel_initializer='normal', activation='relu'),
    Dropout(0.2),
    Dense(hidden_units2, kernel_initializer='normal', activation='relu'),
    Dropout(0.2),
    Dense(hidden_units3, kernel_initializer='normal', activation='relu'),
    Dense(1, kernel_initializer='normal', activation='linear')
  ])
  return model
# build the model
model = build_model_using_sequential()


# In[52]:


# loss function to examine errors 
msle = MeanSquaredLogarithmicError()


# In[1]:


#Compile and run model
model.compile(
    loss=msle, 
    optimizer=Adam(learning_rate=learning_rate), 
    metrics=[msle]
)
# train the model
model_run1 = model.fit(
    x_train_scaled.values, 
    y_train.values, 
    epochs=10, 
    batch_size=64,
    validation_split=0.2
)


# In[29]:


def plot_history(model_run1, key):
  plt.plot(model_run1.model_run1[key])
  plt.plot(model_run1.model_run1['val_'+key])
  plt.xlabel("Epochs")
  plt.ylabel(key)
  plt.legend([key, 'val_'+key])
  plt.show()
# Plot the history
plot_history(model_run1, 'mean_squared_logarithmic_error')


# In[32]:


#Attempt to account for problemeatic hyperparameters to ensure best practices
#This attempt only included 3 parameters to test, and took several hours. Plan to run in background
#Or on remote computing cluster in future.
batch_size_list=[8, 64, 128, 256]
epoch_list=[5, 10, 25, 50, 100]
SearchResultsData=pd.DataFrame(columns=['TrialNumber','Parameters','Accuracy'])

TrialNumber=0
for batch_size_trial in batch_size_list:
    for epochs_trial in epoch_list:
        TrialNumber+=1
        model = build_model_using_sequential()
        model.compile(
        loss=msle, 
        optimizer=Adam(learning_rate=learning_rate), 
        metrics=[msle]
        )
        model.fit(x_train, y_train ,batch_size = batch_size_trial, epochs = epochs_trial, verbose=0)
        MAPE = np.mean(100 * (np.abs(y_test-model.predict(x_test))/y_test))
            
        # printing the results of the current iteration
        print(TrialNumber, 'Parameters:','batch_size:', batch_size_trial,'-', 'epochs:',epochs_trial, 'Accuracy:', 100-MAPE)
        SearchResultsData=SearchResultsData.append(pd.DataFrame(data=[[TrialNumber, str(batch_size_trial)+'-'+str(epochs_trial), 100-MAPE]],
                                                                    columns=['TrialNumber', 'Parameters', 'Accuracy']))


# In[9]:


#Additional attempts to optimize hyperparameters with a function this time and Grid Search CV
def make_regression_ann(Optimizer_trial):
    #from keras.models import Sequential
    #from keras.layers import Dense
    model = build_model_using_sequential()
    model.compile(
    loss=msle, 
    optimizer=Optimizer_trial,#Adam(learning_rate=learning_rate), 
        metrics=[msle]
    )
    return model


# In[18]:


###########################################
from sklearn.model_selection import GridSearchCV
from keras.wrappers.scikit_learn import KerasRegressor

# Listing all the parameters to try
Parameter_Trials={'batch_size':[8, 64, 128],
                  'epochs':[10, 50, 100],
                  'Optimizer_trial':['adam', 'rmsprop']
                 }


# In[19]:


# Creating the regression ANN model
RegModel=KerasRegressor(make_regression_ann, verbose=0)


# In[20]:


###########################################
from sklearn.metrics import make_scorer
#Defining a custom function to calculate accuracy
def Accuracy_Score(orig,pred):
    MAPE = np.mean(100 * (np.abs(orig-pred)/orig))
    print('#'*70,'Accuracy:', 100-MAPE)
    return(100-MAPE)

custom_Scoring=make_scorer(Accuracy_Score, greater_is_better=True)


# In[21]:


#########################################
# Creating the Grid search space
# See different scoring methods by using sklearn.metrics.SCORERS.keys()
grid_search=GridSearchCV(estimator=RegModel, 
                         param_grid=Parameter_Trials, 
                         scoring=custom_Scoring, 
                         cv=5)


# In[24]:


#########################################
# Measuring how much time it took to find the best params
import time
StartTime=time.time()

# Running Grid Search for different paramenters
grid_search.fit(x_test,y_test, verbose=1)

EndTime=time.time()
print("########## Total Time Taken: ", round((EndTime-StartTime)/60), 'Minutes')

print('### Printing Best parameters ###')
grid_search.best_params_


# In[25]:


grid_search.best_params_


# In[27]:


grid_search.cv_results_['params']


# In[29]:


def plot_grid_search(cv_results, grid_param_1, grid_param_2, name_param_1, name_param_2):
    # Get Test Scores Mean and std for each grid search
    scores_mean = cv_results['mean_test_score']
    scores_mean = np.array(scores_mean).reshape(len(grid_param_2),len(grid_param_1))

    scores_sd = cv_results['std_test_score']
    scores_sd = np.array(scores_sd).reshape(len(grid_param_2),len(grid_param_1))

    # Plot Grid search scores
    _, ax = plt.subplots(1,1)

    # Param1 is the X-axis, Param 2 is represented as a different curve (color line)
    for idx, val in enumerate(grid_param_2):
        ax.plot(grid_param_1, scores_mean[idx,:], '-o', label= name_param_2 + ': ' + str(val))

    ax.set_title("Grid Search Scores", fontsize=20, fontweight='bold')
    ax.set_xlabel(name_param_1, fontsize=16)
    ax.set_ylabel('CV Average Score', fontsize=16)
    ax.legend(loc="best", fontsize=15)
    ax.grid('on')


# In[30]:


# Calling Method 
plot_grid_search(grid_search.cv_results_, n_estimators, max_features, 'N Estimators', 'Max Features')


# In[31]:


grid_search.cv_results_['mean_test_score']


# In[32]:


plot.grid_search(grid_search.cv_results_, change='batch_size', kind='bar')


# In[33]:


from sklearn_evaluation import plot


# In[34]:


def GridSearch_table_plot(grid_clf, param_name,
                          num_results=15,
                          negative=True,
                          graph=True,
                          display_all_params=True):

    '''Display grid search results

    Arguments
    ---------

    grid_clf           the estimator resulting from a grid search
                       for example: grid_clf = GridSearchCV( ...

    param_name         a string with the name of the parameter being tested

    num_results        an integer indicating the number of results to display
                       Default: 15

    negative           boolean: should the sign of the score be reversed?
                       scoring = 'neg_log_loss', for instance
                       Default: True

    graph              boolean: should a graph be produced?
                       non-numeric parameters (True/False, None) don't graph well
                       Default: True

    display_all_params boolean: should we print out all of the parameters, not just the ones searched for?
                       Default: True

    Usage
    -----

    GridSearch_table_plot(grid_clf, "min_samples_leaf")

                          '''
    from matplotlib      import pyplot as plt
    from IPython.display import display
    import pandas as pd

    clf = grid_clf.best_estimator_
    clf_params = grid_clf.best_params_
    if negative:
        clf_score = -grid_clf.best_score_
    else:
        clf_score = grid_clf.best_score_
    clf_stdev = grid_clf.cv_results_['std_test_score'][grid_clf.best_index_]
    cv_results = grid_clf.cv_results_

    print("best parameters: {}".format(clf_params))
    print("best score:      {:0.5f} (+/-{:0.5f})".format(clf_score, clf_stdev))
    if display_all_params:
        import pprint
        pprint.pprint(clf.get_params())

    # pick out the best results
    # =========================
    scores_df = pd.DataFrame(cv_results).sort_values(by='rank_test_score')

    best_row = scores_df.iloc[0, :]
    if negative:
        best_mean = -best_row['mean_test_score']
    else:
        best_mean = best_row['mean_test_score']
    best_stdev = best_row['std_test_score']
    best_param = best_row['param_' + param_name]

    # display the top 'num_results' results
    # =====================================
    display(pd.DataFrame(cv_results)             .sort_values(by='rank_test_score').head(num_results))

    # plot the results
    # ================
    scores_df = scores_df.sort_values(by='param_' + param_name)

    if negative:
        means = -scores_df['mean_test_score']
    else:
        means = scores_df['mean_test_score']
    stds = scores_df['std_test_score']
    params = scores_df['param_' + param_name]

    # plot
    if graph:
        plt.figure(figsize=(8, 8))
        plt.errorbar(params, means, yerr=stds)

        plt.axhline(y=best_mean + best_stdev, color='red')
        plt.axhline(y=best_mean - best_stdev, color='red')
        plt.plot(best_param, best_mean, 'or')

        plt.title(param_name + " vs Score\nBest Score {:0.5f}".format(clf_score))
        plt.xlabel(param_name)
        plt.ylabel('Score')
        plt.show()


# In[37]:


GridSearch_table_plot(grid_search, "Optimizer_trial", negative=False)


# In[64]:


model = build_model_using_sequential()
model.compile(
        loss=msle, 
        optimizer=Adam(learning_rate=learning_rate), 
        metrics=[msle]
        )

#model.fit(x_train, y_train ,batch_size = batch_size_trial, epochs = epochs_trial, verbose=0)
# train the model
#history = 
model.fit(
    x_train.values, 
    y_train.values, 
    epochs=100, 
    batch_size=64,
    validation_split=0.2
)


# In[41]:


model.predict(x_test).flatten()


# In[54]:


MAPE = np.mean(100 * (np.abs(y_test-model.predict(x_test).flatten())/y_test))


# In[43]:


MAPE


# In[55]:


mean_absolute_error(model.predict(x_test).flatten(), y_test)


# In[56]:


np.sqrt(mean_squared_error(model.predict(x_test).flatten(), y_test))


# In[65]:


model.predict(x_test)


# In[69]:


y_test


# In[ ]:




