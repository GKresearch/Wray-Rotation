#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from sklearn.svm import SVC
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import make_column_transformer
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error


# In[2]:


import pandas as pd


# In[3]:


from sklearn.svm import SVC


# In[4]:


import numpy as np


# In[5]:


import matplotlib.pyplot as plt
import seaborn as sns


# In[6]:


from sklearn.ensemble import RandomForestClassifier


# In[7]:


import tensorflow


# In[8]:


import keras


# In[9]:


df2=pd.read_csv('C:/Users/Ictinike/Documents/WrayLab/raw_data/OCRs_inTADs_unique.csv')
df2.isnull().sum()


# In[10]:


test=df2#.iloc[1:5000,::]


# In[11]:


# here we'll start by using wgCERES_score_nosig as the response vector,
x = test[["DHS_prop_repeat", 
                    "DHS_prop_GC", "DHS_length", "n_SNV_Zhou_per_bp", 
                    "distanceToTSS", "zeta.human", "zeta.chimp", "PP_con", "PP_acc", 
                    "PhastCons",
                    "chromHMM_cat_longest", 
                    "annotation", "PhyloP_primates_score"]]
y = test["wgCERES_score_nosig"]#.values
x = pd.get_dummies(x, columns = ['chromHMM_cat_longest','annotation'])


# In[12]:


random_state=None


# In[13]:


x_train, x_test, y_train, y_test=train_test_split(x,y,test_size=0.2, random_state=0)


# In[14]:


### Sandardization of data ###
#from sklearn.preprocessing import StandardScaler
#PredictorScaler=StandardScaler()
#TargetVarScaler=StandardScaler()
# Storing the fit object for later reference
#PredictorScalerFit=PredictorScaler.fit(x)
#TargetVarScalerFit=TargetVarScaler.fit(y)
## Generating the standardized values of X and y
#x=PredictorScalerFit.transform(x)
#y=TargetVarScalerFit.transform(y)
y_test


# In[15]:


from keras.models import Sequential
from keras.layers import Dense


# In[16]:


# create ANN model
model = Sequential()


# In[17]:


# Defining the Input layer and FIRST hidden layer, both are same!
model.add(Dense(units=5, input_dim=32, kernel_initializer='normal', activation='relu'))


# In[18]:


# Defining the Second layer of the model
# after the first layer we don't have to specify input_dim as keras configure it automatically
model.add(Dense(units=5, kernel_initializer='normal', activation='tanh'))


# In[19]:


# The output neuron is a single fully connected node 
# Since we will be predicting a single number
model.add(Dense(1, kernel_initializer='normal'))


# In[20]:


model.add(Dense(256, input_shape=(x_train.shape[1],), activation='sigmoid'))


# In[21]:


# Compiling the model
model.compile(loss='mean_squared_error', optimizer='adam')


# In[22]:


#y_train


# In[23]:


# Fitting the ANN to the Training set
#model.fit(x_train, y_train ,batch_size = 20, epochs = 50)


# In[24]:


def FunctionFindBestParams(X_train, y_train, X_test, y_test):
    # Defining the list of hyper parameters to try
    batch_size_list=[5, 10, 15, 20]
    epoch_list=[5, 10, 50, 100]
    
    SearchResultsData=pd.DataFrame(columns=['TrialNumber','Parameters','Accuracy'])
    # initializing the trials
    TrialNumber=0
    for batch_size_trial in batch_size_list:
        for epochs_trial in epoch_list:
            TrialNumber+=1
            # create ANN model
            model = Sequential()
            # Defining the first layer of the model
            model.add(Dense(units=5, input_dim=X_train.shape[1], kernel_initializer='normal', activation='relu'))
            
            # Defining the Second layer of the model
            model.add(Dense(units=5, kernel_initializer='normal', activation='relu'))
            # The output neuron is a single fully connected node 
            # Since we will be predicting a single number
            model.add(Dense(1, kernel_initializer='normal'))
            # Compiling the model
            model.compile(loss='mean_squared_error', optimizer='adam')
            
            # Fitting the ANN to the Training set
            model.fit(X_train, y_train ,batch_size = batch_size_trial, epochs = epochs_trial, verbose=0)
            MAPE = np.mean(100 * (np.abs(y_test-model.predict(X_test))/y_test))
            
            # printing the results of the current iteration
            print(TrialNumber, 'Parameters:','batch_size:', batch_size_trial,'-', 'epochs:',epochs_trial, 'Accuracy:', 100-MAPE)
            SearchResultsData=SearchResultsData.append(pd.DataFrame(data=[[TrialNumber, str(batch_size_trial)+'-'+str(epochs_trial), 100-MAPE]],
                                                                    columns=['TrialNumber', 'Parameters', 'Accuracy']))
    return(SearchResultsData)


# In[25]:


# Calling the function
ResultsData=FunctionFindBestParams(x_train, y_train, x_test, y_test)


# In[31]:


#Examine single model fit
model.fit(x_train, y_train ,batch_size = 5, epochs = 10, verbose=0)
#MAPE = np.mean(100 * (np.abs(y_test-model.predict(X_test))/y_test))


# In[34]:


y_train


# In[142]:


#View parameter accuracy
get_ipython().run_line_magic('matplotlib', 'inline')
ResultsData.plot(x='Parameters', y='Accuracy', figsize=(15,4), kind='line')


# In[112]:


model = Sequential()
# Defining the first layer of the model
model.add(Dense(units=5, input_dim=x_train.shape[1], kernel_initializer='normal', activation='relu'))
            
# Defining the Second layer of the model
model.add(Dense(units=5, kernel_initializer='normal', activation='relu'))
# The output neuron is a single fully connected node 
# Since we will be predicting a single number
model.add(Dense(1, kernel_initializer='normal'))
# Compiling the model
model.compile(loss='mean_squared_error', optimizer='adam')
            


# In[143]:


# Fitting the ANN to the Training set
model.fit(x_train, y_train, batch_size = 20, epochs = 100, verbose=0)


# In[144]:


# Generating Predictions on testing data
Predictions=model.predict(x_test)


# In[164]:


### Sandardization of data ###
from sklearn.preprocessing import StandardScaler
PredictorScaler=StandardScaler()
TargetVarScaler=StandardScaler()
# Storing the fit object for later reference
PredictorScalerFit=PredictorScaler.fit(x)
TargetVarScalerFit=TargetVarScaler.fit(x)
# Generating the standardized values of X and y
x=PredictorScalerFit.transform(x)
y=TargetVarScalerFit.transform(y)


# In[2]:


# Scaling the predicted Price data back to original price scale
#Predictions=TargetVarScalerFit.inverse_transform(Predictions)


# In[153]:


# Scaling the y_test Price data back to original price scale
#y_test_orig=TargetVarScalerFit.inverse_transform(y_test)
# Scaling the test data back to original scale
#Test_Data=PredictorScalerFit.inverse_transform(X_test)
TestingData=pd.DataFrame(data=x_test, columns=x_train.columns)
TestingData['Price']=y_test
TestingData['PredictedPrice']=Predictions
TestingData.head()


# In[ ]:


# printing the results of the current iteration
print(TrialNumber, 'Parameters:','batch_size:', 5,'-', 'epochs:',5, 'Accuracy:', 100-MAPE)
SearchResultsData=SearchResultsData.append(pd.DataFrame(data=[[TrialNumber, str(5)+'-'+str(5), 100-MAPE]],
                                                                    columns=['TrialNumber', 'Parameters', 'Accuracy']))

