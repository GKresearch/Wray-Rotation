```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Documentation:
Final Updated Machine Learning Harness script for a rotation project in the Wray Lab
author=GKennedy
previous_author=jamesonblount
    (see below for details) 
Date=2022 December 12th

Note: this was developed from another script authored by Jameson Blount, see jb621-star/Wray-Rotation/ML-Harness-K562.py for the original
"""


# Importing required packages
#Importing basic packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#Importing sklearn modules
import rfpimp
from sklearn.linear_model import RANSACRegressor
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
from sklearn.svm import SVC, SVR
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import make_column_transformer
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
```


```python
#Loading the data and checking for missing values
dataset=pd.read_csv('C:/Users/ictinike/Documents/WrayLab/raw_data/x_0011_df_phyloP.csv')
dataset.isnull().sum()

datasetv2 = dataset.dropna(axis=1)
datasetv2.isnull().sum()
# Checking the data set for any NULL values is very essential, as MLAs can not 
# handle NULL values. We have to either eliminate the records with NULL values 
# or replace them with the mean/median of the other values. we can see each of 
# the variables are printed with number of null values. This data set has no null 
# values so all are zero here.
```




    seqnames                  0
    start                     0
    end                       0
    width                     0
    strand                    0
                             ..
    gene.y                    0
    dTSS                      0
    PhastCons                 0
    PhyloP_primates_score     0
    PhyloP_placental_score    0
    Length: 90, dtype: int64




```python
test=datasetv2#[1:50,::]
```


```python
#Read in the standard dataframe with TADs 
df2=pd.read_csv('C:/Users/Ictinike/Documents/WrayLab/raw_data/OCRs_inTADs.csv')
df2.isnull().sum()
```




    seqnames                  0
    start                     0
    end                       0
    width                     0
    strand                    0
                             ..
    PhyloP_placental_score    0
    TAD_chr                   0
    TAD_start                 0
    TAD_end                   0
    TAD_ID                    0
    Length: 145, dtype: int64




```python
#Read in the Gene-Element Tad Confirmed, and Removal of off TAD gene-element pairs dataframe
df3=pd.read_csv('C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/x0011_Gene-Element_TAD_Domain_Check_No-OFF-TAD_2022-12-13.csv')
df3.isnull().sum()
```




    seqnames         0
    start            0
    end              0
    width            0
    strand           0
                    ..
    TAD_ID           0
    CHRPOS           0
    element_check    0
    gene_check       0
    Relationship     0
    Length: 149, dtype: int64




```python
#Read in the Gene_Element TAD Confirmed Dataframe
df4=pd.read_csv('C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/x0011_Gene-Element_TAD_Domain_Check_2022-12-13.csv')
df4.isnull().sum()
```




    seqnames         0
    start            0
    end              0
    width            0
    strand           0
                    ..
    TAD_ID           0
    CHRPOS           0
    element_check    0
    gene_check       0
    Relationship     0
    Length: 149, dtype: int64




```python
#Read in the Gene-Element TAD Confirmed Promoter-Enhancer pair existing dataframe
df5=pd.read_csv('C:/Users/Ictinike/Documents/WrayLab/raw_data/pipeline_data/x0011_Gene-Element_TAD_Domain_Check_Confirmed_Promoter-Enhancer_2022-12-14.csv')
df5.isnull().sum()
```




    chromHMM_cat_longest    0
    start                   0
    end                     0
    seqnames                0
    gene_id                 0
                           ..
    TAD_chr.y               0
    TAD_start.y             0
    TAD_end.y               0
    TAD_ID.y                0
    CHRPOS.y                0
    Length: 157, dtype: int64




```python
#Assign whichever dataframe you are working with to be tested
test=df5#.iloc#[1:5000,::]
len(test)
```




    26721




```python
# with wgCERES_score_nosig as the response vector,
x = test[["DHS_prop_repeat", "DHS_prop_GC", "DHS_length", "n_SNV_Zhou_per_bp", 
                    "distanceToTSS", "zeta.human", "zeta.chimp", "PP_con", "PP_acc", 
                    "PhastCons",
                    "chromHMM_cat_longest", 
                    "annotation", "PhyloP_primates_score","Relationship"]]
y = test["wgCERES_score_nosig"]
x = pd.get_dummies(x, columns = ['chromHMM_cat_longest','annotation',"Relationship"])
```


```python
random_state=None
```


```python
#Split the Dataset
x_train, x_test, y_train, y_test=train_test_split(x,y,test_size=0.2, random_state=0)
```


```python
len(x_test)
```




    5345




```python
#Set up intended models to run
lm  = LinearRegression()
gbm = GradientBoostingRegressor()
#lr = LogisticRegression()
rf = RandomForestRegressor()
sv = SVR()
```


```python
ransac = RANSACRegressor(LinearRegression(),
                        max_trials=20,
                        min_samples=50,
                        residual_threshold=5.0,
                        random_state=0)
```


```python
#Fit Training Data
ransac.fit(x_train, y_train)
lm.fit(x_train, y_train)
rf.fit(x_train, y_train)
gbm.fit(x_train, y_train)
sv.fit(x_train, y_train)
```




    SVR()




```python
#Predict Test data
rf_pred=rf.predict(x_test)
ran_pred=ransac.predict(x_test)
sv_pred=sv.predict(x_test)
lm_pred=lm.predict(x_test)
gbm_pred=gbm.predict(x_test)
```


```python
#Calculate Mean absolute and squared errors

ran_mae = mean_absolute_error(ran_pred, y_test)
ran_rmse = np.sqrt(mean_squared_error(ran_pred, y_test))
print("RAN MAE: {:.2f}".format(round(ran_mae, 2)))
print("RAN RMSE: {:.2f}".format(round(ran_rmse, 2)))
lm_mae = mean_absolute_error(lm_pred, y_test)
lm_rmse = np.sqrt(mean_squared_error(lm_pred, y_test))
print("LM MAE: {:.2f}".format(round(lm_mae, 2)))
print("LM RMSE: {:.2f}".format(round(lm_rmse, 2)))
sv_mae = mean_absolute_error(sv_pred, y_test)
sv_rmse = np.sqrt(mean_squared_error(sv_pred, y_test))
print("SV MAE: {:.2f}".format(round(sv_mae, 2)))
print("SV RMSE: {:.2f}".format(round(sv_rmse, 2)))
gbm_mae = mean_absolute_error(gbm_pred, y_test)
gbm_rmse = np.sqrt(mean_squared_error(gbm_pred, y_test))
print("GBM MAE: {:.2f}".format(round(gbm_mae, 2)))
print("GBM RMSE: {:.2f}".format(round(gbm_rmse, 2)))
rf_mae = mean_absolute_error(rf_pred, y_test)
rf_rmse = np.sqrt(mean_squared_error(rf_pred, y_test))
print("RF MAE: {:.2f}".format(round(rf_mae, 2)))
print("RF RMSE: {:.2f}".format(round(rf_rmse, 2)))
```

    RAN MAE: 1.29
    RAN RMSE: 1.92
    LM MAE: 1.26
    LM RMSE: 1.91
    SV MAE: 1.27
    SV RMSE: 1.94
    GBM MAE: 1.27
    GBM RMSE: 1.91
    RF MAE: 1.31
    RF RMSE: 1.97
    


```python
from tabulate import tabulate
```


```python
#Tabulate the data
table = [['Linear Regression', 'MAE', round(lm_mae,2)], 
         ['-', 'RMSE', round(lm_rmse,2)], 
         ['Gradient Boosting Machine', "MAE", round(gbm_mae,2)], 
         ['-', 'RMSE', round(gbm_rmse,2)], 
         ['Support Vector Regression', "MAE", round(sv_mae,2)], 
         ['-', 'RMSE', round(sv_rmse,2)], 
         ['RANSAC Regression', "MAE", round(ran_mae,2)], 
         ['-', 'RMSE', round(ran_rmse,2)], 
         ['Random Forest Regression', "MAE", round(rf_mae,2)],
         ['-', 'RMSE', round(rf_rmse,2)]]
```


```python
print(tabulate(table))
```

    -------------------------  ----  ----
    Linear Regression          MAE   1.26
    -                          RMSE  1.91
    Gradient Boosting Machine  MAE   1.27
    -                          RMSE  1.91
    Support Vector Regression  MAE   1.27
    -                          RMSE  1.94
    RANSAC Regression          MAE   1.29
    -                          RMSE  1.92
    Random Forest Regression   MAE   1.31
    -                          RMSE  1.97
    -------------------------  ----  ----
    


```python
print(tabulate(table, tablefmt='fancy_grid'))
```

    ╒═══════════════════════════╤══════╤══════╕
    │ Linear Regression         │ MAE  │ 1.26 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.91 │
    ├───────────────────────────┼──────┼──────┤
    │ Gradient Boosting Machine │ MAE  │ 1.27 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.91 │
    ├───────────────────────────┼──────┼──────┤
    │ Support Vector Regression │ MAE  │ 1.27 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.94 │
    ├───────────────────────────┼──────┼──────┤
    │ RANSAC Regression         │ MAE  │ 1.29 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.92 │
    ├───────────────────────────┼──────┼──────┤
    │ Random Forest Regression  │ MAE  │ 1.31 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.97 │
    ╘═══════════════════════════╧══════╧══════╛
    


```python
print(tabulate(table, tablefmt='fancy_grid'))
```

    ╒═══════════════════════════╤══════╤══════╕
    │ Linear Regression         │ MAE  │ 1.26 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.91 │
    ├───────────────────────────┼──────┼──────┤
    │ Gradient Boosting Machine │ MAE  │ 1.27 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.91 │
    ├───────────────────────────┼──────┼──────┤
    │ Support Vector Regression │ MAE  │ 1.27 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.94 │
    ├───────────────────────────┼──────┼──────┤
    │ RANSAC Regression         │ MAE  │ 1.29 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.92 │
    ├───────────────────────────┼──────┼──────┤
    │ Random Forest Regression  │ MAE  │ 1.31 │
    ├───────────────────────────┼──────┼──────┤
    │ -                         │ RMSE │ 1.97 │
    ╘═══════════════════════════╧══════╧══════╛
    


```python

```
