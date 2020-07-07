import math
from random import uniform
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
%matplotlib inline
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LinearRegression
import sklearn
from sklearn import linear_model
from sklearn.neighbors import KNeighborsRegressor

from sklearn.neural_network import MLPRegressor
from sklearn.datasets import make_regression
from sklearn.model_selection import GridSearchCV, learning_curve
from sklearn import svm, datasets

import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import scipy.stats
from scipy.stats import multivariate_normal
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import GridSearchCV
 
from sklearn.svm import SVC, LinearSVC

from sklearn import svm

import scipy.stats

from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix

from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import mean_squared_error, r2_score


# ------------------------------------------------------------------------
# ------------------------DATA ADJUSTMENT----------------------------------
# ------------------------------------------------------------------------

# drop column, because it's dosent metter for learning
def removeBadColumns(df: pd.DataFrame) -> None:
    del df['Case #']
    del df['Cow']
    del df['Sire']
    del df['Birth year']
    del df['Dam']
    del df['inbreeding']
	
def boxplotOfData(df: pd.DataFrame, columnName: str) -> None:
    plt.figure(figsize=(15,6))
    
    sns.boxplot(x=columnName, data=df)
    plt.ticklabel_format(style='plain', axis='x')
    plt.title('Boxplot of roh-cow')
    plt.show()
    
def removeRowsNotInRange(df: pd.DataFrame, columnName, startRange, endRange) -> pd.DataFrame:
    df = df[(df[columnName] >= startRange) & (df[columnName] <= endRange)]
    return df
	
def correlationMatrix(df: pd.DataFrame) -> None:
    # Compute the correlation matrix
    d = df[['ROH-Cow' , 'ROH-Sire', 'ROH-DAM' , 'IBS', 'IBD', 'HA', 'expect-roh'  ]]
    corr = d.corr()

    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, mask=mask,
                square=True, linewidths=.5, annot=True, cmap=cmap)
    plt.yticks(rotation=0)
    plt.title('Correlation Matrix of all Numerical Variables')
    plt.show()

def histogramGragh(df: pd.DataFrame, columnName: str) -> None:
 
	# Set the size of the plot
	plt.figure(figsize=(10, 6))
	
	# Plot the data and configure the settings
	sns.distplot(df['ROH-Cow'])
	plt.title('Histogram of ROH-Cow')
	plt.show()


# ------------------------------------------------------------------------
# ----------------------------"MAIN"--------------------------------------
# ------------------------------------------------------------------------

df = pd.read_excel('Measures.xlsx')
removeBadColumns(df)
#df = shuffle(df)
boxplotOfData(df, 'ROH-Cow')

df = removeRowsNotInRange(df, 'ROH-Cow', startRange=600, endRange=950)

# print(len(df))
# print(df.head())
# print(df.describe())
histogramGragh(df, 'ROH-Cow')

correlationMatrix(df)

data_target = df['ROH-Cow']
df = df.drop('ROH-Cow', 1)



df = df.drop('ROH-DAM' ,1 )
df = df.drop('ROH-Sire' ,1 )
#df = df.drop('IBS' , 1)
#df = df.drop('IBD' , 1)
#df = df.drop('HA' , 1)
#df = df.drop('expect-roh' , 1)

# ------------------------------------------------------------------------
# --------------------split to train & test-------------------------------
# ------------------------------------------------------------------------

# Split data into training and testing set with 70% of the data going into training
#x_train, x_test, y_train , y_test = train_test_split(df.values , data_target.values, test_size=0.265, random_state=0)

# Split data into training and testing 
x_train = df.values[:1103]
x_test = df.values[1103:]
y_train = data_target.values[:1103]
y_test = data_target.values[1103:]

# normalize x_train and x_test matrix
means = np.mean(x_train, axis=0)
stds = np.std(x_train, axis=0)
x_train = (x_train - means) / stds
x_test = (x_test - means) / stds

x_train = np.array(x_train)
y_train = np.array(y_train)

# ------------------------------------------------------------------------
# ----------------------Linear Regression---------------------------------
# ------------------------------------------------------------------------

print("Linear Regression")
regr = linear_model.LinearRegression()
parameters = {'fit_intercept':[True,False], 'normalize':[True,False], 'copy_X':[True, False] , 'n_jobs' : [1,-1]}
grid_cv = GridSearchCV(regr,parameters, cv=10)
grid_cv.fit(x_train, y_train)
print("tuned hpyerparameters :(best parameters) ",grid_cv.best_params_)

y_pred_train = grid_cv.predict(x_train)
print('Mean Squared Error - train:', mean_squared_error(y_train, y_pred_train))  
print('Root Mean Squared Error - train:', np.sqrt(mean_squared_error(y_train, y_pred_train)))

# Make predictions using the testing set
y_pred = grid_cv.predict(x_test)


# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f'
      % r2_score(y_test, y_pred))


print('Mean Squared Error:', mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(mean_squared_error(y_test, y_pred)))

# Compute 5-fold cross-validation scores: cv_scores
cv_scores_linreg = cross_val_score(grid_cv, x_train, y_train, cv=5)
print("Average 5-Fold CV Score: {}".format(np.mean(cv_scores_linreg)))
# Print the 5-fold cross-validation scores
print(cv_scores_linreg)

train_sizes, train_scores, test_scores = learning_curve(grid_cv, x_train, y_train , n_jobs=5)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)


plt.fill_between(train_sizes, train_scores_mean - train_scores_std,train_scores_mean + train_scores_std, alpha=0.1, color="#ff9124", )
plt.fill_between(train_sizes, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.1, color="#2492ff")
plt.plot(train_sizes, train_scores_mean, 'o-', color="#ff9124", label="Training score")
plt.plot(train_sizes, test_scores_mean, 'o-', color="#2492ff", label="validation score")
plt.title("Linear Regression Learning Curve", fontsize=14)
plt.xlabel('Training size (m)')
plt.ylabel('Score')
plt.grid(True)
plt.legend(loc="best")
# ------------------------------------------------------------------------
# ----------------------Support Vector Regression---------------------------------
# ------------------------------------------------------------------------

parameters = {'kernel': ('linear', 'rbf','poly' , 'sigmoid'), 'C':[1,5,1.5, 10],'gamma': [1e-7, 1e-4],'epsilon':[0.1,0.2,0.5,0.3]}
svr = svm.SVR()
clf = GridSearchCV(svr, parameters)
clf.fit(x_train,y_train)
print(clf.best_params_)


y_pred = clf.predict(x_test)
# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f' % r2_score(y_test, y_pred))

y_pred_train = clf.predict(x_train)
print('Mean Squared Error - train:', mean_squared_error(y_train, y_pred_train))  
print('Root Mean Squared Error - train:', np.sqrt(mean_squared_error(y_train, y_pred_train)))

# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f'
      % r2_score(y_test, y_pred))


print('Mean Squared Error:', mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(mean_squared_error(y_test, y_pred)))

# Compute 5-fold cross-validation scores: cv_scores
cv_scores_linreg = cross_val_score(clf, x_train, y_train, cv=5)
print("Average 5-Fold CV Score: {}".format(np.mean(cv_scores_linreg)))
# Print the 5-fold cross-validation scores
print(cv_scores_linreg)


#plt.ylim(0.5, 1.01)
train_sizes, train_scores, test_scores = learning_curve(clf, x_train, y_train)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)


plt.fill_between(train_sizes, train_scores_mean - train_scores_std,train_scores_mean + train_scores_std, alpha=0.1, color="#ff9124", )
plt.fill_between(train_sizes, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.1, color="#2492ff")
plt.plot(train_sizes, train_scores_mean, 'o-', color="#ff9124", label="Training score")
plt.plot(train_sizes, test_scores_mean, 'o-', color="#2492ff", label="validation score")
plt.title("svm.svr Learning Curve", fontsize=14)
plt.xlabel('Training size (m)')
plt.ylabel('Score')
plt.grid(True)
plt.legend(loc="best")

# ------------------------------------------------------------------------
# ----------------------XGB Regressor---------------------------------
# ------------------------------------------------------------------------

xgb1 = XGBRegressor()
parameters = {'nthread':[4], #when use hyperthread, xgboost may become slower
              'objective':['reg:linear'],
              'learning_rate': [0.3, 0.5, 0.7], #so called `eta` value
              'max_depth': [1,2,3,4,5],   
              'colsample_bytree': [0.7,0.3],
              'n_estimators': [500,10], 
              'alpha': [5]
              }

xgb_grid = GridSearchCV(xgb1,parameters,cv = 2,n_jobs = 5,verbose=True)

xgb_grid.fit(x_train,y_train)


print(xgb_grid.best_params_)
y_pred = xgb_grid.predict(x_test)
#rmse = np.sqrt(mean_squared_error(y_test, preds))
#print("RMSE: %f" % (rmse))


y_pred_train = xgb_grid.predict(x_train)
print('Mean Squared Error - train:', mean_squared_error(y_train, y_pred_train))  
print('Root Mean Squared Error - train:', np.sqrt(mean_squared_error(y_train, y_pred_train)))

# Make predictions using the testing set
y_pred = xgb_grid.predict(x_test)

# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f'
      % r2_score(y_test, y_pred))


print('Mean Squared Error:', mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(mean_squared_error(y_test, y_pred)))

# Compute 5-fold cross-validation scores: cv_scores
cv_scores_linreg = cross_val_score(xgb_grid, x_train, y_train, cv=5)
print("Average 5-Fold CV Score: {}".format(np.mean(cv_scores_linreg)))
# Print the 5-fold cross-validation scores
print(cv_scores_linreg)

#plt.ylim(0.5, 1.01)
train_sizes, train_scores, test_scores = learning_curve(xgb_grid, x_train, y_train)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)


plt.fill_between(train_sizes, train_scores_mean - train_scores_std,train_scores_mean + train_scores_std, alpha=0.1, color="#ff9124", )
plt.fill_between(train_sizes, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.1, color="#2492ff")
plt.plot(train_sizes, train_scores_mean, 'o-', color="#ff9124", label="Training score")
plt.plot(train_sizes, test_scores_mean, 'o-', color="#2492ff", label="validation score")
plt.title("XGB Regressor Learning Curve", fontsize=14)
plt.xlabel('Training size (m)')
plt.ylabel('Score')
plt.grid(True)
plt.legend(loc="best")

# ------------------------------------------------------------------------
# ------------------------Random Forest-----------------------------------
# ------------------------------------------------------------------------
from sklearn.model_selection import GridSearchCV
# Create the parameter grid based on the results of random search 
param_grid = {
    'bootstrap': [True],
    'max_depth': [20,80, 90, 100, 110],
    'max_features': [2, 3 , 'sqrt'],
    'min_samples_leaf': [2, 3, 4, 5],
    'min_samples_split': [5, 8, 10, 12],
    'n_estimators': [100, 200, 300, 1000,1200]
}
# Create a based model
rf = RandomForestRegressor()
# Instantiate the grid search model
grid_search = GridSearchCV(estimator = rf, param_grid = param_grid, 
                          cv = 3, n_jobs = -1, verbose = 2)


grid_search.fit(x_train,y_train)
# Apply The Full Featured Classifier To The Test Data
y_pred = grid_search.predict(x_test)

# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f' % r2_score(y_test, y_pred))


y_pred_train = grid_search.predict(x_train)
print('Mean Squared Error - train:', mean_squared_error(y_train, y_pred_train))  
print('Root Mean Squared Error - train:', np.sqrt(mean_squared_error(y_train, y_pred_train)))

print(grid_search.best_params_)

# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f'
      % r2_score(y_test, y_pred))


print('Mean Squared Error:', mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(mean_squared_error(y_test, y_pred)))

# Compute 5-fold cross-validation scores: cv_scores
#cv_scores_linreg = cross_val_score(grid_search, x_train, y_train, cv=5)
#print("Average 5-Fold CV Score: {}".format(np.mean(cv_scores_linreg)))
# Print the 5-fold cross-validation scores
#print(cv_scores_linreg)

from sklearn.model_selection import GridSearchCV, learning_curve

import matplotlib.pyplot as plt
#plt.ylim(0.5, 1.01)
train_sizes, train_scores, test_scores = learning_curve(grid_search, x_train, y_train , n_jobs=5)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)


plt.fill_between(train_sizes, train_scores_mean - train_scores_std,train_scores_mean + train_scores_std, alpha=0.5, color="#ff9124", )
plt.fill_between(train_sizes, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.5, color="#2492ff")
plt.plot(train_sizes, train_scores_mean, 'o-', color="#ff9124", label="Training score")
plt.plot(train_sizes, test_scores_mean, 'o-', color="#2492ff", label="validation score")
plt.title("Random Forest Learning Curve", fontsize=14)
plt.xlabel('Training size (m)')
plt.ylabel('Score')
plt.grid(True)
plt.legend(loc="best")

   
    
# ------------------------------------------------------------------------
# ----------------------K-Nearest Neighbors---------------------------------
# ------------------------------------------------------------------------

params = {'n_neighbors':[1,2,3,4,5,6,7,8,9 ,10,11,12,13,14, 15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]}
knn = KNeighborsRegressor()
model = GridSearchCV(knn, params, cv=5)
model.fit(x_train,y_train)
print(model.best_params_)
y_pred_train = model.predict(x_train)
print('Mean Squared Error - train:', mean_squared_error(y_train, y_pred_train))  
print('Root Mean Squared Error - train:', np.sqrt(mean_squared_error(y_train, y_pred_train)))
# Make predictions using the testing set
y_pred = model.predict(x_test)
# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f'
      % r2_score(y_test, y_pred))
print('Mean Squared Error:', mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(mean_squared_error(y_test, y_pred)))

train_sizes, train_scores, test_scores = learning_curve(model, x_train, y_train , n_jobs=5)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)


plt.fill_between(train_sizes, train_scores_mean - train_scores_std,train_scores_mean + train_scores_std, alpha=0.1, color="#ff9124", )
plt.fill_between(train_sizes, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.1, color="#2492ff")
plt.plot(train_sizes, train_scores_mean, 'o-', color="#ff9124", label="Training score")
plt.plot(train_sizes, test_scores_mean, 'o-', color="#2492ff", label="validation score")
plt.title("KNN Learning Curve", fontsize=14)
plt.xlabel('Training size (m)')
plt.ylabel('Score')
plt.grid(True)
plt.legend(loc="best")



# ------------------------------------------------------------------------
# ----------------------Neural Networks---------------------------------
# ------------------------------------------------------------------------

mlp = MLPRegressor()
parameter_space = {'hidden_layer_sizes': [1,2,3,4,5,6,7,8,9,10,20,30,40],
'activation': ['relu'],
'solver':['lbfgs'], 'alpha':[0.0001],
'batch_size':['auto'], 'learning_rate':['constant'],
'learning_rate_init':[0.001], 'max_iter':[500]}
clf = GridSearchCV(mlp, parameter_space, n_jobs=-1, cv=3)
clf.fit(x_train, y_train)
print("tuned hpyerparameters :(best parameters) ",clf.best_params_)

y_pred_train = clf.predict(x_train)
print('Mean Squared Error - train:', mean_squared_error(y_train, y_pred_train))  
print('Root Mean Squared Error - train:', np.sqrt(mean_squared_error(y_train, y_pred_train)))
# Make predictions using the testing set
y_pred = clf.predict(x_test)
# The coefficient of determination: 1 is perfect prediction
print('Coefficient of determination: %.2f'
      % r2_score(y_test, y_pred))
print('Mean Squared Error:', mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(mean_squared_error(y_test, y_pred)))

train_sizes, train_scores, test_scores = learning_curve(clf, x_train, y_train , n_jobs=5)
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)


plt.fill_between(train_sizes, train_scores_mean - train_scores_std,train_scores_mean + train_scores_std, alpha=0.1, color="#ff9124", )
plt.fill_between(train_sizes, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.1, color="#2492ff")
plt.plot(train_sizes, train_scores_mean, 'o-', color="#ff9124", label="Training score")
plt.plot(train_sizes, test_scores_mean, 'o-', color="#2492ff", label="validation score")
plt.title("Neural Networks Learning Curve", fontsize=14)
plt.xlabel('Training size (m)')
plt.ylabel('Score')
plt.grid(True)
plt.legend(loc="best")

#####################################################################################

df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df

df1 = df.head(25)
df1.plot(kind='bar',figsize=(16,10))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()


# ------------------------------------------------------------------------
# -------------------------------plot_params_effect------------------------
# ------------------------------------------------------------------------


def plot_params_effect(lr_cv_results, svc_cv_results):
    f, (ax2, ax3) = plt.subplots(1, 2, figsize=(30, 10))
   
    lr_params = [str(value) for value in lr_cv_results['params']]
    ax2.plot(range(len(lr_params)), lr_cv_results['mean_test_score'], 'o-', color='blue')
    ax2.set_xticklabels(lr_params, rotation='vertical', fontsize=10)
    ax2.set_xlabel('Parameters')
    ax2.set_ylabel('Score')

    svc_params = [str(value) for value in svc_cv_results['params']]
    ax3.plot(range(len(svc_params)), svc_cv_results['mean_test_score'], 'o-', color='blue')
    ax3.set_xticklabels(svc_params,  rotation='vertical', fontsize=10)
    ax3.set_xlabel('Parameters')
    ax3.set_ylabel('Score')
    
    plt.show()
    
plot_params_effect(xgb_grid.cv_results_, clf.cv_results_)

