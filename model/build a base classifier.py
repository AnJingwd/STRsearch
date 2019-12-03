#!/usr/bin/env python
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: build_a_base_classifier.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch

'''
used k-fold cross validation and selected the best model with the best parameters & hyper-parameters
Optimal model score (accuracy):0.9612698412698413
'''

import argparse
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import GridSearchCV
from sklearn.externals import joblib

parser = argparse.ArgumentParser(description='Build a base classifier with quality control matrix')
parser.add_argument('--qc_matrix', help='(must) assign quality control matrix file from STRsearch',required=True)
args = parser.parse_args()


def preprocessing_data(file,features):
    df = pd.read_csv(file,header=0)
    df["Stutter_ratio"] = df["Supp_reads2"] / df["Supp_reads1"]
    # split data into X and y
    X = df[features+["Stutter_ratio"]]
    Y = df["Label"]
    # set missing values to 0
    X_new= X.fillna(0)
    return X_new,Y

features_name = ["Total_bases","Num_reads","Q30","Dis1_min_5","Dis1_max_5","Dis1_min_3","Dis1_max_3","Dis1_mean_5", \
            "Dis1_mean_3","Dis2_min_5","Dis2_max_5","Dis2_min_3","Dis2_max_3","Dis2_mean_5","Dis2_mean_3","Supp_reads1","Supp_reads2"]
X,y = preprocessing_data(args.qc_matrix,features_name)


#############    Manually select the best superparameter step by step   ################

#cv_params = {'n_estimators': [300, 310, 320, 330, 340]}
#cv_params = {'max_depth': [3, 4, 5, 6, 7, 8, 9, 10], 'min_child_weight': [1, 2, 3, 4, 5, 6]}
#cv_params = {'gamma': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7]}
#cv_params = {'subsample': [0.6, 0.7, 0.8, 0.9], 'colsample_bytree': [0.6, 0.7, 0.8, 0.9]}
#cv_params = {'reg_alpha': [0.05, 0.1, 1, 2, 3], 'reg_lambda': [0.05, 0.1, 1, 2, 3]}
cv_params = {'learning_rate': [0.01, 0.05, 0.07, 0.1, 0.2,0.25]}

## original_params = {'learning_rate': 0.1, 'n_estimators': 500, 'max_depth': 5, 'min_child_weight': 1, 'seed': 0,
##                    'subsample': 0.8, 'colsample_bytree': 0.8, 'gamma': 0, 'reg_alpha': 0, 'reg_lambda': 1}
other_params = {'learning_rate': 0.1, 'n_estimators': 320, 'max_depth': 7, 'min_child_weight': 1, 'seed': 0,
                    'subsample': 0.9, 'colsample_bytree': 0.7, 'gamma': 0.6, 'reg_alpha': 0.1, 'reg_lambda': 0.1}

kFold = 5
model = xgb.XGBClassifier(**other_params)
optimized_GBM = GridSearchCV(estimator=model, param_grid=cv_params, scoring='accuracy', cv=kFold, verbose=1, n_jobs=4)
optimized_GBM.fit(X, y)
evalute_result = optimized_GBM.cv_results_
print('results of each epoch:{0}'.format(evalute_result))
print('The best value of the hyper-parameters：{0}'.format(optimized_GBM.best_params_))
print('Optimal model score (accuracy):{0}'.format(optimized_GBM.best_score_))

################   save best model   #####################
joblib.dump(optimized_GBM.best_estimator_, 'Ion_S5.pkl', compress = 1)
