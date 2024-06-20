
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
import scipy.stats as stats
from sklearn.calibration import CalibratedClassifierCV
from sklearn.calibration import calibration_curve
pd.options.mode.chained_assignment = None  # default='warn'

def get_pred_probs(tmp_f, mydf, target_y, fold_id_lst, my_params, col_name):
    eid_lst, region_lst = [], []
    y_test_lst, y_pred_lst = [], []
    for fold_id in fold_id_lst:
        train_idx = mydf['Region_code'].index[mydf['Region_code'] != fold_id]
        test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
        X_train, X_test = mydf.iloc[train_idx][tmp_f], mydf.iloc[test_idx][tmp_f]
        y_train, y_test = mydf[target_y].iloc[train_idx], mydf[target_y].iloc[test_idx]
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2020)
        my_lgb.set_params(**my_params)
        calibrate = CalibratedClassifierCV(my_lgb, method='isotonic', cv=5)
        calibrate.fit(X_train, y_train)
        y_pred_prob = calibrate.predict_proba(X_test)[:, 1].tolist()
        y_pred_lst += y_pred_prob
        y_test_lst += mydf[target_y].iloc[test_idx].tolist()
        eid_lst += mydf.eid.iloc[test_idx].tolist()
        region_lst += mydf.Region_code.iloc[test_idx].tolist()
    myout_df = pd.DataFrame([eid_lst, region_lst, y_test_lst, y_pred_lst]).T
    myout_df.columns = ['eid', 'Region_code', target_y, 'y_pred_' + col_name]
    myout_df[['eid', 'Region_code']] = myout_df[['eid', 'Region_code']].astype('int')
    return myout_df

my_params = {'n_estimators': 300,
             'max_depth': 10,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}

target = 'ACD'
food = 'FoodGroupsResidual'
dpath = '/Volumes/JasonWork/Projects/OtherProjects/ChenSijia/Diet/'
outfile = dpath + 'Results/' + food + '/' + target + '/PredProbs.csv'

food_df = pd.read_csv(dpath + 'Data/' + food + '.csv')
target_df = pd.read_csv(dpath + 'Data/Dementia.csv', usecols = ['eid', target, 'time'])
reg_df = pd.read_csv(dpath + 'Data/Eid_info_data.csv', usecols = ['eid', 'Region_code'])

mydf = pd.merge(target_df, reg_df, how = 'inner', on = ['eid'])
mydf = pd.merge(mydf, food_df, how = 'inner', on = ['eid'])
fold_id_lst = [i for i in range(10)]

my_f_df = pd.read_csv(dpath + 'Results/' + food + '/' + target + '/FoodImportance.csv')
my_f_df.sort_values(by = 'TotalGain_cv', ascending=False, inplace = True)
food_TG_lst = my_f_df.Food.tolist()[:8]

pred_TG_df = get_pred_probs(food_TG_lst, mydf, target, fold_id_lst, my_params, 'TG_8foods')

myout_df = pd.merge(target_df, pred_TG_df, how = 'inner', on = ['eid'])

myout_df.to_csv(outfile, index = False)


