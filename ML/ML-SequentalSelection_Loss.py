
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, log_loss
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
import scipy.stats as stats
pd.options.mode.chained_assignment = None  # default='warn'

target = 'ACD'
food = 'FoodGroupsResidual'
ImpMethod = 'TotalGain'

dpath = '/Volumes/JasonWork/Projects/OtherProjects/ChenSijia/Diet/'
outfile = dpath + 'Results/' + food + '/' + target + '/AccLOSS_' + ImpMethod + '.csv'

food_df = pd.read_csv(dpath + 'Data/' + food + '.csv')
food_f_lst = food_df.columns.tolist()[1:]
target_df = pd.read_csv(dpath + 'Data/Dementia.csv', usecols = ['eid', target, 'time'])
reg_df = pd.read_csv(dpath + 'Data/Eid_info_data.csv', usecols = ['eid', 'Region_code'])

mydf = pd.merge(target_df, reg_df, how = 'inner', on = ['eid'])
mydf = pd.merge(mydf, food_df, how = 'inner', on = ['eid'])

my_f_df = pd.read_csv(dpath + 'Results/' + food + '/' + target + '/FoodImportance.csv')
my_f_df.sort_values(by = ImpMethod + '_cv', ascending=False, inplace = True)
food_f_lst = my_f_df.Food.tolist()

fold_id_lst = [i for i in range(10)]

my_params = {'n_estimators': 300,
             'max_depth': 10,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}

y_test_full = np.zeros(shape = [1,1])
for fold_id in fold_id_lst:
    test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
    y_test_full = np.concatenate([y_test_full, np.expand_dims(mydf.iloc[test_idx][target], -1)])

y_pred_full_prev = y_test_full
tmp_f, loss_cv_lst= [], []

for f in food_f_lst:
    tmp_f.append(f)
    my_X = mydf[tmp_f]
    loss_cv = []
    y_pred_full = np.zeros(shape = [1,1])
    for fold_id in fold_id_lst:
        train_idx = mydf['Region_code'].index[mydf['Region_code'] != fold_id]
        test_idx = mydf['Region_code'].index[mydf['Region_code'] == fold_id]
        X_train, X_test = mydf.iloc[train_idx][tmp_f], mydf.iloc[test_idx][tmp_f]
        y_train, y_test = mydf.iloc[train_idx][target], mydf.iloc[test_idx][target]
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True,  n_jobs=4, verbosity=-1, seed=2020)
        my_lgb.set_params(**my_params)
        my_lgb.fit(X_train, y_train)
        y_pred_prob = my_lgb.predict_proba(X_test)[:, 1]
        loss_cv.append(np.round(log_loss(y_test, y_pred_prob), 3))
        y_pred_full = np.concatenate([y_pred_full, np.expand_dims(y_pred_prob, -1)])
    y_pred_full_prev = y_pred_full
    loss_all = log_loss(y_test_full, y_pred_full)
    tmp_out = np.array([np.round(np.mean(loss_cv), 5), np.round(np.std(loss_cv), 5), np.round(loss_all, 5)] + loss_cv)
    loss_cv_lst.append(tmp_out)
    print((f, np.round(np.mean(loss_cv), 5), np.round(loss_all, 5)))

loss_df = pd.DataFrame(loss_cv_lst, columns = ['loss_cv_mean', 'loss_cv_sd', 'loss_all'] + ['loss_' + str(i) for i in range(10)])

loss_df = pd.concat((pd.DataFrame({'Food':tmp_f}), loss_df), axis = 1)
loss_df.to_csv(outfile, index = False)

print('finished')



