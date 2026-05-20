#!/usr/bin/env python3

import os
import sys
import os.path as osp
import argparse as arp

import time
import numpy as np
import pysindy as ps
import csv
import pandas as pd
from rich.progress import track
from sklearn.linear_model import Ridge
from sklearn.model_selection import RepeatedKFold,GridSearchCV

from multiprocessing import Process
import warnings
warnings.filterwarnings('ignore')

"""
Usage:
python WISpR_new.py -r single_cell_data -s Spatial_Data -o /Directory/To/Save
"""


## stlsq optimizer: Its default model is Ridge Regression
def stlsq_w(A,y,best_res):
    """
    Inputs
    A: DEG vs Cell Type matrix
    y: Capture spots
    best_res: Output of GridSearch algorithm
    Output:
    x: Coefficients to be used in initial prediction
    """
    stlsq = ps.STLSQ(threshold=best_res.best_params_['threshold'], alpha= best_res.best_params_['alpha'], fit_intercept=True, copy_X=False)
    model = ps.SINDy(feature_library=ps.IdentityLibrary(), optimizer=stlsq, discrete_time=True)
    model.fit(A, x_dot=y)
    x = model.coefficients()
    return x


def weighting(A,y):
    """
    Inputs
    A: DEG vs Cell Type matrix
    y: Capture spots
    Output:
    x: Coefficients to be used in weight calculation
    """
  #  ridge_optimizer = Ridge(alpha=0.1, max_iter=1000, fit_intercept=True, copy_X=False, solver="auto", positive = True) #orginal
    ridge_optimizer = Ridge(alpha=0.1, max_iter=1000, fit_intercept=True, copy_X=False, solver="auto", positive = True, normalize=True) #for subcluster weighting
    solution= ridge_optimizer.fit(X=A, y=y)
    x = solution.coef_
    return x

def deconvolve(sc_cnt, st_cnt, sc_cnt_df_t):
    """
    Inputs
    sc_cnt: DEG vs Cell Type matrix
    st_cnt: DEG vs Spot Matrix
    sc_cnt_df_t: Transpose of dataframe for DEG vs Cell Type matrix
    Output:
    output: Deconvoluted matrix
    """
    start = time.time()
    small_indices = np.zeros((4, sc_cnt.shape[1],st_cnt.shape[1]), float)
    predict_ = np.zeros((3, sc_cnt.shape[1],st_cnt.shape[1]), float)
    predict_dff=pd.DataFrame(np.empty((st_cnt.shape[1], sc_cnt.shape[1])))
    p_weight = np.zeros((sc_cnt.shape[0], st_cnt.shape[1]), float)
    output=np.zeros((st_cnt.shape[1], sc_cnt.shape[1]), float)

    cv = RepeatedKFold(n_splits=5, n_repeats=3, random_state=1)
    model = ps.STLSQ()
    # For more crowded spots use these parameters
    """
    param_grid = {'alpha': np.arange(0.3, 0.5, 0.1),
                 'threshold': np.arange(0.02, 0.05, 0.001),
                 'fit_intercept': ['True'],
                 'normalize_columns': ['True'],
                 'copy_X': ['False']
                 }
                 
    
    """
    #For Visium datasets use these parameters
    param_grid = {'alpha': np.arange(0.0, 0.3, 0.1), #original range (0.0, 0.3, 0.1)
                 'threshold': np.arange(0.001, 0.01, 0.001), #original
               # 'threshold': np.arange(0.01, 0.47, 0.001),
                 'fit_intercept': ['True'],
                 'normalize_columns': ['True'],
                 'copy_X': ['False']
                 }

    gscv = GridSearchCV( \
        model, param_grid, scoring='neg_mean_absolute_error', \
        cv=cv, n_jobs=-1)
        
    def weight(A, y):
        sol = np.empty((1, A.shape[1]), float)
        sol = weighting(A.astype(float), y.astype(float))
        sol_t = sol.transpose()
        w=A.astype(float)@sol_t.astype(float)
        w_weight = 1/w
        w_weight[w != 0] = 1/w[w != 0]; w_weight[w == 0] = 0
        if (np.sum(w_weight) == 0):
            w_weight[y != 0] = 1/y[y != 0]; w_weight[y == 0] = 0
        return(w_weight)

    for i in track(range(st_cnt.shape[1]), description="Deconvoluting data"):
        if (np.sum(st_cnt[:,i]) ==0):
            continue
        else:
            p_weight[:,i] = weight(sc_cnt.astype(float), st_cnt.astype(float)[:,i])
            best_res = gscv.fit(sc_cnt.astype(float), st_cnt.astype(float)[:,i], sample_weight=(p_weight[:,i]))
           # print('MAE: %.5f' % best_res.best_score_)
           # print('Config: %s' % best_res.best_params_)

            ## Nested regression with ridge optimizer with >0 weights
            def stlsq_nested_n(A,y,best_res):
                ridge_optimizer = Ridge(alpha=best_res.best_params_['alpha'], max_iter=1000, fit_intercept=True, copy_X=False, solver="auto", positive = True) #original
              #  ridge_optimizer = Ridge(alpha=best_res.best_params_['alpha'], max_iter=1000, fit_intercept=True, copy_X=False, solver="auto", positive = True, normalize=True)
                solution= ridge_optimizer.fit(X=A, y=y, sample_weight=(p_weight[:,i]))
                x = solution.coef_
                return x

            predict_[0,:,i] = stlsq_w(sc_cnt.astype(float), st_cnt.astype(float)[:,i], best_res) #initial guess
           # threshold = np.arange(0.001, 0.01, 0.001) #original
            threshold = np.arange(0.001, 0.01, 0.001)
            k=0
            for k in range(threshold.shape[0]):
               # print("k: ",k)
                small_indices[0,:,i] = predict_.astype(float)[0,:,i] < best_res.best_params_['threshold']
                small1 = small_indices[0,:,i].astype(bool)
                big_indices = ~small1
                if ((big_indices==0).all(axis=0)):
                    small_indices[1,:,i] = predict_.astype(float)[0,:,i] < threshold[-k]
                    small2 = small_indices[1,:,i].astype(bool)
                    big_indices2 = ~small2
                    predict_[0,small2,i] =0.0
                    predict_dff = pd.DataFrame(predict_[0,:,i]).transpose()
                    predict_dff.columns = sc_cnt_df_t.index
                    output[i,:] = predict_dff
                  #  print("a")
                    if (np.sum(output[i,:])>0):
                        break
                else:
                    predict_[0,small1,i] =0.0
                    predict_[1,big_indices,i] = stlsq_nested_n(sc_cnt.astype(float)[:,big_indices], st_cnt.astype(float)[:,i], best_res)
                    small_indices[1,:,i] = predict_.astype(float)[1,:,i] < best_res.best_params_['threshold']
                    small2 = small_indices[1,:,i].astype(bool)
                    big_indices2 = ~small2
                    if ((big_indices2==0).all(axis=0)):
                        small_indices[2,:,i] = predict_.astype(float)[1,:,i] < threshold[-k]
                        small3 = small_indices[2,:,i].astype(bool)
                        big_indices3 = ~small3
                        predict_[1, small3,i] =0.0
                        predict_dff = pd.DataFrame(predict_[1,:,i]).transpose()
                        predict_dff.columns = sc_cnt_df_t.index
                        output[i,:] = predict_dff
                     #   print("b")
                        if (np.sum(output[i,:])>0):
                            break
                    else:
                        predict_[1, small2,i] =0.0
                        predict_[2,big_indices2,i] = stlsq_nested_n(sc_cnt.astype(float)[:,big_indices2], st_cnt.astype(float)[:,i], best_res)
                        if (predict_.astype(float)[2,:,i] >= 0).all(axis=0):
                            small_indices[3,:,i] = predict_.astype(float)[2,:,i] < best_res.best_params_['threshold']
                            small4 = small_indices[3,:,i].astype(bool)
                            predict_[2, small4,i] =0.0
                            big_indices4 = ~small4
                            predict_dff = pd.DataFrame(predict_[2,:,i]).transpose()
                            predict_dff.columns = sc_cnt_df_t.index
                            #Sum each subcluster members to the corresponding main cluster
                            output[i,:] = predict_dff
                         #   print("c")
                            break

    end = time.time()
    return(pd.DataFrame(output))
    # total time taken
    print(f"Runtime of the program for real mixture is {end - start}")
    
def main():
    prs = arp.ArgumentParser()
    prs.add_argument('-r','--reference',
                     type = str,
                     required = True,
                     help = ' '.join(["path to reference",
                                      "count data",]
                                    )
                    )

    prs.add_argument('-s','--spatial',
                     type = str,
                     required = True,
                     help = ' '.join(["path to spatial",
                                     "data",
                                     ],
                                     )
                    )

    prs.add_argument('-o','--out_dir',
                     default = None,
                     help = 'output directory',
                    )

    args = prs.parse_args()

    if args.out_dir is None:
        out_dir = osp.dirname(args.sc_cnt_df)
    else:
        out_dir = args.out_dir

    if not osp.exists(out_dir):
        os.mkdir(out_dir)

    sc_cnt_pth =  args.reference
    st_cnt_pth = args.spatial

    # read data
    sc_cnt_df = pd.read_csv(sc_cnt_pth,
                         sep = ',',
                         index_col = 0,
                         header = 0)

    st_cnt_df = pd.read_csv(st_cnt_pth,
                         sep = ',',
                         index_col = 0,
                         header = 0)
   # sc_cnt = np.array(pd.read_csv(sc_cnt_pth))
   # st_cnt = np.array(pd.read_csv(st_cnt_pth))

    # match count and label data
    inter = sc_cnt_df.index.intersection(st_cnt_df.index)

    sc_cnt = np.array(sc_cnt_df.loc[inter,:])
    st_cnt = np.array(st_cnt_df.loc[inter,:])
    
    sc_cnt_df_t = sc_cnt_df.transpose()
    
    result = deconvolve(sc_cnt,
                    st_cnt,
                    sc_cnt_df_t)
                    
    result.columns=sc_cnt_df.columns.values
    result.index=st_cnt_df.columns.values
    
    result.to_csv(osp.join(out_dir,
                                   '.'.join(['WISpR_Visium',
                                             osp.basename(sc_cnt_pth)]
                                           )
                                  ),
                          sep = '\t',
                          header = True,
                          index = True)


if __name__ == '__main__':
    main()


