# The elbow method to find number of PCs to use
# Reference: 10.1186/s13059-022-02761-4
# Select K by choosing the point that is the farthest from the diagonal line,
# i.e., the line that passes through the first point and the last point (the PCA scree plot)
# Specifically, the distance from (x0,y0) to the line that passes through (x1,y1) and (x2,y2) is given by
# |(x2 −x1)(y1 −y0)−(x1 −x0)(y2 −y1) | /((x2 −x1)^2 +(y2 −y1)^2)^0.5

import pandas as pd
import numpy as np

def get_n_pcs_by_elbow(df_pca_var_explained):
    '''
    Implement the elbow method to find the optimal number of PCs
    Select K by choosing the point that is the farthest from the line that
    passes through the first point (x1,y1), and the last point (x2,y2).
    Specifically, the distance from (x0,y0) to the line that passes through (x1,y1) and (x2,y2) is given by
    | (x2 −x1)(y1 −y0)−(x1 −x0)(y2 −y1) | /((x2 −x1)^2 +(y2 −y1)^2)^(1/2)
    Params:
    - df_pca_var_explained: a dataframe contains columns explained_variance and PC
                            explained_variance can be values from pca.explained_variance_ or pca.explained_variance_ratio
    '''
    df_pca_var_explained = df_pca_var_explained[['PC', 'explained_variance']].copy()
    x1, y1 = df_pca_var_explained.sort_values(by='explained_variance').iloc[-1, :]
    x2, y2 = df_pca_var_explained.sort_values(by='explained_variance').iloc[0, :]
    df_pca_var_explained['distance'] = np.abs((x2-x1)*(y1-df_pca_var_explained['explained_variance'])-(x1-df_pca_var_explained['PC'])*(y2-y1))/((x2-x1)**2 +(y2-y1)**2)**0.5
    k = df_pca_var_explained.sort_values(by='distance')['PC'].iloc[-1]
    return k, df_pca_var_explained