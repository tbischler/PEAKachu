import numpy as np
import pandas as pd
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()


class TMM(object):

    def __init__(self, count_df):
        r("suppressMessages(library(edgeR))")
        self.count_df = count_df

    def calc_size_factors(self):
        # Convert pandas dataframe to R dataframe
        r_dge = r.DGEList(self.count_df)
        # Calculate normalization factors
        r_dge = r.calcNormFactors(r_dge, method="TMM")
        size_factors = (np.array(r_dge.rx2('samples')["lib.size"]) *
                        np.array(r_dge.rx2("samples")["norm.factors"]))
        # convert to pandas series
        size_factors = pd.Series(size_factors, index=self.count_df.columns)
        # adjust size factors so that the maximum is 1.0
        size_factors = size_factors/size_factors.max()
        return size_factors
