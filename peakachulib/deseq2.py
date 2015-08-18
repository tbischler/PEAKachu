import numpy as np
import pandas as pd
from rpy2 import robjects
robjects.pandas2ri.activate()


class RunDESeq2(object):

    def __init__(self, count_df, exp_lib_list, ctr_lib_list):
        robjects.r("suppressMessages(library(DESeq2))")
        self._count_df = count_df
        self._exp_lib_list = exp_lib_list
        self._ctr_lib_list = ctr_lib_list
        
    def run_deseq2(self):
        self._count_df = np.round(self._count_df, decimals=0)
        self._count_df = self._count_df.astype(int)
        libs = self._exp_lib_list + self._ctr_lib_list
        conds = ["exp"] * len(self._exp_lib_list) + ["ctr"] * len(
            self._ctr_lib_list)
        colData = robjects.DataFrame({"conditions": robjects.StrVector(conds)})
        colData.rownames = libs
        design = robjects.r("design <- ~ conditions")
        dds = robjects.r.DESeqDataSetFromMatrix(countData=self._count_df,
                                                colData=colData, design=design)
        dds = robjects.r.DESeq(dds, minReplicatesForReplace=np.inf)
        size_factors = pd.Series(robjects.r.sizeFactors(dds),
                                 index=self._count_df.columns)
        results = robjects.r.results(dds, cooksCutoff=False,
                                     independentFiltering=False,
                                     contrast=robjects.StrVector(
            ("conditions", "exp", "ctr")))
        results_df = robjects.pandas2ri.ri2py_dataframe(robjects.r[
            'as.data.frame'](results))
        results_df.index = self._count_df.index
        return(results_df, size_factors)
