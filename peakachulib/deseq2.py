import numpy as np
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import r, Formula, pandas2ri
pandas2ri.activate()


class RunDESeq2(object):

    def __init__(self, count_df, exp_lib_list, ctr_lib_list, size_factors):
        r("suppressMessages(library(DESeq2))")
        self._count_df = count_df
        self._exp_lib_list = exp_lib_list
        self._ctr_lib_list = ctr_lib_list
        self._size_factors = size_factors
        
    def run_deseq2(self):
        self._count_df = np.round(self._count_df, decimals=0)
        self._count_df = self._count_df.astype(int)
        libs = self._exp_lib_list + self._ctr_lib_list
        conds = ["exp"] * len(self._exp_lib_list) + ["ctr"] * len(
            self._ctr_lib_list)
        colData = robjects.DataFrame({"conditions": robjects.StrVector(conds)})
        colData.rownames = libs
        #design = r("design <- ~ conditions")
        design = Formula('~ conditions')
        dds = r.DESeqDataSetFromMatrix(countData=self._count_df,
                                       colData=colData, design=design)
        #dds = r.DESeq(dds)
        if self._size_factors is None:
            dds = r.estimateSizeFactors(dds)
        else:
            robjects.globalenv["dds"] = dds
            robjects.globalenv["sf"] = robjects.FloatVector(self._size_factors)
            r("sizeFactors(dds) <- sf")
            dds = robjects.globalenv["dds"]
            r.rm("dds")
        dds = r.estimateDispersions(dds)
        dds = r.nbinomWaldTest(dds)
        size_factors = pd.Series(r.sizeFactors(dds),
                                 index=self._count_df.columns)
        results = r.results(dds, contrast=robjects.StrVector(
            ("conditions", "exp", "ctr")))
        results_df = pandas2ri.ri2py_dataframe(r['as.data.frame'](results))
        results_df.index = self._count_df.index
        return(results_df, size_factors)
