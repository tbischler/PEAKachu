import numpy as np
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import r, Formula, pandas2ri
pandas2ri.activate()


class RunDESeq2(object):

    def __init__(self, count_df, exp_lib_list, ctr_lib_list, size_factors,
                 pairwise_replicates):
        r("suppressMessages(library(DESeq2))")
        self._count_df = count_df
        self._exp_lib_list = exp_lib_list
        self._ctr_lib_list = ctr_lib_list
        self._size_factors = size_factors
        self._pairwise_replicates = pairwise_replicates

    def run_deseq2(self):
        self._count_df = np.round(self._count_df, decimals=0)
        self._count_df = self._count_df.astype(int)
        conds = ["exp"] * len(self._exp_lib_list) + ["ctr"] * len(
            self._ctr_lib_list)
        if self._pairwise_replicates:
            samples = list(range(1, len(self._exp_lib_list) + 1)) + list(
                range(1, len(self._ctr_lib_list) + 1))
            colData = robjects.DataFrame({
                    "conditions": robjects.StrVector(conds),
                    "samples": robjects.StrVector(samples)})
            design = Formula('~ samples + conditions')
        else:
            colData = robjects.DataFrame(
                    {"conditions": robjects.StrVector(conds)})
            design = Formula('~ conditions')
        r_count_df = robjects.DataFrame(self._count_df)
        r_count_df.colnames = robjects.rinterface.NULL
        dds = r.DESeqDataSetFromMatrix(countData=r_count_df,
                                       colData=colData, design=design)
        if self._size_factors is None:
            dds = r.estimateSizeFactors(dds)
        else:
            assign_sf = r["sizeFactors<-"]
            dds = assign_sf(object=dds, value=robjects.FloatVector(
                self._size_factors))
        dds = r.estimateDispersions(dds, quiet=True)
        dds = r.nbinomWaldTest(dds, quiet=True)
        size_factors = pd.Series(r.sizeFactors(dds),
                                 index=self._count_df.columns)
        results = r.results(dds, contrast=robjects.StrVector(
            ("conditions", "exp", "ctr")), altHypothesis="greater")
        results_df = pandas2ri.ri2py_dataframe(r['as.data.frame'](results))
        results_df.index = self._count_df.index
        return(results_df, size_factors)
