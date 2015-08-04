import pandas as pd
import numpy as np
from scipy.stats import chisqprob

class GTest(object):
    
    def __init__(self, ctr_rep_counts, tagged_rep_counts,
            pairwise_replicates=False):
        self._rep_df = pd.DataFrame({"ctr_counts": ctr_rep_counts,
                                     "tagged_counts": tagged_rep_counts})
        self._pairwise_replicates = pairwise_replicates
        self._pooled_g_res = {}
        self._total_g_res = {}
        self._heterogenous_g_res = {}
        self._single_g_res = {}
        
    def _olnf(self, obs, exp):
        return obs * np.log(obs/exp) if obs > 0.1 else 0
    
    def _g_to_p_value(self, g_value, dof):
         return chisqprob(g_value, dof)
    
    def _gtest(self, values):
        ctr_obs = values[0]
        tagged_obs = values[1]
        ctr_exp = 0.5 * (ctr_obs + tagged_obs)
        tagged_exp = 0.5 * (ctr_obs + tagged_obs)
        return 2 * (self._olnf(ctr_obs, ctr_exp) + 
                    self._olnf(tagged_obs, tagged_exp))
    
    def _pooled_gtest(self):
        self._pooled_g_res["dof"] = 1
        pooled_ctr = self._rep_df.loc[:, "ctr_counts"].sum()
        pooled_tagged = self._rep_df.loc[:, "tagged_counts"].sum()
        self._pooled_g_res["g_value"] = self._gtest((pooled_ctr, pooled_tagged))
        self._pooled_g_res["p_value"] = self._g_to_p_value(
                self._pooled_g_res["g_value"],
                self._pooled_g_res["dof"])
    
    def _total_gtest(self):
        '''Use maximum over all g-value sums for all possible replicate
           combinations. Instead using the mean would also be possible!
        '''
        g_value_dict = {}
        ctr_counts = self._rep_df.loc[:, "ctr_counts"]
        for comb_iter in range(len(self._rep_df.index)):
            tagged_counts = self._rep_df.loc[:, "tagged_counts"]
            comb_df = pd.concat([ctr_counts.reset_index(drop=True),
                    tagged_counts.iloc[comb_iter:].append(
                    tagged_counts.iloc[0:comb_iter]).reset_index(drop=True)],
                    axis=1)
            rep_g_values = comb_df.apply(self._gtest, axis=1)
            rep_p_values = rep_g_values.apply(self._g_to_p_value, args=(1,))
            tot_g_value = rep_g_values.sum()
            g_value_dict[tot_g_value] = max_rep_p_value
        self._total_g_res["g_value"] = max(g_value_dict.keys())
        self._replicate_p_values = g_value_dict[
                self._total_g_res["g_value"]]
        # Degrees of freedom = replicate number
        self._total_g_res["dof"] = len(self._rep_df.index)
        # Calculate p-value
        self._total_g_res["p_value"] = self._g_to_p_value(
                self._total_g_res["g_value"],
                self._total_g_res["dof"])
        
    def _total_gtest_pairwise(self):
        '''Use only replicate pairs according to library input order
        '''
        rep_g_values = self._rep_df.apply(
                self._gtest, axis=1)
        self._replicate_p_values = rep_g_values.apply(
            self._g_to_p_value, args=(1,))
        self._total_g_res["g_value"] = rep_g_values.sum()
        # Degrees of freedom = replicate number
        self._total_g_res["dof"] = len(self._rep_df.index)
        # Calculate p-value
        self._total_g_res["p_value"] = self._g_to_p_value(
                self._total_g_res["g_value"],
                self._total_g_res["dof"])
        
    def _heterogenous_gtest(self):
        self._heterogenous_g_res["g_value"] = (self._total_g_res["g_value"] -
                                              self._pooled_g_res["g_value"])
        self._heterogenous_g_res["dof"] = (self._total_g_res["dof"] -
                                          self._pooled_g_res["dof"])
        self._heterogenous_g_res["p_value"] = self._g_to_p_value(
                self._heterogenous_g_res["g_value"],
                self._heterogenous_g_res["dof"])
    
    def run_with_repl(self):
        self._pooled_gtest()
        (self._total_gtest_pairwise() if self._pairwise_replicates
                                    else self._total_gtest())
        self._heterogenous_gtest()
        return {"replicate_G_p_values": "/".join([str(value) for value in list(
                self._replicate_p_values)]),
                "pooled_G_p_value": self._pooled_g_res["p_value"],
                "total_G_p_value": self._total_g_res["p_value"],
                "heterogenous_G_p_value": self._heterogenous_g_res["p_value"]}
    
    def run_without_repl(self):
        self._single_g_res["dof"] = 1
        ctr = self._rep_df.loc[:, "ctr_counts"].sum()
        tagged = self._rep_df.loc[:, "tagged_counts"].sum()
        self._single_g_res["g_value"] = self._gtest((ctr, tagged))
        self._single_g_res["p_value"] = self._g_to_p_value(
                self._single_g_res["g_value"],
                self._single_g_res["dof"])
        return {"single_G_p_value": self._single_g_res["p_value"]}