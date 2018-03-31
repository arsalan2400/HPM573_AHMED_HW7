#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 19:44:39 2018

@author: Aslan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 17:21:22 2018
@author: Aslan
"""

from enum import Enum
import scipy.stats as stat
import numpy as np
import InOutFunctions as InOutSupport
import StatisticalClasses as StatSupport
import FormatFunctions as FormatSupport
import SurvivalModelClasses as SurvivalCls
from scipy.stats import binom


####THIS IS FROM THE CALIBRATIONCLASS.PY GITHUB#########
class CalibrationColIndex(Enum):
    """ indices of columns in the calibration results cvs file  """
    ID = 0          # cohort ID
    W = 1  # likelihood weight
    MORT_PROB = 2   # mortality probability

#####THE SAME AS CALIBRATION PAGE#####
class Calibration:
    def __init__(self):
        """ initializes the calibration object"""
        np.random.seed(1)   # specifying the seed of the numpy random number generator
        self._cohortIDs = range(POST_N)   # IDs of cohorts to simulate
        self._mortalitySamples = []      # values of mortality probability at which the posterior should be sampled
        self._mortalityResamples = []    # resampled values for constructing posterior estimate and interval
        self._weights = []               # likelihood weights of sampled mortality probabilities
        self._normalizedWeights = []     # normalized likelihood weights (sums to 1)
        self._csvRows = \
            [['Cohort ID', 'Likelihood Weights' ,'Mortality Prob']]  # list containing the calibration results

    def sample_posterior(self):
        """ sample the posterior distribution of the mortality probability """

        # find values of mortality probability at which the posterior should be evaluated
        self._mortalitySamples = np.random.uniform(
            low=POST_L,
            high=POST_U,
            size=POST_N)

        # create a multi cohort
        multiCohort = SurvivalCls.MultiCohort(
            ids=self._cohortIDs,
            mortality_probs=self._mortalitySamples,
            pop_sizes=[SIM_POP_SIZE]*POST_N)
#####THE SAME AS CALIBRATION PAGE#####
        
        # simulate the multi cohort
        multiCohort.simulate(SUPERTIMESTEP)
    
        # calculate the likelihood of each simulated cohort
        for cohort_id in self._cohortIDs:
            mean = multiCohort.get_cohort_mean_survival(cohort_id)
            # construct a gaussian likelihood
            # with mean calculated from the simulated data and standard deviation from the clinical study.
            # evaluate this pdf (probability density function) at the mean reported in the clinical study.
            weight = binom._pmf(x = 800, n = 1146, p=0.5)
            # store the weight
            self._weights.append(weight)

        # normalize the likelihood weights
        sum_weights = np.sum(self._weights)
        self._normalizedWeights = np.divide(self._weights, sum_weights)


#####THE SAME AS CALIBRATION PAGE#####
        # re-sample mortality probability (with replacement) according to likelihood weights
        self._mortalityResamples = np.random.choice(
            a=self._mortalitySamples,
            size=NUM_SIM_COHORTS,
            replace=True,
            p=self._normalizedWeights)

        # produce the list to report the results
        for i in range(0, len(self._mortalitySamples)):
            self._csvRows.append(
                [self._cohortIDs[i], self._normalizedWeights[i], self._mortalitySamples[i]])

        # write the calibration result into a csv file
        InOutSupport.write_csv('CalibrationResults.csv', self._csvRows)
        



#####THE SAME AS CALIBRATION PAGE#####
    def get_mortality_resamples(self):
        """
        :return: mortality resamples
        """
        return self._mortalityResamples

    def get_mortality_estimate_credible_interval(self, alpha, deci):
        """
        :param alpha: the significance level
        :param deci: decimal places
        :returns text in the form of 'mean (lower, upper)' of the posterior distribution"""

        # calculate the credible interval
        sum_stat = StatSupport.SummaryStat('Posterior samples', self._mortalityResamples)

        estimate = sum_stat.get_mean()  # estimated mortality probability
        credible_interval = sum_stat.get_PI(alpha)  # credible interval

        return FormatSupport.format_estimate_interval(estimate, credible_interval, deci)

    def get_effective_sample_size(self):
        """
        :returns: the effective sample size
        """
        return 1 / np.sum(self._normalizedWeights ** 2)

#####THE SAME AS CALIBRATION PAGE#####
class CalibratedModel:
    """ to run the calibrated survival model """

    def __init__(self, cvs_file_name, drug_effectiveness_ratio=1):
        """ extracts seeds, mortality probabilities and the associated likelihood from
        the csv file where the calibration results are stored
        :param cvs_file_name: name of the csv file where the calibrated results are stored
        :param calibrated_model_with_drug: calibrated model simulated when drug is available
        """

        # read the columns of the csv files containing the calibration results
        cols = InOutSupport.read_csv_cols(
            file_name=cvs_file_name,
            n_cols=3,
            if_ignore_first_row=True,
            if_convert_float=True)

        # store likelihood weights, cohort IDs and sampled mortality probabilities
        self._cohortIDs = cols[CalibrationColIndex.ID.value].astype(int)
        self._weights = cols[CalibrationColIndex.W.value]
        self._mortalityProbs = cols[CalibrationColIndex.MORT_PROB.value] * drug_effectiveness_ratio
        self._multiCohorts = None  # multi-cohort

    def simulate(self, num_of_simulated_cohorts, cohort_size, timestep, cohort_ids=None):
        """ simulate the specified number of cohorts based on their associated likelihood weight
        :param num_of_simulated_cohorts: number of cohorts to simulate
        :param cohort_size: the population size of cohorts
        :param timestep: simulation length
        :param cohort_ids: ids of cohort to simulate
        """
        # resample cohort IDs and mortality probabilities based on their likelihood weights
        # sample (with replacement) from indices [0, 1, 2, ..., number of weights] based on the likelihood weights
        sampled_row_indices = np.random.choice(
            a=range(0, len(self._weights)),
            size=num_of_simulated_cohorts,
            replace=True,
            p=self._weights)
        # use the sampled indices to populate the list of cohort IDs and mortality probabilities
        resampled_ids = []
        resampled_probs = []
        for i in sampled_row_indices:
            resampled_ids.append(self._cohortIDs[i])
            resampled_probs.append(self._mortalityProbs[i])

        # simulate the desired number of cohorts
        if cohort_ids is None:
            # if cohort ids are not provided, use the ids stored in the calibration results
            self._multiCohorts = SurvivalCls.MultiCohort(
                ids=resampled_ids,
                pop_sizes=[cohort_size] * num_of_simulated_cohorts,
                mortality_probs=resampled_probs)
        else:
            # if cohort ids are provided, use them instead of the ids stored in the calibration results
            self._multiCohorts = SurvivalCls.MultiCohort(
                ids=cohort_ids,
                pop_sizes=[cohort_size] * num_of_simulated_cohorts,
                mortality_probs=resampled_probs)

        # simulate all cohorts
        self._multiCohorts.simulate(timestep)


#####THE SAME AS CALIBRATION PAGE#####
    def get_all_mean_survival(self):
        """ :returns a list of mean survival time for all simulated cohorts"""
        return self._multiCohorts.get_all_mean_survival()

    def get_mean_survival_time_proj_interval(self, alpha, deci):
        """
        :param alpha: the significance level
        :param deci: decimal places
        :returns text in the form of 'mean (lower, upper)' of projection interval
        """

        mean = self._multiCohorts.get_overall_mean_survival()
        proj_interval = self._multiCohorts.get_PI_mean_survival(alpha)

        return FormatSupport.format_estimate_interval(mean, proj_interval, deci)
####END OF CALIBRATION PAGE COPY######
        
    
####NOW WE JUST NEED TO SET VALUES FOR ALL THESE INPUTS!!!!!#####
#'Things we need to know....')
#'Step 1. define the SIM_POP_SIZE in multi cohort = 1000.')
SIM_POP_SIZE = 1000     # population size of the simulated cohort
#'Step 2. Set a value for TIMESTEP multiCohort.simulate(CalibSets.SUPERTIMESTEP =100.')
#I'm just using SUPERTIMESTEP instead of TIME_STEPS here because it's from my Question 7.1 
SUPERTIMESTEP = 100
#Step 3. Set the alpha =0.05 ')
ALPHA = 0.05 #sig level
#'Step 4. Remember, there is a re-sample mortality probability (with replacement) according to likelihood weights aka size=CalibSets.NUM_SIM_COHORTS = 300 because why not')
NUM_SIM_COHORTS = 300   # number of simulated cohorts used to calculate prediction intervals
#'Step 5. Set upper and lower limit. 0.05 and 0.25 is reasonable')
POST_L, POST_U, POST_N = 0.05, 0.25, 1000
#'Step 6. We need to put in our values for our clinical study. We need Obs_mean survival time and Obs_STDEV')


print('The estimate of new mortality probability ({:.{prec}%} credible interval):'.format(1-ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(ALPHA, 4))
#an answer:
#Estimate of new mortality probability (95% credible interval): 0.1500 (0.0533, 0.2477)
#note that the range is smaller. 
print('The credible interval of the estimated annual mortality is smaller. Possibly due to more confidence w larger #s.')

this_calibrated_model= CalibratedModel('CalibrationResults.csv')
this_calibrated_model.simulate(NUM_SIM_COHORTS,SIM_POP_SIZE, SUPERTIMESTEP)
print('The new mean survival time and {:.{prec}%} projection interval:'.format(1 - ALPHA, prec=0),
      this_calibrated_model.get_mean_survival_time_proj_interval(ALPHA, deci=4))
print('The proj interval the mean survival time is smaller. Possibly due to more confidence w larger #s.')
#an answer: 
#E.g. The new mean survival time and 95% projection interval: 8.1106 (4.0057, 18.1507).

#a complete answer:
#The estimate of new mortality probability (95% credible interval): 0.1500 (0.0533, 0.2477)
#The credible interval of the estimated annual mortality is smaller. Possibly due to more confidence w larger #s.
#The new mean survival time and 95% projection interval: 7.9504 (4.1620, 17.9836)
#The proj interval the mean survival time is slightly smaller. Possibly due to more confidence w larger #s.
