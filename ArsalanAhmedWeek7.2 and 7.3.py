print('Question 7.2')
print('We know, based on the >5 year survivors, this is Binomial Dist with parameters: N= 573 patients, p= about 0.42, q = about 0.58')
print('\n')
print('Question 7.3')
import scipy.stats as stat
prob=stat.binom._pmf(400, 573, 0.58)
#I got about 0.58 for the prob from my single trial.
print('The chances of 400/573 people surviving @ 5 year end if 50% of the patients in our simulated cohort survived beyond 5 years is...', prob)

