print('Question 7.2')
print('We know, based on the >5 year survivors, this is Binomial Dist with parameters: n=573, p=0.41)
      #q = one run gave me 0.5863874345549738; we'll go with it. 
print('\n')
print('Question 7.3')
import scipy.stats as stat
prob=stat.binom._pmf(400, 573, 0.5)
print('The chances of 400/573 people surviving @ 5 year end if 50% of the patients in our simulated cohort survived beyond 5 years is...', prob)

