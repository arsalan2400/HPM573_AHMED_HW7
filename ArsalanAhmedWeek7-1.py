####HOMEWORK 7.1#######

from enum import Enum
import numpy as np

#####note this is just sourced from the survival model in the HPM573 folder#####
class HealthStat(Enum):
    """ health status of patients  """
    ALIVE = 1
    DEAD = 0
####the variables have been renamed a bit because I wanted to play around with it a little
class Patient(object):
    def __init__(self, id, deathrisk):
        # every patient has an id, assorted randomly.
        self._id = id
        self._rnd = np.random
        # There's a risk of death, but all are alive @ Time = 0.
        self._healthcondition = Hospital.ALIVE
        self._timesurvive = 0
        self._risk_of_death = deathrisk

    def simulate(self, timestep):
        t = 0
        while self._healthcondition == Hospital.ALIVE and t < timestep:
            # if you live an extra year, you add 1 yr to survival time.
            # otherwise, you're dead
            if self._rnd.sample() < self._risk_of_death:
                self._healthcondition = Hospital.DEAD
                self._timesurvive = t + 1
            # increments of +1
            t += 1

    def findsurvivaltime(self):
        # shift them out if they die
        if self._healthcondition == Hospital.DEAD:
            return self._timesurvive
        else:
            return None


class Cohort:
    def __init__(self, id, pop_size, deathrisk):
        self._supertimesurvivelist = []
        self._superpatientslist = []
        for i in range(pop_size):
            thispatient = Patient(pop_size * id + i, deathrisk)
            self._superpatientslist.append(thispatient)

    # throw the patients in a simulation
    def simulate(self, timestep):
        for thispatient in self._superpatientslist:
            thispatient.simulate(timestep)
            # this is how we record surv. time with variable 'output'
            output = thispatient.findsurvivaltime()
            if not (output is None):
                self._supertimesurvivelist.append(output)

    #####the above is pretty much cut-n-paste from the survival model in class#####
    def above5years(self):
        startcount = 0
        for x in self._supertimesurvivelist:
            MAXYEAR = 5
            if x > MAXYEAR:
                startcount += 1
        return startcount

    def get_average_survival_time(self):
        return sum(self._supertimesurvivelist) / len(self._supertimesurvivelist)
    ##^^^easier just to use this since I've used it before. Could've used get.mean.
    ##just thought it'd be nice to have this.


SUPERTIMESTEP = 100  # this is somewhat arbitrary, how often the exp repeats
# If I put timesteps to 5, no patients would survive >5 yrs.
POPULATION = 573  # this was requested for

thisCohort = Cohort(id=1, pop_size=POPULATION, deathrisk=0.1)
#note that we've put a 10% likelihood of dying every year! Quite high!!
thisCohort.simulate(SUPERTIMESTEP)
print("% Patients who survived >5 years with 100 timesteps: ", thisCohort.above5years() / POPULATION*100)
print('The avg survival time in yrs is (for context) with 100 timesteps... ', thisCohort.get_average_survival_time())
print("# of Patients who survived >5 years with 100 timesteps: ", thisCohort.above5years())
