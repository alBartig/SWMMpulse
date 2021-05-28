from datetime import timedelta

class Distributions:
    BRISTOL = 'Bristol'

class Discretization:
    TIMESTEPLENGTH = timedelta(seconds=10) #timesteplength in s
    SERIESLENGTH = 86400 #length of the timeseries to be calculated in s

class Loading:
    CONSTITUENTS = 'constituents'
    DEGRADATION = 'specific_degradation'
    LOAD = 'specific_loading'
    CLASSIFICATIONS = 'specific_classifications'
    DISPERSION = 'dispersion'
    FECAL = 'Fecal_Matter'
    COV = 'Cov_RNA'
    PEP = 'Pepper_virus'