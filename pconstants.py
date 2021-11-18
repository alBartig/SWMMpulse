from datetime import timedelta

class Distributions:
    BRISTOL = 'Bristol'

class Discretization:
    TIMESTEPLENGTH = timedelta(seconds=10) #timesteplength in s
    SERIESLENGTH = 86400 #length of the timeseries to be calculated in s
    REFINE = 10

class Loading:
    CONSTITUENTS = 'constituents'
    DEGRADATION = 'specific_degradation'
    LOAD = 'specific_loading'
    CLASSIFICATIONS = 'specific_classifications'
    DISPERSION = 'dispersion'
    FECAL = 'Fecal-Matter'
    COV = 'Cov-RNA'
    PEP = 'Pepper-virus'
    FRACTION = "Load-fraction"
    SKEWEDNESS = "skewedness"
    MAX_CONTINUITY_ERROR = 0.02
    FRACTIONS = {COV : [{FRACTION: 0.5, DISPERSION: 0.1, SKEWEDNESS: 1},
                        {FRACTION:0.5, DISPERSION:0.1, SKEWEDNESS:8}],
                 FECAL : [{FRACTION:0.5, DISPERSION:0.1, SKEWEDNESS:8},
                          {FRACTION:0.5, DISPERSION:0.1, SKEWEDNESS:8}]}