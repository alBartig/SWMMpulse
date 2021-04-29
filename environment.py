import copy

class Units:
    GRAM = 'g'
    MILLIGRAM = 'mg'
    MICROGRAM = 'ug'
    COUNT = '#'

class Constituent:
    def __init__(self, name, specific_load, unit, degradation_coefficient=0, probability=1):
        self.name = name
        self.specific_load = specific_load
        self.unit = unit
        self.degradation_coefficient = degradation_coefficient
        self.probability = probability

    def __repr__(self):
        print(f'{self.name}')

class DefaultConstituents:
    FECAL = Constituent('Fecal Matter', 200, Units.GRAM)
    COV = Constituent('Cov RNA', 10000, Units.COUNT)
    PEP = Constituent('Pepper virus', 10000, Units.COUNT)

class Pattern:
    def __init__(self, name, values):
        self.name = name
        self.values = values

    def __repr__(self):
        print(self.name)

class DefaultPatterns:
    _BRISTOLmen = [1.4,0.3,0.1,0.0,0.3,1.7,9.1,21,13,9,6.9,4.9,1.9,3.6,2.5,2,2.9,2.3,4.1,4.0,2.7,2.1,2.2,2.0]
    _BRISTOLwomen = [2.7,0.1,0.1,0.08,0.02,0.2,3.3,16.9,20.1,12.4,8.8,5.0,2.5,2.6,3.1,2.5,3.0,2.2,4.3,3.3,2.1,1.5,2.0,1.2]
    BRISTOL = Pattern('Bristol',[round((w*0.5 + v*0.5),3) for w,v in zip(_BRISTOLmen,_BRISTOLwomen)])

class Group:
    def __init__(self, name, constituents, weight=None, pattern=DefaultPatterns.BRISTOL):
        self.name = name
        self.constituents = constituents
        self.weight = weight
        self.pattern = pattern

    def __repr__(self):
        print(f'{self.name}, Constituents: {" ".join(self.constituents)}, Weight: {self.weight}, Pattern: {self.pattern}')

    def set(self, **kwargs):
        [setattr(self,k,v) for k,v in kwargs.items()]

class DefaultGroups:
    HEALTHY = Group('Healthy', [DefaultConstituents.FECAL,DefaultConstituents.PEP])
    INFECTED = Group('Infected', [DefaultConstituents.FECAL,DefaultConstituents.PEP,DefaultConstituents.COV])

class Environment:
    def __init__(self, groups, dispersion=1.6):
        self.DISPERSION = dispersion
        self.GROUPS = groups