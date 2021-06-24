class StepsizeError(BaseException):
    def __init__(self):
        print('Timeseries length and timestep result in an uneven timestep')

class PacketdictError(BaseException):
    def __init__(self):
        print('Could not append input to timeseries')

class RoutingError(BaseException):
    def __init__(self, packet=None, stops=None):
        print('Error routing a packet - Implausible outcome')

class PlausibilityError(BaseException):
    pass

class MassbalanceError(BaseException):
    def __init__(self):
        print('Expected total load does not match calculated total load')