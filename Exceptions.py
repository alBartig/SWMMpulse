class StepsizeError(BaseException):
    def __init__(self):
        print('Timeseries length and timestep result in an uneven timestep')

class PacketdictError(BaseException):
    def __init__(self):
        print('Could not append input to timeseries')

class RoutingError(BaseException):
    print('Error routing a packet - Implausible outcome')
    def __init__(self, packet=None, stops=None):
        print(packet)
        print(stops)

class PlausibilityError(BaseException):
    pass