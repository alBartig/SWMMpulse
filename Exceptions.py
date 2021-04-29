class StepsizeError:
    def __init__(self):
        print('Timeseries length and timestep result in an uneven timestep')

class PacketdictError:
    def __init__(self):
        print('Could not append input to timeseries')
