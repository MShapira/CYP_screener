from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from pandas import *
from mathtools import baseline_als


class Well:
    def __init__(self, position: str = '', name: str = '', type: str = '', raw_data: list = list()):
        self.position = position
        self.name = name
        if self.name == '':
            self.name = self.position
        self.type = type
        self.spectrum = dict()
        self.raw_data = raw_data

    def __str__(self):
        return 'Name: {0}\n'.format(self.name) + \
               'Position: {0}\n'.format(self.position) + \
               'Type: {0}\n'.format(self.type) + \
               '\n'

    # generate graph
    def generate_graph(self, prot: str, save: bool = False):
        if self.type != 'blank':
            x = np.fromiter(self.spectrum.keys(), dtype=float)
            y = np.fromiter(self.spectrum.values(), dtype=float)
            # smooth the curve
            y_smoothed = savgol_filter(y, 15, 3)  # window size 15, polynomial order 3
            bsl = baseline_als(y_smoothed)
            fig = plt.figure()
            plt.plot(x, y, color='blue')
            plt.plot(x, y_smoothed, color='red')
            plt.plot(x, y_smoothed - bsl, color='green')
            plt.xlabel('Wavelength, nm')
            plt.ylabel('Optical Density, A')
            plt.title('{0}_{1}'.format(prot, self.name))
            if save:
                fig.savefig('{0}_{1}.png'.format(prot, self.name))
            else:
                plt.show()


class Plate:
    def __init__(self, number: int = 0, name: str = '', wells: list = list(), raw_data: list = list()):
        self.name = name
        self.number = number
        self.wells = wells
        self.raw_data = raw_data
        self.protein_well = ''
        self.working_area = list()

    def __str__(self):
        return 'Plate number: {0}\n'.format(self.number) + \
               'Plate name: {0}\n'.format(self.name) + \
               'Plate protein well : {0}\n'.format(self.protein_well) + \
               '\n'

    # get filed area
    def get_area(self):
        # if two lists of column names are equal then the range is common
        if list(self.raw_data[0].columns) == list(self.raw_data[1].columns):
            self.working_area =  list(self.raw_data[0].columns)[1:]

    # calculate only one half (sample part) of the difference spectrum equation
    def calculate_difference_for_well(self, well_position: str)-> list():
        difference = list()

        for i in range(0, len(self.raw_data[1][well_position].get_values())):
            mid_data = self.raw_data[1][well_position].get_values()[i] - self.raw_data[0][well_position].get_values()[i]
            difference.append(mid_data)

        return difference
