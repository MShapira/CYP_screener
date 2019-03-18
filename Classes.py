from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from pandas import *
from mathtools import baseline_als


class Well:
    def __init__(self, position: str = '', name: str = '', type: str = '', spectrum: dict = {}):
        self.position = position
        self.name = name
        if self.name == '':
            self.name = self.position
        self.type = type
        self.spectrum = spectrum

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