import numpy as np
from pandas import read_csv
from scipy import sparse
from scipy.signal import savgol_filter, argrelextrema
from scipy.sparse.linalg import spsolve
from numpy.fft import rfft, irfft
from pprint import pprint
import csv
import math
import sys


# get extremum on range
def get_extremum(arr: np.array, type: str = 'min')->np.array:

    answer = None
    coords = None

    if type == 'min':
        answer = np.amin(arr, axis=1)[1]
    elif type == 'max':
        answer = np.amax(arr, axis=1)[1]

    for x in arr[1]:
        if x == answer:
            index = list(arr[1]).index(x)
            coords = [arr[0][index], arr[1][index]]

    return coords


# write the csv files
def write_csv_files(bins: dict):
    for bin in bins:
        with open(f'results/1/{bin}.csv', 'w') as file:
            writer = csv.writer(file, delimiter=',', quoting=csv.QUOTE_MINIMAL)
            samples_dict = bins[bin]['samples']
            for sample in samples_dict:
                writer.writerow(str(sample))
                writer.writerow([samples_dict[sample]['answer_range'][0][0], samples_dict[sample]['answer_range'][1][0], samples_dict[sample]['positive_answer'], samples_dict[sample]['minimum'], samples_dict[sample]['maximum']])
                for x in range(0, len(samples_dict[sample]['answer_range'][0])):
                    writer.writerow([samples_dict[sample]['answer_range'][0][x], samples_dict[sample]['answer_range'][1][x]])
                writer.writerow('-'*15)
            file.close()


# find index of nearest value to x in array
def find_nearest(array: np.array, value: float):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# get Deckard distance
def distance(p0: list, p1: list)->float:
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)


# check the extrema
def check_the_extrema(sample: dict):
    # check the peaks

    # define which extremum is closer to infra-red range and which to ultra-violet
    right_extremum = None
    if 435.0 > sample['minimum'][0] > 405.0 and sample['minimum'][0] > sample['maximum'][0]: # 405 - is additional check for the
        right_extremum = sample['minimum']
    elif 405.0 < sample['maximum'][0] < 435.0 and sample['minimum'][0] < sample['maximum'][0]:
        right_extremum = sample['maximum']

    # check if the right extremum exists
    if right_extremum is not None:

        # make the fft (real), collect only first 3 frequencies and make ifft (real)
        rft = rfft(sample['spectrum'])
        rft[3:] = 0
        irft = irfft(rft)

        # get the polynomial fit of the ifft curve with appropriate range
        polyfit_range = 3
        r_squar = 0
        coeffs = None
        p = None
        while r_squar < 0.99:
            coeffs = np.polyfit(wavelength[:-1], irft, polyfit_range)
            # r-squared
            p = np.poly1d(coeffs)
            # fit values, and mean
            yhat = p(wavelength)  # or [p(z) for z in x]
            ybar = np.sum(irft) / len(irft)  # or sum(y)/len(y)
            ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
            sstot = np.sum((irft - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
            r_squar = ssreg / sstot

            polyfit_range += 1

        if coeffs is not None:
            coeffs = list(coeffs)
            coeffs.reverse()

            deriv_poly = [coeffs[i] * i for i in range(1, len(coeffs))]
            deriv_poly.reverse()
            roots = np.roots(deriv_poly)
            print(f'roots: {roots}')

            test_points = [p(root) for root in roots]
            test_points.insert(0, p(wavelength[0]))
            test_points.append(p(wavelength[-1]))

            # test roots and remove ones that do not illustrate real extrema
            confirmed_roots = []
            for root_index in range(0, len(roots)):
                left_y = test_points[root_index]    #
                root_y = test_points[root_index+1]  # because test points contain additional elements at the beginning and the end
                right_y = test_points[root_index+2] #

                if left_y > root_y < right_y or left_y < root_y > right_y:
                    confirmed_roots.append(roots[root_index])

            print(f'confirmed roots: {confirmed_roots}')
            if len(confirmed_roots) != len(roots):
                print('false extrema removed', file=sys.stderr)

        pprint(f"answer maximum: {sample['maximum']}")
        pprint(f"answer minimum: {sample['minimum']}")
        print('-'*10)



# baseline correction
# y - 1D array of y_coords, lam - lambda for smoothness (10^2 <= λ <= 10^9), p for asymmetry (0.001 <= p <= 0.1)
def baseline_correction(y, lam: int = 10000, p: float = 0.01, niter: int =15):
  L = len(y)
  D = sparse.csc_matrix(np.diff(np.eye(L), 2))
  w = np.ones(L)
  for i in range(niter):
    W = sparse.spdiags(w, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = spsolve(Z, w*y)
    w = p * (y > z) + (1-p) * (y < z)

  return z


if __name__ == '__main__':
    data = read_csv('Data.csv', header=0)

    # get the wavelength array
    wavelength = [int(x) for x in data['number'][1:]]

    # collect sets in bins
    bins = {}
    index = 0
    for column in data.columns[1:]:
        index += 1
        weight = data[column][0]
        if weight == '?':
            continue
        weight = float(weight)
        if weight not in bins:
            bins[weight] = {'samples': {}}

        if sum([x*x for x in data[column][1:]]) == 0.0:
            print(f'column {column} skipped')
            continue
        else:
            # create an key for each sample and description for it
            bins[weight]['samples'][index] = {'spectrum': [], 'smoothed_spectrum': [], 'minimum': float,
                                              'maximum': float, 'answer_range': np.array, 'positive_answer': False}
            sample = bins[weight]['samples'][index]

            # fill the spectrum after savgol smothering
            sample['spectrum'] = [x for x in data[column][1:]]
            baseline = baseline_correction(sample['spectrum'])
            corrected_spectrum = []
            for i in range(0, len(sample['spectrum'])):
                corrected_spectrum.append(sample['spectrum'][i] - baseline[i])
            sample['smoothed_spectrum'] = [x for x in savgol_filter(corrected_spectrum, 15, 3)]

            # cut the spectra
            first = wavelength.index(380)
            last = wavelength.index(450)

            # rotate the answer range and get extrema
            sample['answer_range'] = np.array([[float(w) for w in wavelength[first:last]], [b for b in sample['smoothed_spectrum'][first:last]]])

            sample['minimum'] = get_extremum(sample['answer_range'])
            sample['maximum'] = get_extremum(sample['answer_range'], type='max')

            # check the extrema of spectrum
            print(index)
            check_the_extrema(sample)

    write_csv_files(bins)