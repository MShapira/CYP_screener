from pandas import DataFrame
import numpy as np
from numpy.fft import rfft, irfft
import sys
from pprint import pprint


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


# make all checks
def get_the_answer(answer_range):

    spectra_minimum = get_extremum(answer_range)
    spectra_maximum = get_extremum(answer_range, type='max')

    # define which extremum is closer to infra-red range and which to ultra-violet
    right_extremum = None
    if 435.0 > spectra_minimum[0] > 405.0 and spectra_minimum[0] > spectra_maximum[
        0]:  # 405 - is additional check for the
        right_extremum = spectra_minimum
    elif 405.0 < spectra_minimum[0] < 435.0 and spectra_minimum[0] < spectra_maximum[0]:
        right_extremum = spectra_maximum

    # check if the right extremum exists
    if right_extremum is not None:

        # make the fft (real), collect only first 3 frequencies and make ifft (real)
        rft = rfft(answer_range[1])
        rft[3:] = 0
        irft = irfft(rft)

        # get the polynomial fit of the ifft curve with appropriate range
        polyfit_range = 3
        r_squar = 0
        coeffs = None
        p = None
        while r_squar < 0.999:
            coeffs = np.polyfit(answer_range[0], irft, polyfit_range)
            # r-squared
            p = np.poly1d(coeffs)
            # fit values, and mean
            yhat = p(answer_range[0])  # or [p(z) for z in x]
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

            test_points = [p(root) for root in roots]
            test_points.insert(0, p(answer_range[0][0]))
            test_points.append(p(answer_range[0][-1]))

            # test roots and remove ones that do not illustrate real extrema
            confirmed_roots = {}
            for root_index in range(0, len(roots)):
                left_y = test_points[root_index]  #
                root_y = test_points[
                    root_index + 1]  # because test points contain additional elements at the beginning and the end
                right_y = test_points[root_index + 2]  #

                if left_y > root_y < right_y or left_y < root_y > right_y:
                    if answer_range[0][0] <= roots[root_index] <= answer_range[0][-1]:
                        confirmed_roots[roots[root_index]] = test_points[root_index + 1]


            # get distances between the extrema
            distances = {}
            start = 0
            for root in confirmed_roots:
                for i in range(start, len(confirmed_roots)):
                    key = list(confirmed_roots.keys())[i]
                    if key != root:
                        distances[f'{root}-{key}'] = abs(confirmed_roots[root] - confirmed_roots[key])
                start += 1

            a = [x / sorted(distances.values())[-1] for x in sorted(distances.values())]
            bins = {}
            bins[1] = [x for x in a if x < 0.9]
            bins[2] = [x for x in a if 0.9 < x]

            if len(bins[2]) == 1:

                # find pair of extrema that have maximal distance by Y between them
                left_extremum = None
                right_extremum = None
                dist = sys.float_info.min
                xs = list(confirmed_roots.keys())
                for first_root_index in range(0, len(xs)):
                    left_root_x = xs[first_root_index]
                    left_root_y = confirmed_roots[left_root_x]
                    for second_root_index in range(first_root_index + 1, len(xs)):
                        right_root_x = xs[second_root_index]
                        right_root_y = confirmed_roots[right_root_x]
                        current_distance = abs(right_root_y - left_root_y)
                        if current_distance > dist:
                            left_extremum = left_root_x
                            right_extremum = right_root_x
                            dist = current_distance
                left_extremum, right_extremum = right_extremum, left_extremum  # as our xs are somehow sorted in reverse order, flip them back

                if left_extremum is not None:
                    # cut the graph by left and right extrema
                    shortened_xs = []
                    shortened_ys = []
                    for element_index in range(0, len(answer_range[0]) - 1):
                        x = answer_range[0][element_index]
                        if left_extremum < x < right_extremum:
                            y = irft[element_index]
                            shortened_xs.append(x)
                            shortened_ys.append(y)

                    # find R^2 value for 3rd polynome approximation
                    coeffs = np.polyfit(shortened_xs, shortened_ys, 3)
                    # r-squared
                    p = np.poly1d(coeffs)
                    # fit values, and mean
                    yhat = p(shortened_xs)  # or [p(z) for z in x]
                    ybar = np.sum(shortened_ys) / len(shortened_ys)
                    ssreg = np.sum((yhat - ybar) ** 2)
                    sstot = np.sum((shortened_ys - ybar) ** 2)
                    r_squar = ssreg / sstot
                    if r_squar > 0.98:
                        return True
                    else:
                        return False
            else:
                return False
        else:
            return False
    else:
        return False


def check_the_spectra(corrected_spectrum: list, differential_spectrum: DataFrame):

    # define the spectra edges
    first =[x for x in differential_spectrum['Wavelength']].index(380)
    last =[x for x in differential_spectrum['Wavelength']].index(450)

    # rotate the answer range and get extrema
    answer_range = np.array(
        [[float(w) for w in [x for x in differential_spectrum['Wavelength']][first:last]], [b for b in corrected_spectrum[first:last]]])

    return get_the_answer(answer_range)