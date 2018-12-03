from pandas import *
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from sklearn import linear_model


# read the input and generate out csv
def read_input(filename: str)->tuple:
    csv_files= []
    with open(filename, encoding="utf-16", errors='ignore') as file:
        lines = file.readlines()

        # check the number of plates
        blocks = int(lines[0].split(' ')[1])
        if blocks % 2 != 0:
            print('I have an odd number of plates!')

        for line in range(1, len(lines)):
            if lines[line].split('\t')[0] == 'Plate:':
                name = lines[line].split('\t')[1]
                csv_files.append(name)
                with open('{0}.csv'.format(name), 'w') as out:
                    line += 1
                    for i in range(line, len(lines)):
                        if lines[i].split('\t')[0].rstrip('\n') != '~End':
                            out.write(lines[i].replace('\t', ','))
                        else:
                            line = line + i
                            break
                    out.close()

        file.close()

        return tuple(csv_files)


# get the data as dictionary
def scv_parser(csv_filename: str)->pandas.core.frame.DataFrame:
    with open(csv_filename, newline='') as csvfile:
        data = pandas.read_csv(csvfile)
        csvfile.close()

    return data


# construct the final data lists
def calculate_final_data(first_plate: pandas.core.frame.DataFrame, second_plate: pandas.core.frame.DataFrame, protein_well: str, sample_well: str)->dict:

    # generate calculations to form final data set for graph. Data is collecting as a dict with wavelenght as a key and OD as a value
    final_data = {}
    for i in range(0, len(second_plate[sample_well].get_values())):
        mid_data = second_plate[sample_well].get_values()[i] - second_plate[protein_well].get_values()[i] + \
                   first_plate[protein_well].get_values()[i] - first_plate[sample_well].get_values()[i]
        final_data[first_plate['Wavelength'].get_values()[i]] = mid_data

    return final_data


# generate graph
def generate_graph(data: dict, sample_well: str):
    x = np.fromiter(data.keys(), dtype = float)
    y = np.fromiter(data.values(), dtype = float)
    # smooth the curve
    y_smoothed = savgol_filter(y, 15, 3)  # window size 15, polynomial order 3

    # build the linear regression
    X = pandas.DataFrame(data=x, index=np.array(range(0, len(x))))
    Y = pandas.DataFrame(data=y, index=np.array(range(0, len(y))))
    lm = linear_model.LinearRegression()
    model = lm.fit(X, Y)
    predictions = lm.predict(X)
    new_y = [x[0] for x in predictions]
    plt.plot(x, new_y, color='green')

    plt.plot(x, y, color='blue')
    plt.plot(x, y_smoothed, color='red')
    plt.xlabel('Wavelength, nm')
    plt.ylabel('Optical Density, A')
    plt.title(sample_well)
    plt.show()


file1, file2 = read_input(filename='test.txt')
first_plate = scv_parser('{0}.csv'.format(file1))
second_plate = scv_parser('{0}.csv'.format(file2))
final_data = calculate_final_data(first_plate, second_plate, protein_well='A1', sample_well='C12')
generate_graph(final_data, 'C12')