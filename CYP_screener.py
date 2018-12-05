from pandas import *
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from sklearn import linear_model
import json
import scipy
import sys


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

#baseline calculation by Asymmetric Least Squares Smoothing
def baseline_als(y, lam=1e+7, p=0.001, niter=10):
  L = len(y)
  D = scipy.sparse.csc_matrix(np.diff(np.eye(L), 2))
  w = np.ones(L)
  for i in range(niter):
    W = scipy.sparse.spdiags(w, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = scipy.sparse.linalg.spsolve(Z, w*y)
    w = p * (y > z) + (1-p) * (y < z)
    
  return z

# generate graph
def generate_graph(data: dict, sample_well: str, prot: str):
    x = np.fromiter(data.keys(), dtype = float)
    y = np.fromiter(data.values(), dtype = float)
    # smooth the curve
    y_smoothed = savgol_filter(y, 15, 3)  # window size 15, polynomial order 3
    bsl = baseline_als(y_smoothed)
    fig = plt.figure()
#    plt.plot(x, y, color='blue')
#    plt.plot(x, y_smoothed, color='red')
    plt.plot(x, y_smoothed - bsl, color='green')
    plt.xlabel('Wavelength, nm')
    plt.ylabel('Optical Density, A')
    plt.title('{0}_{1}'.format(prot, sample_well))
    plt.show()
    fig.savefig('{0}_{1}.png'.format(prot, sample_well))

def get_cells(cellRange): #obtain cells sequence list [A3, A4, A5 .... H10, H11]
    cells = []
    if cellRange[0][0] == cellRange[1][0]:
        for num in range(int(cellRange[0][1:]), int(cellRange[1][1:])+1):
            cells.append(cellRange[0][0]+str(num))
    else:
        for letter in range(ord(cellRange[0][0]), ord(cellRange[1][0]) + 1):
            if letter == ord(cellRange[0][0]):
                for num in range(int(cellRange[0][1:]), 13):
                    cells.append(cellRange[0][0]+str(num))
            elif letter == ord(cellRange[1][0]):
                for num in range(1, int(cellRange[1][1:])+1):
                    cells.append(cellRange[1][0]+str(num))
            else:
                for num in range(1, 13):
                    cells.append(chr(letter)+str(num))
    return cells

parmFile = 'input.json' #sys.argv[2] #json parameters file
inputFile = 'CYP_7A_7B_19_21_LjG.txt' #sys.argv[1] #input file from Spectramax

files = read_input(filename = inputFile) #in one file it can be not only 2 plates, but 4, 6, 8 etc.
with open(parmFile) as f:
    inputData = json.load(f)
for protein in inputData.keys():
    file1 = inputData[protein]['first_measurement'] #plate with first measurement
    file2 = inputData[protein]['second_measurement'] #plate with second measurement
    blank = inputData[protein]['blank'] #blank cell
    cellRan = inputData[protein]['range'].split('-')
    cells = get_cells(cellRan)
    if inputData[protein]['substrate']:
        substrate = inputData[protein]['substrate'] #substrate cell (if any)
    if inputData[protein]['inhibitor']:
        inhibitor = inputData[protein]['inhibitor'] #inhibitor cell (if any)
     # if input file has compound names, graphs will be with compounds names, in other case with cells names 
    if inputData[protein]['compounds']:
        compounds = {inputData[protein]['compounds'][i]:cells[i] for i in range(len(cells))}
    else:
        compounds = {cells[i]:cells[i] for i in range(len(cells))}
    if substrate:
        compounds['substrate'] = substrate
    if inhibitor:
        compounds['inhibitor'] = inhibitor

    first_plate = scv_parser('{0}.csv'.format(file1))
    second_plate = scv_parser('{0}.csv'.format(file2))
    
    comp = "LjG 40-58" #current compound
    final_data = calculate_final_data(first_plate, second_plate, protein_well = blank, sample_well=compounds.get(comp))

    generate_graph(final_data, comp, protein)