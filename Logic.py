from pandas import *
import io
from Classes import Well
from tkinter import filedialog
from pprint import pprint


# parsing the file into plates
def parse_the_input()->dict:
    # open and decode the file
    data_filename = select_file()
    with open(data_filename, encoding="utf-16", errors='ignore') as file:
        lines = file.readlines()

        # make a dict of arrays by way of plates
        plates = dict()
        index = 1

        # get the number of plates
        number_of_plates = int(lines[0].split('=')[1])
        if number_of_plates % 2 != 0:
            print('I have an odd number of plates!')

        # find the end of each plate and divide file on them
        for line in lines[1:]:
            if line.rstrip().split('\t')[0] != '~End':
                if not index in plates:
                    plates[index] = []
                plates[index].append(line)
            else:
                index += 1

        # get rid of additional "plate"
        if len(plates) != number_of_plates:
            plates.pop(index)

    # create a final dict "experiments" that contains a couples of plates
    experiments = dict()
    dict_index = 0
    for plate_index in range(1, len(plates)+1):
        if dict_index not in experiments.keys():
            experiments[dict_index] = []
        df = read_csv(io.StringIO('\n'.join(plates[plate_index])), sep='\t', header=1)

        # delete all columns with NaN and Temperature column
        df = df.dropna(how='all', axis=1)
        if 'Temperature(¡C)' in list(df.columns):
            df = df.drop('Temperature(¡C)', axis=1)

        experiments[dict_index].append(df)
        # pprint(experiments[dict_index])

        if plate_index % 2 == 0:
            dict_index += 1

    pprint(type(experiments[0][1]))

    return experiments


# construct the list of wells with calculated spectra
def get_filled_wells(experiment: list, protein_well: str, area: list)->list:

    wells = []

    for position in area:
        if position != protein_well:

            # generate calculations to form final data set for graph. Data is collecting as a dict
            # with Wavelength as a key and OD as a value
            spectrum = {}
            for i in range(0, len(experiment[1][position].get_values())):
                mid_data = experiment[1][position].get_values()[i] - experiment[1][protein_well].get_values()[i] + \
                           experiment[0][protein_well].get_values()[i] - experiment[0][position].get_values()[i]
                spectrum[experiment[0]['Wavelength'].get_values()[i]] = mid_data

            well = Well(position=position, spectrum=spectrum)
            wells.append(well)
        else:
            well = Well(position=position, type='blank', name='Protein')
            wells.append(well)

    return wells


# get the name of selected file
def select_file():
    file_name = filedialog.askopenfilename(filetypes=([('All files', '*.*'), ('Text files', '*.txt')]))

    return file_name


# initialize the button
def initialize_the_button():
    button = tkinter.Button(plate_tab, width=5, text=button_text)
    button.bind('<ButtonRelease-3>', show_well_button_popup)
    button.label = button_text
    button.grid(row=row_index, column=column_index)
