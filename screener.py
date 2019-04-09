import io
import tkinter
import numpy
import matplotlib.pyplot as pyplot
from tkinter import ttk, filedialog, font
from pandas import concat, read_csv, DataFrame
from scipy import sparse
from scipy.signal import savgol_filter
from scipy.sparse.linalg import spsolve
from Spectrum_check import check_the_spectra


# TODO:
# - implement saving graphs for wells (is it still needed as a separate action? yes)
# - what about well types (like inhibitor, substrate, sample)?
# - what does 'inactivate' mean for well? if well is inactivated, do not compute answer for it
# - implement own text style for inhibitor well button


well_row_captions = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
well_column_count = 12


def construct_well_position(row_index: int, column_index: int) -> str:
    '''
        helper function converting indices like (2, 0) to text like 'C1'
    '''
    return f'{well_row_captions[row_index]}{column_index+1}'

class Well:
    def __init__(self, plate_name: str, position: str, first_half_spectrum: DataFrame, second_half_spectrum: DataFrame):
        self.plate_name = plate_name
        self.position = position
        self.first_half_spectrum = first_half_spectrum
        self.second_half_spectrum = second_half_spectrum
        self.differential_spectrum = None
        self.corrected_spectrum = None

    def construct_differential_spectrum(self, wavelengths: DataFrame, protein_well):
        # spectrum = (U2_x - U1_x) + (U1_p - U2_p)
        spectrum = self.second_half_spectrum - self.first_half_spectrum + (protein_well.first_half_spectrum - protein_well.second_half_spectrum)
        self.differential_spectrum = concat([wavelengths, spectrum], axis=1, keys=['Wavelength', 'Difference'])
        print(type(self.differential_spectrum))

    def construct_corrected_spectrum(self):

        # baseline correction
        # y - 1D array of y_coords, lam - lambda for smoothness (10^2 <= λ <= 10^9), p for asymmetry (0.001 <= p <= 0.1)
        def baseline_correction(y, lam: int = 10000, p: float = 0.01, niter: int = 15):
            L = len(y)
            D = sparse.csc_matrix(numpy.diff(numpy.eye(L), 2))
            w = numpy.ones(L)
            for i in range(niter):
                W = sparse.spdiags(w, 0, L, L)
                Z = W + lam * D.dot(D.transpose())
                z = spsolve(Z, w * y)
                w = p * (y > z) + (1 - p) * (y < z)

            return z

        baseline = baseline_correction(self.differential_spectrum['Difference'])
        corrected_spectrum = []
        for i in range(0, len(self.differential_spectrum['Difference'])):
            corrected_spectrum.append(self.differential_spectrum['Difference'][i] - baseline[i])
        self.corrected_spectrum = [x for x in savgol_filter(corrected_spectrum, 15, 3)]

    def show_differential_spectrum_graph(self):
        if self.differential_spectrum is None:
            print('WARNING: attempt to show differential spectrum graph for well that does not have differential spectrum')
            return
        
        x = self.differential_spectrum['Wavelength']
        y_raw = self.differential_spectrum['Difference']
        y_smoothed = self.corrected_spectrum[30:120]
        fig = pyplot.figure()
        pyplot.plot(x, y_raw, color='blue')
        pyplot.plot(x[30:120], y_smoothed, color='red')
        pyplot.xlabel('Wavelength, nm')
        pyplot.ylabel('Optical Density, A')
        pyplot.legend(['raw', 'smoothed'])
        title = f'{self.position} [{self.plate_name}]'  # TODO: pass here well name when implemented
        fig.canvas.set_window_title(title)
        pyplot.title(title)
        pyplot.show()


class PlateRawData:
    def __init__(self, index: int, first_half_header: str, first_half_data: DataFrame, second_half_header: str, second_half_data: DataFrame):
        self.index = index
        self.wavelengths = first_half_data['Wavelength']
        self.first_half_header = first_half_header
        self.first_half_data = first_half_data
        self.second_half_header = second_half_header
        self.second_half_data = second_half_data

    def get_combined_name(self) -> str:
        first_half_name = self.first_half_header.split('\tTimeFormat')[0].split('Plate:\t')[1]
        second_half_name = self.second_half_header.split('\tTimeFormat')[0].split('Plate:\t')[1]
        
        # cut strange right part with numbers
        first_half_name = first_half_name.split('\t')[0]
        second_half_name = second_half_name.split('\t')[0]

        return f'{first_half_name} / {second_half_name}'


class Plate:
    def __init__(self, raw_data: PlateRawData):
        self.raw_data = raw_data
        self.protein_well_position = None
        self.substrate_well_position = None
        self.wells = Plate.construct_wells(raw_data)

    @staticmethod
    def construct_wells(raw_data: PlateRawData) -> list:
        wells = list()

        for row_index in range(len(well_row_captions)):
            for column_index in range(well_column_count):
                well_position = construct_well_position(row_index, column_index)
                # skip well that is not presented in data
                if not well_position in raw_data.first_half_data.keys():
                    continue
                well = Well(raw_data.get_combined_name(), well_position, raw_data.first_half_data[well_position], raw_data.second_half_data[well_position])
                wells.append(well)

        return wells

    def get_well(self, position: str) -> Well:
        for well in self.wells:
            if well.position == position:
                return well
        
        # if we are here then well was not found
        print(f'WARNING: attempt to get well {position} for plate #{self.index+1} failed - no such well')
        return None


def get_vertical_width_for_spectrum(spectrum: list) -> float:
    lowest_value = min(spectrum[30:120])
    highest_value = max(spectrum[30:120])
    return highest_value - lowest_value


class Screener:
    def __init__(self):
        # setup main UI
        self.window = Screener.construct_window()
        self.construct_plate_tab_control()
        self.construct_bottom_panel()

        # add one empty plate just to have some visible buttons
        self.add_plate_tab()

        # setup data
        self.plates = list()

    def run(self):
        self.window.mainloop()

    @staticmethod
    def construct_window() -> tkinter.Tk:
        window = tkinter.Tk()
        window.title('[no file selected] - screener')
        return window

    def construct_plate_tab_control(self):
        # construct panel that will hold tab control widget and all tabs
        panel = ttk.Frame(self.window, name='plate_panel')
        panel.pack()

        tab_control = ttk.Notebook(panel, name='plate_tab_control')
        tab_control.pack(expand=1, fill='both')

    def add_plate_tab(self, plate: Plate = None):
        plate_tab_control = self.get_plate_tab_control()
        if plate is not None:
            # construct normal plate tab
            plate_tab = ttk.Frame(plate_tab_control, name=f'plate_{plate.raw_data.index}')
            plate_tab_control.add(plate_tab, text=plate.raw_data.get_combined_name())
            for row_index in range(len(well_row_captions)):
                for column_index in range(well_column_count):
                    button_position = construct_well_position(row_index, column_index)
                    disabled = button_position not in plate.raw_data.first_half_data.keys()
                    self.construct_well_button(plate_tab, row_index, column_index, disabled)
        else:
            # construct plate tab with disabled buttons only
            plate_tab = ttk.Frame(plate_tab_control, name='plate_0')
            plate_tab_control.add(plate_tab, text='[no data]')
            for row_index in range(len(well_row_captions)):
                for column_index in range(well_column_count):
                    self.construct_well_button(plate_tab, row_index, column_index)

    def construct_well_button(self, tab: ttk.Frame, row_index: int, column_index: int, disabled: bool = True):
        position = construct_well_position(row_index, column_index)
        name = f'button_{position}'
        button = tkinter.Button(tab, width=5, text=position, name=name)
        Screener.set_well_button_style(button, False)
        if disabled:
            button.config(state='disabled')
        else:
            button.bind('<ButtonRelease-1>', self.show_well_graph)
            button.bind('<ButtonRelease-3>', self.show_well_button_popup)
        button.grid(row=row_index, column=column_index)

    def show_well_graph(self, event):
        button = event.widget

        plate = self.get_plate_for_button(button)
        # we cannot show graph if protein well is not set yet for plate
        if plate.protein_well_position is None:
            print(f'WARNING: cannot show differential spectrum for well from plate #{plate.raw_data.index+1} - protein well is not set')
            return

        well = self.get_corresponding_well(button)
        # construct differential spectrum for well
        well.construct_differential_spectrum(plate.raw_data.wavelengths, plate.get_well(plate.protein_well_position))
        well.show_differential_spectrum_graph()

    def show_well_button_popup(self, event):
        button = event.widget
        window = button.master.master.master

        plate = self.get_plate_for_button(button)
        well = self.get_corresponding_well(button)

        # construct pop-up menu
        popup = tkinter.Menu(window, tearoff=0)
        # add 'set protein well' command only if corresponding well is not already set as protein well
        if plate.protein_well_position != well.position:
            popup.add_command(label='set as protein / blank well',
                              command=lambda: self.set_selected_well_as_protein_well(button))
        else:
            popup.add_command(label='[this well is already set as protein / blank well]', state='disabled')
        # add 'set substrate well' command only if corresponding well is not already set as substrate well
        if plate.substrate_well_position != well.position:
            popup.add_command(label='set as substrate / inhibitor well',
                              command=lambda: self.set_selected_well_as_substrate_well(button))
        else:
            popup.add_command(label='[this well is already set as substrate / inhibitor well]', state='disabled')
        popup.add_separator()
        popup.add_command(label='set name...', state='disabled')
        popup.add_separator()
        popup.add_command(label='save graph...', state='disabled')
        popup.add_separator()
        popup.add_command(label='inactivate', state='disabled')

        # show pop-up
        try:
            popup.tk_popup(event.x_root, event.y_root, 0)
        finally:
            popup.grab_release()

    def get_corresponding_well(self, button: ttk.Button) -> Well:
        # find corresponding plate
        plate = self.get_plate_for_button(button)

        # then find well
        well_position = Screener.get_position_for_button(button)
        return plate.get_well(well_position)

    def set_selected_well_as_protein_well(self, button: ttk.Button):
        plate = self.get_plate_for_button(button)

        # un-highlight previous protein well button for active plate (if any)
        if plate.protein_well_position is not None:
            previous_protein_well_button = self.get_button_for_well(plate.raw_data.index, plate.protein_well_position)
            self.set_well_button_style(previous_protein_well_button, False)

        # store protein well position for active plate
        plate.protein_well_position = Screener.get_position_for_button(button)
        print(f'plate #{plate.raw_data.index+1}: set protein well position as {plate.protein_well_position}')

        # highlight new protein well button
        self.set_well_button_style(button, True)

        # when any well marked as protein well, enable 'run' button
        self.get_run_button().configure(state='normal')

    def set_selected_well_as_substrate_well(self, button: ttk.Button):
        plate = self.get_plate_for_button(button)

        # un-highlight previous substrate well button for active plate (if any)
        if plate.substrate_well_position is not None:
            previous_substrate_well_button = self.get_button_for_well(plate.raw_data.index, plate.substrate_well_position)
            self.set_well_button_style(previous_substrate_well_button, False)

        # store substrate well position for active plate
        plate.substrate_well_position = Screener.get_position_for_button(button)
        print(f'plate #{plate.raw_data.index+1}: set substrate well position as {plate.substrate_well_position}')

        # highlight new substrate well button
        self.set_well_button_style(button, True)

    def get_plate_for_button(self, button: ttk.Button) -> Plate:
        plate_tab_name = button.winfo_parent().split('.')[-1]
        plate_index = int(plate_tab_name.split('_')[1])
        return self.plates[plate_index]

    @staticmethod
    def get_position_for_button(button: ttk.Button) -> str:
        return button.config('text')[-1]  # that's just a kind of TKinter magic, so trust me

    @staticmethod
    def set_well_button_style(button: ttk.Button, highlighted: bool):
        if highlighted:
            button.configure(font=font.Font(size=10, weight='bold', slant='roman', underline=1))
        else:
            button.configure(font=font.Font(size=10, weight='normal', slant='italic'))

    def construct_bottom_panel(self):
        panel = ttk.Frame(self.window, height=20, name='bottom_panel')
        panel.pack(side='bottom')

        # 'select file...' button
        select_file_button = ttk.Button(panel, width=20, text='select file...', name='select_file_button', command=self.open_data_file)
        select_file_button.pack(side=tkinter.LEFT, padx=5, pady=5)

        # 'run' button
        run_button = ttk.Button(panel, width=20, text='run', name='run_button', state='disabled', command=self.calculate_well_answers)
        run_button.pack(side=tkinter.RIGHT, padx=5)

    def get_plate_tab_control(self) -> ttk.Notebook:
        return self.window.children['plate_panel'].children['plate_tab_control']

    def get_button_for_well(self, plate_index: int, well_position: str) -> ttk.Button:
        plate_tab_control = self.get_plate_tab_control()

        # find corresponding plate tab at first
        plate_tab_name = f'plate_{plate_index}'
        plate_tab = plate_tab_control.children[plate_tab_name]

        # then find the button
        button_name = f'button_{well_position}'
        return plate_tab.children[button_name]

    def open_data_file(self):
        # let user select file
        data_filename = filedialog.askopenfilename(filetypes=([('all files', '*.*'), ('text files', '*.txt')]))
        # do nothing if no file selected
        if data_filename == () or data_filename == '':
            return

        self.window.title(f'[{data_filename}] - screener')

        # load its contents
        with open(data_filename, encoding='utf-16') as file:
            lines = file.readlines()

            # get number of plates
            plate_count = int(lines[0].split('=')[1])
            if plate_count % 2 != 0:
                print(f'WARNING: file has odd number of plates - {plate_count}')
            lines = lines[1:]

            # load raw plate lines at first
            plate_lines = list()
            current_plate_lines = list()
            for line in lines:
                if line.rstrip().startswith('~End'):
                    plate_lines.append(current_plate_lines.copy())
                    current_plate_lines = list()
                else:
                    current_plate_lines.append(line)

            # construct data frames for plates
            plates_raw_data = list()
            for plate_index in range(plate_count // 2):
                # extract header for first half                
                first_half_lines = plate_lines[plate_index * 2 + 0]
                first_half_header = first_half_lines[0][:-1]
                
                # construct data for first half
                first_half_data_line_buffer = io.StringIO('\n'.join(first_half_lines[1:]))
                first_half_data = read_csv(first_half_data_line_buffer, sep='\t', header=0)
                # remove columns with NaN
                first_half_data = first_half_data.dropna(how='all', axis=1)
                # remove temperature column
                first_half_data = first_half_data.drop('Temperature(¡C)', axis=1)
                
                # extract header for second half               
                second_half_lines = plate_lines[plate_index * 2 + 1]
                second_half_header = second_half_lines[0][:-1]

                # construct data for second half
                second_half_data_line_buffer = io.StringIO('\n'.join(second_half_lines[1:]))
                second_half_data = read_csv(second_half_data_line_buffer, sep='\t', header=0)
                # remove columns with NaN
                second_half_data = second_half_data.dropna(how='all', axis=1)
                # remove temperature column
                second_half_data = second_half_data.drop('Temperature(¡C)', axis=1)

                raw_data = PlateRawData(plate_index, first_half_header, first_half_data,
                                        second_half_header, second_half_data)
                plates_raw_data.append(raw_data)

        self.plates = [Plate(raw_data) for raw_data in plates_raw_data]

        # update the UI
        plate_tab_control = self.get_plate_tab_control()
        # remove all existing plate tabs
        for i in range(len(plate_tab_control.tabs())):
            plate_tab_control.forget(0)
        # add plate tabs with data
        for plate in self.plates:
            self.add_plate_tab(plate)
        # disable 'run' button
        self.get_run_button().configure(state='disabled')

    def get_run_button(self):
        return self.window.children['bottom_panel'].children['run_button']

    def calculate_well_answers(self):
        for plate in self.plates:
            # skip plates where protein well position is not set
            if plate.protein_well_position is None:
                continue

            # process wells
            protein_well = plate.get_well(plate.protein_well_position)
            min_vertical_width = None
            max_vertical_width = None
            for well in plate.wells:
                # skip protein well
                if well.position == plate.protein_well_position:
                    continue

                well.construct_differential_spectrum(plate.raw_data.wavelengths, protein_well)
                well.construct_corrected_spectrum()

                # TODO: run logic from Friday to get answer for well
                well.answer = check_the_spectra(well.corrected_spectrum, well.differential_spectrum)
                well.vertical_width = get_vertical_width_for_spectrum(well.corrected_spectrum)

                if min_vertical_width is None:
                    min_vertical_width = well.vertical_width
                else:
                    min_vertical_width = min(min_vertical_width, well.vertical_width)

                if plate.substrate_well_position is None:
                    if max_vertical_width is None:
                        max_vertical_width = well.vertical_width
                    else:
                        max_vertical_width = max(max_vertical_width, well.vertical_width)

            if plate.substrate_well_position is not None:
                substrate_well = [w for w in plate.wells if w.position == plate.substrate_well_position][0]
                max_vertical_width = substrate_well.vertical_width

            # update correpsonding UI
            for well in plate.wells:
                if well.position == plate.protein_well_position:
                    self.update_well_button_appearance(plate.raw_data.index, well.position, None, None)
                else:
                    relative_vertical_width = (well.vertical_width - min_vertical_width) / (max_vertical_width - min_vertical_width)
                    self.update_well_button_appearance(plate.raw_data.index, well.position, well.answer, relative_vertical_width)

    def update_well_button_appearance(self, plate_index: int, well_position: str, answer: bool, relative_vertical_width: float):
        button = self.get_button_for_well(plate_index, well_position)
        if answer is None:
            button['bg'] = 'grey'
        elif answer:
            if relative_vertical_width < 0.2:
                button['bg'] = 'grey'
            elif relative_vertical_width < 0.6:
                button['bg'] = 'yellow'
            elif relative_vertical_width < 0.8:
                button['bg'] = '#77ff77'
            else:
                button['bg'] = 'green'
        else:
            button['bg'] = 'red'


if __name__ == '__main__':
    screener = Screener()
    screener.run()
