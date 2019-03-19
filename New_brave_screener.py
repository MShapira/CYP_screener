import tkinter
from tkinter import ttk
from Logic import *
from Classes import Well, Plate
from guitools import initialize_button, window


# run the logic on gui
def tuning_the_gui():
    experiments = parse_the_input()

    inner_popups_constructor(experiments)


# construct the inner popups
def inner_popups_constructor(data: dict = None):

    # 2. setup panel with wells
    try:
        window.nametowidget('well_panel')
        ttk.Frame(window, name='well_panel').destroy()
        well_panel = ttk.Frame(window, name='well_panel')
    except KeyError:
        well_panel = ttk.Frame(window, name='well_panel')
    well_panel.pack()

        # 2.1 tab control widget
    plate_tab_control = ttk.Notebook(well_panel)
    plate_tab_control.pack(expand=1, fill='both')
        # 2.2 tab for plates
            # 2.2.1 frame for plate tab
    if data:
        # generate as much plates as data size
        for plate_index in range(0, len(data)):
            plate_tab = ttk.Frame(plate_tab_control)
            plate_tab_control.add(plate_tab, text='Plate {0}'.format(plate_index + 1))
            plate = Plate(number=plate_index, raw_data=data[plate_index])
            plate.get_area()
            # 2.2.2 well buttons
            plate_well_buttons = []
            for row_index in range(len(row_captions)):
                row_buttons = []
                for column_index in range(column_count):
                    position = f'{row_captions[row_index]}{column_index+1}'
                    if position in plate.working_area:
                        well = Well(position=position, type='sample', raw_data=plate.calculate_difference_for_well(position))
                        plate.wells.append(well)
                    state = 'enabled'
                    if plate.working_area and position not in plate.working_area:
                        state = 'disabled'
                    button = initialize_button(plate_tab=plate_tab, button_text=position, row_index=row_index,
                                               column_index=column_index, state=state)
                    row_buttons.append(button)
                plate_well_buttons.append(row_buttons)
    else:
        plate_tab = ttk.Frame(plate_tab_control)
        plate_tab_control.add(plate_tab, text='Plate {0}'.format(1))
        # 2.2.2 well buttons
        plate_well_buttons = []
        for row_index in range(len(row_captions)):
            row_buttons = []
            for column_index in range(column_count):
                button_text = f'{row_captions[row_index]}{column_index+1}'
                button = initialize_button(plate_tab=plate_tab, button_text=button_text, row_index=row_index,
                                           column_index=column_index)
                row_buttons.append(button)
            plate_well_buttons.append(row_buttons)



if __name__ == '__main__':

    row_captions = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    column_count = 12

    # 1. construct window
    window = window
    window.title("Brave new screener")

    inner_popups_constructor()

    # 3. setup bottom panel with buttons
    bottom_panel = ttk.Frame(window, height=20)
    bottom_panel.pack(side='bottom')
        # 3.1 button for file selection
    select_file_button = ttk.Button(bottom_panel, width=20, text='select file...', command=tuning_the_gui)
    select_file_button.pack(side=tkinter.LEFT, padx=5, pady=5)
        # 3.2 RUN button
    run_button = ttk.Button(bottom_panel, width=20, text='run')
    run_button.pack(side=tkinter.RIGHT, padx=5)

	# N. bring window to life
    window.mainloop()
