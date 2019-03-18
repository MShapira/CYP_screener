import tkinter
from tkinter import ttk, filedialog
from Logic import *
from pprint import pprint


# this function is called after right click on well button
def show_well_button_popup(event):
    button = event.widget
    print(button.label) # as you can see, we can access to button that was pressed

    type_labels = ['blank / protein', 'inhibitor', 'substrate', 'type I', 'reverse type I', 'type II', 'reverse type II(?)', 'sample']

	# construct type variable that will hold well type and share it across radiobuttons
    type_variable = tkinter.IntVar()
    type_variable.set(3)

    # construct type sub-menu
    type_menu = tkinter.Menu(window, tearoff=0)
    for i in range(len(type_labels)):
        type_menu.add_radiobutton(label=type_labels[i], value=i, variable=type_variable)

    # construct pop-up menu
    popup = tkinter.Menu(window, tearoff=0)
    popup.add_command(label='set name...')
    popup.add_cascade(label='type', menu=type_menu)
    popup.add_separator()
    popup.add_command(label='save graph...')
    popup.add_separator()
    popup.add_command(label='inactivate')

    # show pop-up
    try:
        popup.tk_popup(event.x_root, event.y_root, 0)
    finally:
        popup.grab_release()


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
        for plate in data:
            print(plate)
            working_area = get_area(data[plate])
            plate_tab = ttk.Frame(plate_tab_control)
            plate_tab_control.add(plate_tab, text='Plate {0}'.format(plate + 1))
                # 2.2.2 well buttons
            plate_well_buttons = []
            for row_index in range(len(row_captions)):
                row_buttons = []
                for column_index in range(column_count):
                    button_text = f'{row_captions[row_index]}{column_index+1}'
                    button = tkinter.Button(plate_tab, width=5, text=button_text)
                    if working_area and button_text not in working_area:
                        button.config(state='disabled')
                    button.bind('<ButtonRelease-3>', show_well_button_popup)
                    button.label = button_text
                    button.grid(row=row_index, column=column_index)
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
                button = tkinter.Button(plate_tab, width=5, text=button_text)
                button.bind('<ButtonRelease-3>', show_well_button_popup)
                button.label = button_text
                button.grid(row=row_index, column=column_index)
                row_buttons.append(button)
            plate_well_buttons.append(row_buttons)



if __name__ == '__main__':

    row_captions = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    column_count = 12

    # 1. construct window
    window = tkinter.Tk()  
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
