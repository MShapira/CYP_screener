import tkinter


def set_selected_well_as_protein_well(button):
    plate = button.plate
    well_position = button.label
    plate.set_proten_well_position(well_position)


# this function is called after right click on well button
def show_well_button_popup(event):
    button = event.widget

    type_labels = ['blank / protein', 'inhibitor', 'substrate', 'type I', 'reverse type I', 'type II', 'reverse type II(?)', 'sample']

	# construct type variable that will hold well type and share it across radiobuttons
    type_variable = tkinter.IntVar()
    type_variable.set(3)

    # construct pop-up menu
    popup = tkinter.Menu(window, tearoff=0)
    popup.add_command(label='set as protein / blank well', command=lambda: set_selected_well_as_protein_well(button))
    popup.add_separator()
    popup.add_command(label='set name...')
    popup.add_separator()
    popup.add_command(label='save graph...')
    popup.add_separator()
    popup.add_command(label='inactivate')

    # show pop-up
    try:
        popup.tk_popup(event.x_root, event.y_root, 0)
    finally:
        popup.grab_release()


# initialize the sample button
def initialize_button(plate_tab, plate, button_text, row_index, column_index, state='normal'):
    button = tkinter.Button(plate_tab, width=5, text=button_text)
    button.config(state=state)
    button.bind('<ButtonRelease-3>', show_well_button_popup)
    button.label = button_text
    button.grid(row=row_index, column=column_index)
    button.plate = plate

    return button


window = tkinter.Tk()