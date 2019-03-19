import tkinter


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


# initialize the sample button
def initialize_button(plate_tab, button_text, row_index, column_index, state='enabled'):
    button = tkinter.Button(plate_tab, width=5, text=button_text)
    if state == 'disabled':
        button.config(state=state)
    button.bind('<ButtonRelease-3>', show_well_button_popup)
    button.label = button_text
    button.grid(row=row_index, column=column_index)

    return button


window = tkinter.Tk()