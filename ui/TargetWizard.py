import os
import tkinter as tk
from tkinter import ttk

import ttkbootstrap

import meta
import util

os.makedirs(meta.homedir, exist_ok=True)

win = tk.Tk()
win.title(meta.name)
win.iconphoto(True, tk.PhotoImage(file=util.get_content(f"{meta.name}.png", shared=True)))
win.resizable(False, False)

main = ttk.Frame(win)
main.pack(padx=16, pady=8)
util.add_headline(main, meta.server)[0].pack()
notebook = ttk.Notebook(main)
notebook.pack(fill="x")
util.add_console(main).pack(fill="x")
ttk.Label(main, text=meta.copyright, justify="center").pack()

import TargetSelect, TargetReport, TargetView
notebook.add(TargetSelect.main, text="Target Selection")
notebook.add(TargetReport.main, text="Report Generation")
notebook.add(TargetView.main, text="Visualization")

util.bind_exit(win, [TargetSelect, TargetReport, TargetView])
util.center_window(win)
tk.mainloop()
