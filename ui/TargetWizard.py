import threading
import tkinter as tk
from tkinter import ttk

import ttkbootstrap

import meta
import util

main = ttk.Frame(meta.win)
main.pack(padx=16, pady=8)
var = tk.StringVar()
threading.Thread(target=lambda: util.show_headline(var, meta.server)).start()
ttk.Label(main, textvariable=var, justify="center").pack()
notebook = ttk.Notebook(main)
notebook.pack(fill="x")
util.add_console(main).pack(fill="x")
ttk.Label(main, text=meta.copyright, justify="center").pack()

import TargetSelect, TargetBind, TargetReport, TargetView, extra
notebook.add(TargetSelect.main, text="Target Selection")
notebook.add(TargetBind.main, text="Target Binding")
notebook.add(TargetReport.main, text="Report Generation")
notebook.add(TargetView.main, text="Visualization")
notebook.add(extra.main, text="Extra Configuration")


util.bind_exit(meta.win, [TargetSelect, TargetReport, TargetView])
util.center_window(meta.win)
tk.mainloop()
