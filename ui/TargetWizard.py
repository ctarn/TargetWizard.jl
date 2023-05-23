import os
import sys
import threading
import tkinter as tk
from tkinter import ttk, messagebox

import ttkbootstrap

import meta
import util

os.makedirs(meta.homedir, exist_ok=True)

win = tk.Tk()
win.title(meta.name)
win.iconphoto(True, tk.PhotoImage(file=util.get_content(f"{meta.name}.png", shared=True)))
win.resizable(False, False)
util.center_window(win)

main = ttk.Frame(win)
main.pack(padx=16, pady=8)

headline = tk.StringVar()
ttk.Label(main, textvariable=headline, justify="center").pack()

notebook = ttk.Notebook(main)
notebook.pack(fill="x")

console = tk.Text(main, height=12, state="disabled")
console.pack(fill="x")

ttk.Label(main, text=meta.copyright, justify="center").pack()

sys.stdout = util.Console(console)
sys.stderr = util.Console(console)

threading.Thread(target=lambda: util.show_headline(headline, meta.server)).start()

import TargetSelect
notebook.add(TargetSelect.main, text="Target Selection")

import TargetReport
notebook.add(TargetReport.main, text="Report Generation")

import TargetView
notebook.add(TargetView.main, text="Visualization")

mods = [TargetSelect, TargetReport]

def on_exit():
    if (not any([m.running for m in mods]) or
        messagebox.askokcancel("Quit", "Task running. Quit now?")):
        [m.do_stop() for m in mods]
        win.destroy()

win.protocol("WM_DELETE_WINDOW", on_exit)

util.center_window(win)

tk.mainloop()
