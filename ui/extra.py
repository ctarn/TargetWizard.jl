import os
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

V = meta.vars
path = os.path.join(meta.homedir, "extra.cfg")

util.init_form(main)
I = 0
util.add_entry(main, I, "pLink Cfg. Directory:", V["cfg_pl"], "Select", util.askdir(V["cfg_pl"]))
I += 1
util.add_entry(main, I, "pFind Cfg. Directory:", V["cfg_pf"], "Select", util.askdir(V["cfg_pf"]))
I += 1
util.add_entry(main, I, "TargetSelect:", V["targetselect"], "Select", util.askfile(V["targetselect"]))
I += 1
util.add_entry(main, I, "TargetSelectXL:", V["targetselectxl"], "Select", util.askfile(V["targetselectxl"]))
I += 1
util.add_entry(main, I, "TargetBind:", V["targetbind"], "Select", util.askfile(V["targetbind"]))
I += 1
util.add_entry(main, I, "Generators:", V["generators"], "Select", util.askdir(V["generators"]))
I += 1
util.add_entry(main, I, "Viewers:", V["viewers"], "Select", util.askdir(V["viewers"]))
I += 1
f = ttk.Frame(main)
f.grid(column=0, row=I, columnspan=3)
ttk.Button(f, text="Load Configuration", command=lambda: util.load_task(path, V)).pack(side="left", padx=16, pady=8)
ttk.Button(f, text="Save Configuration", command=lambda: util.save_task(path, V)).pack(side="left", padx=16, pady=8)
util.load_task(path, V)
