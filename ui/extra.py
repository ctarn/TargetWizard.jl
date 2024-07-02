import os
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

V = meta.vars
path = os.path.join(meta.homedir, "extra.cfg")

util.init_form(main)
I = util.Counter()
util.add_entry(main, I.next(), "pLink Cfg. Directory:", V["cfg_pl"], "Select", util.askdir(V["cfg_pl"]))
util.add_entry(main, I.next(), "pFind Cfg. Directory:", V["cfg_pf"], "Select", util.askdir(V["cfg_pf"]))
util.add_entry(main, I.next(), "TargetSelect:", V["targetselect"], "Select", util.askfile(V["targetselect"]))
util.add_entry(main, I.next(), "TargetSelectXL:", V["targetselectxl"], "Select", util.askfile(V["targetselectxl"]))
util.add_entry(main, I.next(), "TargetBind:", V["targetbind"], "Select", util.askfile(V["targetbind"]))
util.add_entry(main, I.next(), "Generators:", V["generators"], "Select", util.askdir(V["generators"]))
util.add_entry(main, I.next(), "Viewers:", V["viewers"], "Select", util.askdir(V["viewers"]))
f = ttk.Frame(main)
f.grid(column=0, row=I.next(), columnspan=3)
ttk.Button(f, text="Load Configuration", command=lambda: util.load_task(path, V)).pack(side="left", padx=16, pady=8)
ttk.Button(f, text="Save Configuration", command=lambda: util.save_task(path, V)).pack(side="left", padx=16, pady=8)
util.load_task(path, V)
