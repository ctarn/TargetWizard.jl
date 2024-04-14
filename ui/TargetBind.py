import os
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

PrecModeIC = "by isolation center"
PrecModeIW = "by isolation window"
PrecModeExIW = "by extended isolation window"
PrecModes = {PrecModeIC: "center", PrecModeIW: "window", PrecModeExIW: "extended_window"}
fmts = ["csv", "tsv", "ms2", "mgf", "pf2"]
vars_spec = {
    "data": {"type": tk.StringVar, "value": ""},
    "target": {"type": tk.StringVar, "value": ""},
    "out": {"type": tk.StringVar, "value": ""},
    "mode": {"type": tk.StringVar, "value": PrecModeExIW},
    "error_rt": {"type": tk.StringVar, "value": "16"},
    "error_mz": {"type": tk.StringVar, "value": "10.0"},
}
for fmt in fmts: vars_spec[f"fmt_{fmt}"] = {"type": tk.IntVar, "value": fmt in ["csv"]}
task = util.Task("TargetBind", vars_spec, path=meta.homedir, shared_vars_spec=meta.vars_spec, shared_vars=meta.vars)
V = task.vars

def run():
    task.call(V["targetbind"].get(), *(V["data"].get().split(";")),
        "--target", V["target"].get(),
        "--out", V["out"].get(),
        "--mode", PrecModes[V["mode"].get()],
        "--error_rt", V["error_rt"].get(),
        "--error_mz", V["error_mz"].get(),
        "--fmt", ",".join(filter(lambda x: V[f"fmt_{x}"].get(), fmts)),
    )

util.init_form(main)
I = 0
t = (("UMS", "*.ums"), ("MS2", "*.ms2"), ("All", "*.*"))
util.add_entry(main, I, "Targeted MS Data:", V["data"], "Select", util.askfiles(V["data"], V["out"], filetypes=t))
I += 1
t = (("Target List", "*.csv"), ("All", "*.*"))
util.add_entry(main, I, "Target List:", V["target"], "Select", util.askfile(V["target"], filetypes=t))
I += 1
c = ttk.Combobox(main, textvariable=V["mode"], values=list(PrecModes.keys()), state="readonly", justify="center")
util.add_entry(main, I, "Binding Mode:", c)
I += 1
util.add_entry(main, I, "Max. RTime Error:", V["error_rt"], "sec")
I += 1
util.add_entry(main, I, "Max. Mass Error:", V["error_mz"], "ppm")
I += 1
_, f, _ = util.add_entry(main, I, "Output Format", ttk.Frame(main))
for x in fmts: ttk.Checkbutton(f, text=x.upper(), variable=V[f"fmt_{x}"]).pack(side="left", expand=True)
I += 1
util.add_entry(main, I, "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
I += 1
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I, columnspan=3)
