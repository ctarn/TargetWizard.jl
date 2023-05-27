import os
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

if util.is_darwin:
    path_mono = "/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono"
else:
    path_mono = "mono"

tds = ["TT", "TD", "DD"]
pts = ["Inter", "Intra"]
fmts = ["TW", "TmQE", "TmFu"]
fmts_name = ["TargetWizard", "Thermo Q Exactive", "Thermo Fusion"]
vars_spec = {
    "data": {"type": tk.StringVar, "value": ""},
    "psm": {"type": tk.StringVar, "value": ""},
    "name": {"type": tk.StringVar, "value": "TargetWizard"},
    "fdr_min": {"type": tk.StringVar, "value": "-Inf"},
    "fdr_max": {"type": tk.StringVar, "value": "Inf"},
    "fdr_ge": {"type": tk.StringVar, "value": "≤"},
    "fdr_le": {"type": tk.StringVar, "value": "≤"},
    "batch": {"type": tk.StringVar, "value": "Inf"},
    "rtime": {"type": tk.StringVar, "value": "240"},
    "lc": {"type": tk.StringVar, "value": "Inf"},
    "out": {"type": tk.StringVar, "value": ""},
    "targetselect": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetSelect")},
    "thermorawread": {"type": tk.StringVar, "value": util.get_content("ThermoRawRead", "ThermoRawRead.exe", shared=True)},
    "mono": {"type": tk.StringVar, "value": path_mono},
}
for t in tds: vars_spec[f"td_{t}"] = {"type": tk.IntVar, "value": 1}
for t in pts: vars_spec[f"pt_{t}"] = {"type": tk.IntVar, "value": 1}
for t in fmts: vars_spec[f"fmt_{t}"] = {"type": tk.IntVar, "value": t == "TW"}
task = util.Task("TargetSelect", vars_spec, path=meta.homedir)
V = task.vars

def run_thermorawread(data, out):
    task.call(*([] if util.is_windows else [V["mono"].get()]), V["thermorawread"].get(), data, out)
    return os.path.join(out, os.path.splitext(os.path.basename(data))[0] + ".ms2")

def run():
    paths = []
    for p in V["data"].get().split(";"):
        ext = os.path.splitext(p)[1].lower()
        if ext == ".raw":
            p = run_thermorawread(p, V["out"].get())
        paths.append(p)
    task.call(V["targetselect"].get(), *paths, "--out", V["out"].get(),
        "--name", V["name"].get(),
        "--psm", V["psm"].get(),
        "--fdr_min", V["fdr_min"].get(),
        "--fdr_max", V["fdr_max"].get(),
        *(["--fdr_ge"] if V["fdr_ge"].get() == "≤" else []),
        *(["--fdr_le"] if V["fdr_le"].get() == "≤" else []),
        "--td", ",".join([t for t in tds if V[f"td_{t}"].get()]),
        "--pt", ",".join([t for t in pts if V[f"pt_{t}"].get()]),
        "--batch", V["batch"].get(),
        "--rtime", V["rtime"].get(),
        "--lc", V["lc"].get(),
        "--fmt", ",".join([t for t in fmts if V[f"fmt_{t}"].get()]),
    )

util.init_form(main)
I = 0
t = (("MS2", "*.ms2"), ("RAW", "*.raw"), ("All", "*.*"))
util.add_entry(main, I, "Data:", V["data"], "Select", util.askfiles(V["data"], V["out"], filetypes=t))
t = (("PSM", "*.csv"), ("All", "*.*"))
util.add_entry(main, I, "PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
util.add_entry(main, I, "Task Name:", V["name"])
I += 1
_, f, _ = util.add_entry(main, I, "FDR Range:", ttk.Frame(main), "%")
ttk.Entry(f, textvariable=V["fdr_min"]).pack(side="left", fill="x", expand=True)
ttk.Label(f, text="% ").pack(side="left")
ttk.Combobox(f, textvariable=V["fdr_ge"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Label(f, text=" FDR ").pack(side="left")
ttk.Combobox(f, textvariable=V["fdr_le"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Entry(f, textvariable=V["fdr_max"]).pack(side="left", fill="x", expand=True)
I += 1
_, f, _ = util.add_entry(main, I, "Target / Decoy Type:", ttk.Frame(main))
for t in tds: ttk.Checkbutton(f, text=t, variable=V[f"td_{t}"]).pack(side="left", expand=True)
I += 1
_, f, _ = util.add_entry(main, I, "Inter- / Intra-Protein:", ttk.Frame(main))
for t in pts: ttk.Checkbutton(f, text=t, variable=V[f"pt_{t}"]).pack(side="left", expand=True)
I += 1
util.add_entry(main, I, "Batch Size:", V["batch"])
I += 1
util.add_entry(main, I, "RT Window:", V["rtime"], "sec")
I += 1
util.add_entry(main, I, "LC Length:", V["lc"], "min")
I += 1
_, f, _ = util.add_entry(main, I, "List Format:", ttk.Frame(main))
for t, n in zip(fmts, fmts_name): ttk.Checkbutton(f, text=n, variable=V[f"fmt_{t}"]).pack(side="left", expand=True)
I += 1
util.add_entry(main, I, "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
I += 1
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I, columnspan=3)
I += 1
ttk.Separator(main, orient=tk.HORIZONTAL).grid(column=0, row=I, columnspan=3, sticky="EW")
ttk.Label(main, text="Advanced Configuration").grid(column=0, row=I, columnspan=3)
I += 1
util.add_entry(main, I, "TargetSelect:", V["targetselect"], "Select", util.askfile(V["targetselect"]))
I += 1
util.add_entry(main, I, "ThermoRawRead:", V["thermorawread"], "Select", util.askfile(V["thermorawread"]))
I += 1
if not util.is_windows:
    util.add_entry(main, I, "Mono Runtime:", V["mono"], "Select", util.askfile(V["mono"]))
    I += 1
