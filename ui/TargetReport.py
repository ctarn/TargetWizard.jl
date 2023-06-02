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

BAReport = "Basic Aquisition Report"
TSReport = "Target Selection Report"
TAReport = "Target Aquisition Report"
IDReport = "Identification Report"
XPCReport = "Cross-linked Peptide Coverage Report"
reports = [BAReport, TSReport, TAReport, IDReport, XPCReport]
target_fmts = {"Auto Detect": "auto", "TargetWizard": "TW", "Thermo Q Exactive": "TmQE", "Thermo Fusion": "TmFu"}
ion_types = ["a", "b", "c", "x", "y", "z", "a_NH3", "b_NH3", "y_NH3", "a_H2O", "b_H2O", "y_H2O"]
ion_names = ["a", "b", "c", "x", "y", "z", "a-NH₃", "b-NH₃", "y-NH₃", "a-H₂O", "b-H₂O", "y-H₂O"]
vars_spec = {
    "report": {"type": tk.StringVar, "value": BAReport},
    "ms": {"type": tk.StringVar, "value": ""},
    "error2": {"type": tk.StringVar, "value": "20.0"},
    "target": {"type": tk.StringVar, "value": ""},
    "target_fmt": {"type": tk.StringVar, "value": "Auto Detect"},
    "psm": {"type": tk.StringVar, "value": ""},
    "linker": {"type": tk.StringVar, "value": "DSSO"},
    "cfg_pl": {"type": tk.StringVar, "value": ""},
    "out": {"type": tk.StringVar, "value": ""},
    "generators": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
    "thermorawread": {"type": tk.StringVar, "value": util.get_content("ThermoRawRead", "ThermoRawRead.exe", shared=True)},
    "mono": {"type": tk.StringVar, "value": path_mono},
}
for t in ion_types: vars_spec[f"ion_{t}"] = {"type": tk.IntVar, "value": 1}
task = util.Task("TargetReport", vars_spec, path=meta.homedir)
V = task.vars

def run_basicaquisitionreport():
    task.call(os.path.join(V["generators"].get(), "BasicAquisitionReport"), *(V["ms"].get().split(";")), "--out", V["out"].get())

def run_targetselectionreport():
    task.call(os.path.join(V["generators"].get(), "TargetSelectionReport"),
        *(V["target"].get().split(";")), "--out", V["out"].get(),
        "--fmt", target_fmts[V["target_fmt"].get()],
    )

def run_targetaquisitionreport():
    pass

def run_identificationreport():
    pass

def run_crosslinkedpeptidecoveragereport():
    task.call(os.path.join(V["generators"].get(), "XLCoverageReport"), V["psm"].get(), *(V["ms"].get().split(";")),
        "--out", V["out"].get(),
        "--error", V["error2"].get(),
        "--ion", ",".join([t for t in ion_types if V[f"ion_{t}"].get()]),
        "--linker", V["linker"].get(),
        "--cfg", V["cfg_pl"].get(),
    )

def run():
    if V["report"].get() == BAReport: run_basicaquisitionreport()
    elif V["report"].get() == TSReport: run_targetselectionreport()
    elif V["report"].get() == TAReport: run_targetaquisitionreport()
    elif V["report"].get() == IDReport: run_identificationreport()
    elif V["report"].get() == XPCReport: run_crosslinkedpeptidecoveragereport()
    else: print("Unknown Report Type:", V["report"].get())

def select_report(_=None, verbose=True):
    for r in reports:
        F[r].grid_forget()
    if verbose: print("selected report type:", V["report"].get())
    F[V["report"].get()].grid(column=0, row=I_fs, columnspan=3, sticky="EW")

util.init_form(main)
I = 0
c = ttk.Combobox(main, textvariable=V["report"], values=reports, state="readonly", justify="center")
c.bind("<<ComboboxSelected>>", select_report)
c.grid(column=0, row=I, columnspan=3, sticky="EW", padx=16, pady=4)
I += 1
F = {r: util.init_form(ttk.Frame(main)) for r in reports}
I_fs = I
I += 1
util.add_entry(main, I, "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
I += 1
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I, columnspan=3)
I += 1
ttk.Separator(main, orient=tk.HORIZONTAL).grid(column=0, row=I, columnspan=3, sticky="EW")
ttk.Label(main, text="Advanced Configuration").grid(column=0, row=I, columnspan=3)
I += 1
util.add_entry(main, I, "Generators:", V["generators"], "Select", util.askdir(V["generators"]))
I += 1

t_ms = (("MES file", "*.mes"), ("MS2 file", "*.ms2"), ("All", "*.*"))
t_target = (("Target List", "*.csv"), ("All", "*.*"))
t_psm = (("PSM", "*.psm"), ("All", "*.*"))

f = F[BAReport]
I = 0
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=t_ms))
I += 1

f = F[TSReport]
I = 0
util.add_entry(f, I, "Target List:", V["target"], "Select", util.askfiles(V["target"], V["out"], filetypes=t_target))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["target_fmt"], values=list(target_fmts.keys()), state="readonly", justify="center"))
I += 1

f = F[TAReport]
I = 0
ttk.Label(f, text=f"{TAReport} Not Available").grid(column=0, row=I, columnspan=3)
I += 1
util.add_entry(f, I, "Target List:", V["target"], "Select", util.askfiles(V["target"], V["out"], filetypes=t_target))
I += 1
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=t_ms))
I += 1

f = F[IDReport]
I = 0
ttk.Label(f, text=f"{IDReport} Not Available").grid(column=0, row=I, columnspan=3)
I += 1
t = (("PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1

f = F[XPCReport]
I = 0
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=t_ms))
I += 1
t = (("PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
_, f_ion, _ = util.add_entry(f, I, "Ion Type:", ttk.Frame(f, height=24))
f_ion1 = ttk.Frame(f_ion)
f_ion2 = ttk.Frame(f_ion)
f_ion1.pack(fill="x")
f_ion2.pack(fill="x")
for n, t in zip(ion_names[0:6], ion_types[0:6]): ttk.Checkbutton(f_ion1, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
for n, t in zip(ion_names[6:], ion_types[6:]): ttk.Checkbutton(f_ion2, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
I += 1
util.add_entry(f, I, "Fragment Mass Error:", V["error2"], "ppm")
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "pLink Cfg. Directory:", V["cfg_pl"], "Select", util.askdir(V["cfg_pl"]))
I += 1

select_report(None, False)
