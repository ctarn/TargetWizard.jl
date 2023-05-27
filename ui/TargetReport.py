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
reports = [BAReport, TSReport, TAReport, IDReport]
target_fmts = {"Auto Detect": "auto", "TargetWizard": "TW", "Thermo Q Exactive": "TmQE", "Thermo Fusion": "TmFu"}
vars_spec = {
    "report": {"type": tk.StringVar, "value": BAReport},
    "ms": {"type": tk.StringVar, "value": ""},
    "target": {"type": tk.StringVar, "value": ""},
    "target_fmt": {"type": tk.StringVar, "value": "Auto Detect"},
    "psm": {"type": tk.StringVar, "value": ""},
    "out": {"type": tk.StringVar, "value": ""},
    "generators": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
    "thermorawread": {"type": tk.StringVar, "value": util.get_content("ThermoRawRead", "ThermoRawRead.exe", shared=True)},
    "mono": {"type": tk.StringVar, "value": path_mono},
}
task = util.Task("TargetReport", vars_spec, path=meta.homedir)
V = task.vars

def run_basicaquisitionreport():
    paths = []
    for p in V["ms"].get().split(";"):
        ext = os.path.splitext(p)[1].lower()
        if ext == ".raw": p = run_thermorawread(p, V["out"].get())
        paths.append(p)
    task.call(os.path.join(V["generators"].get(), "BasicAquisitionReport"), *paths, "--out", V["out"].get())

def run_targetselectionreport():
    task.call(os.path.join(V["generators"].get(), "TargetSelectionReport"),
        *(V["target"].get().split(";")), "--out", V["out"].get(),
        "--fmt", target_fmts[V["target_fmt"].get()],
    )

def run_targetaquisitionreport():
    pass

def run_identificationreport():
    pass

def run_thermorawread(data, out):
    task.call(*([] if util.is_windows else [V["mono"].get()]), V["thermorawread"].get(), data, out)
    return os.path.join(out, os.path.splitext(os.path.basename(data))[0] + ".ms2")

def run():
    if V["report"].get() == BAReport: run_basicaquisitionreport()
    elif V["report"].get() == TSReport: run_targetselectionreport()
    elif V["report"].get() == TAReport: run_targetaquisitionreport()
    elif V["report"].get() == IDReport: run_identificationreport()
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
util.add_entry(main, I, "ThermoRawRead:", V["thermorawread"], "Select", util.askfile(V["thermorawread"]))
I += 1
if not util.is_windows:
    util.add_entry(main, I, "Mono Runtime:", V["mono"], "Select", util.askfile(V["mono"]))
    I += 1

t_ms = (("MS2", "*.ms2"), ("RAW", "*.raw"), ("All", "*.*"))
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

select_report(None, False)
