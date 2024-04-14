import os
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

BAReport = "Basic Aquisition Report"
TSReport = "Target Selection Report"
TAReport = "Target Aquisition Report"
IDReport = "Identification Report"
PCReport = "Peptide Coverage Report"
XPCReport = "Cross-linked Peptide Coverage Report"
SNRReportDX = "Comparative Cross-link Signal-to-Noise Ratio Report"
reports = [BAReport, TSReport, TAReport, IDReport, PCReport, XPCReport, SNRReportDX]
target_fmts = {"Auto Detect": "auto", "TargetWizard": "TW", "Thermo Q Exactive": "TmQE", "Thermo Fusion": "TmFu"}
ion_types = ["a", "b", "c", "x", "y", "z", "a_NH3", "b_NH3", "y_NH3", "a_H2O", "b_H2O", "y_H2O"]
ion_names = ["a", "b", "c", "x", "y", "z", "a-NH₃", "b-NH₃", "y-NH₃", "a-H₂O", "b-H₂O", "y-H₂O"]
vars_spec = {
    "report": {"type": tk.StringVar, "value": BAReport},
    "out": {"type": tk.StringVar, "value": ""},
    "ms": {"type": tk.StringVar, "value": ""},
    "ms_": {"type": tk.StringVar, "value": ""},
    "target": {"type": tk.StringVar, "value": ""},
    "target_fmt": {"type": tk.StringVar, "value": "Auto Detect"},
    "psm": {"type": tk.StringVar, "value": ""},
    "psm_": {"type": tk.StringVar, "value": ""},
    "linker": {"type": tk.StringVar, "value": "DSSO"},
    "error2": {"type": tk.StringVar, "value": "20.0"},
    "fdr": {"type": tk.StringVar, "value": "5.0"},
}
for t in ion_types: vars_spec[f"ion_{t}"] = {"type": tk.IntVar, "value": 1}
task = util.Task("TargetReport", vars_spec, path=meta.homedir, shared_vars_spec=meta.vars_spec, shared_vars=meta.vars)
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

def run_peptidecoveragereport():
    task.call(os.path.join(V["generators"].get(), "CoverageReport"), V["psm"].get(),
        "--ms", *(V["ms"].get().split(";")),
        "--out", V["out"].get(),
        "--error", V["error2"].get(),
        "--ion", ",".join([t for t in ion_types if V[f"ion_{t}"].get()]),
        "--cfg", V["cfg_pl"].get(),
    )

def run_crosslinkedpeptidecoveragereport():
    task.call(os.path.join(V["generators"].get(), "CoverageReportX"), V["psm"].get(),
        "--ms", *(V["ms"].get().split(";")),
        "--out", V["out"].get(),
        "--linker", V["linker"].get(),
        "--error", V["error2"].get(),
        "--ion", ",".join([t for t in ion_types if V[f"ion_{t}"].get()]),
        "--cfg", V["cfg_pl"].get(),
    )

def run_snr_report_dx():
    task.call(os.path.join(V["generators"].get(), "NoiseRatioReportDualX"), V["target"].get(),
        "--ms", V["ms"].get(), V["ms_"].get(),
        "--psm", V["psm"].get(), V["psm_"].get(),
        "--out", V["out"].get(),
        "--fmt", target_fmts[V["target_fmt"].get()],
        "--linker", V["linker"].get(),
        "--fdr", V["fdr"].get(),
        "--ion", ",".join([t for t in ion_types if V[f"ion_{t}"].get()]),
        "--error", V["error2"].get(),
        "--cfg", V["cfg_pl"].get(),
    )

def run():
    if V["report"].get() == BAReport: run_basicaquisitionreport()
    elif V["report"].get() == TSReport: run_targetselectionreport()
    elif V["report"].get() == TAReport: run_targetaquisitionreport()
    elif V["report"].get() == IDReport: run_identificationreport()
    elif V["report"].get() == PCReport: run_peptidecoveragereport()
    elif V["report"].get() == XPCReport: run_crosslinkedpeptidecoveragereport()
    elif V["report"].get() == SNRReportDX: run_snr_report_dx()
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

t_ms = (("UMS file", "*.ums"), ("MS2 file", "*.ms2"), ("All", "*.*"))
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

f = F[PCReport]
I = 0
t = (("PSM", "*.csv"), ("PSM", "*.spectra"), ("All", "*.*"))
util.add_entry(f, I, "PSM:", V["psm"], "Select", util.askfile(V["psm"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=t_ms))
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

f = F[XPCReport]
I = 0
t = (("PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "PSM:", V["psm"], "Select", util.askfile(V["psm"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=t_ms))
I += 1
_, f_ion, _ = util.add_entry(f, I, "Ion Type:", ttk.Frame(f, height=24))
f_ion1 = ttk.Frame(f_ion)
f_ion2 = ttk.Frame(f_ion)
f_ion1.pack(fill="x")
f_ion2.pack(fill="x")
for n, t in zip(ion_names[0:6], ion_types[0:6]): ttk.Checkbutton(f_ion1, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
for n, t in zip(ion_names[6:], ion_types[6:]): ttk.Checkbutton(f_ion2, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "Fragment Mass Error:", V["error2"], "ppm")

f = F[SNRReportDX]
I = 0
t = (("Target List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Target List:", V["target"], "Select", util.askfiles(V["target"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["target_fmt"], values=list(target_fmts.keys()), state="readonly", justify="center"))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Data A").grid(column=0, row=I)
I += 1
t = (("UMS file", "*.ums"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=t))
I += 1
t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "XL PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Data B").grid(column=0, row=I)
I += 1
t = (("UMS file", "*.ums"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "MS Data:", V["ms_"], "Select", util.askfile(V["ms_"], filetypes=t))
I += 1
t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "XL PSM:", V["psm_"], "Select", util.askfile(V["psm_"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "FDR:", V["fdr"], "%")
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

select_report(None, False)
