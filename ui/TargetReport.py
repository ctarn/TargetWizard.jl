import os
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

FlowBA = "Basic Acquisition Report"
FlowTS = "Target Selection Report"
FlowTA = "Target Acquisition Report"
FlowTAXL = "Target Acquisition Cross-linked Peptide Report"
FlowPC = "Peptide Coverage Report"
FlowPCXL = "Cross-linked Peptide Coverage Report"
FlowSNRDXL = "Comparative Cross-link Signal-to-Noise Ratio Report"
flows = [FlowBA, FlowTS, FlowTA, FlowTAXL, FlowPC, FlowPCXL, FlowSNRDXL]
ion_types = ["a", "b", "c", "x", "y", "z", "a_NH3", "b_NH3", "y_NH3", "a_H2O", "b_H2O", "y_H2O"]
ion_names = ["a", "b", "c", "x", "y", "z", "a-NH₃", "b-NH₃", "y-NH₃", "a-H₂O", "b-H₂O", "y-H₂O"]
vars_spec = {
    "flow": {"type": tk.StringVar, "value": FlowBA},
    "out": {"type": tk.StringVar, "value": ""},
    "ms": {"type": tk.StringVar, "value": ""},
    "ms_": {"type": tk.StringVar, "value": ""},
    "tg": {"type": tk.StringVar, "value": ""},
    "ft": {"type": tk.StringVar, "value": ""},
    "xl": {"type": tk.StringVar, "value": ""},
    "fmt_tg": {"type": tk.StringVar, "value": "Auto Detect"},
    "psm": {"type": tk.StringVar, "value": ""},
    "psm_": {"type": tk.StringVar, "value": ""},
    "linker": {"type": tk.StringVar, "value": "DSSO"},
    "error": {"type": tk.StringVar, "value": "20.0"},
    "fdr": {"type": tk.StringVar, "value": "5.0"},
    "ms_sim_thres": {"type": tk.StringVar, "value": "0.5"},
}
for t in ion_types: vars_spec[f"ion_{t}"] = {"type": tk.IntVar, "value": 1}
task = util.Task("TargetReport", vars_spec, path=meta.homedir, shared_vars_spec=meta.vars_spec, shared_vars=meta.vars)
V = task.vars

def run_basic_acquisition_report():
    task.call(os.path.join(V["generators"].get(), "BasicAcquisitionReport"), *(V["ms"].get().split(";")), "--out", V["out"].get())

def run_target_selection_report():
    task.call(os.path.join(V["generators"].get(), "TargetSelectionReport"),
        *(V["tg"].get().split(";")), "--out", V["out"].get(),
        "--fmt", meta.fmts_tg[V["fmt_tg"].get()],
    )

def run_target_acquisition_report():
    task.call(os.path.join(V["generators"].get(), "TargetAcquisitionReport"),
        V["tg"].get(),
        "--ms", V["ms"].get(),
        "--ms_old", *(V["ms_"].get().split(";")),
        "--psm", V["psm"].get(),
        "--out", V["out"].get(),
        "--ft", V["ft"].get(),
        "--fmt", meta.fmts_tg[V["fmt_tg"].get()],
        "--error", V["error"].get(),
        "--ms_sim_thres", V["ms_sim_thres"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pf"].get(),
    )

def run_target_acquisition_xl_report():
    task.call(os.path.join(V["generators"].get(), "TargetAcquisitionXLReport"),
        V["tg"].get(),
        "--ms", V["ms"].get(),
        "--ms_old", *(V["ms_"].get().split(";")),
        "--psm", V["psm"].get(),
        "--out", V["out"].get(),
        "--xl", V["xl"].get(),
        "--ft", V["ft"].get(),
        "--psm_pf", V["psm_"].get(),
        "--fmt", meta.fmts_tg[V["fmt_tg"].get()],
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--ms_sim_thres", V["ms_sim_thres"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--cfg_pf", V["cfg_pf"].get(),
    )

def run_peptide_coverage_report():
    task.call(os.path.join(V["generators"].get(), "PeptideCoverageReport"), V["psm"].get(),
        "--ms", *(V["ms"].get().split(";")),
        "--out", V["out"].get(),
        "--error", V["error"].get(),
        "--ion", ",".join([t for t in ion_types if V[f"ion_{t}"].get()]),
        "--cfg", V["cfg_pf"].get(),
    )

def run_peptide_coverage_xl_report():
    task.call(os.path.join(V["generators"].get(), "PeptideCoverageXLReport"), V["psm"].get(),
        "--ms", *(V["ms"].get().split(";")),
        "--out", V["out"].get(),
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--ion", ",".join([t for t in ion_types if V[f"ion_{t}"].get()]),
        "--cfg", V["cfg_pl"].get(),
    )

def run_noise_ratio_dual_xl():
    task.call(os.path.join(V["generators"].get(), "NoiseRatioDualXLReport"), V["tg"].get(),
        "--ms", V["ms"].get(), V["ms_"].get(),
        "--psm", V["psm"].get(), V["psm_"].get(),
        "--out", V["out"].get(),
        "--fmt", meta.fmts_tg[V["fmt_tg"].get()],
        "--linker", V["linker"].get(),
        "--fdr", V["fdr"].get(),
        "--ion", ",".join([t for t in ion_types if V[f"ion_{t}"].get()]),
        "--error", V["error"].get(),
        "--cfg", V["cfg_pl"].get(),
    )

def run():
    if V["flow"].get() == FlowBA: run_basic_acquisition_report()
    elif V["flow"].get() == FlowTS: run_target_selection_report()
    elif V["flow"].get() == FlowTA: run_target_acquisition_report()
    elif V["flow"].get() == FlowTAXL: run_target_acquisition_xl_report()
    elif V["flow"].get() == FlowPC: run_peptide_coverage_report()
    elif V["flow"].get() == FlowPCXL: run_peptide_coverage_xl_report()
    elif V["flow"].get() == FlowSNRDXL: run_noise_ratio_dual_xl()
    else: print("Unknown Flow Type:", V["flow"].get())

def select_flow(_=None, verbose=True):
    for f in flows: F[f].grid_forget()
    if verbose: print("selected flow type:", V["flow"].get())
    F[V["flow"].get()].grid(column=0, row=I_fs, columnspan=3, sticky="EW")

util.init_form(main)
I = util.Counter()
c = ttk.Combobox(main, textvariable=V["flow"], values=flows, state="readonly", justify="center")
c.bind("<<ComboboxSelected>>", select_flow)
c.grid(column=0, row=I.next(), columnspan=3, sticky="EW", padx=16, pady=4)
F = {f: util.init_form(ttk.Frame(main)) for f in flows}
I_fs = I.next()
util.add_entry(main, I.next(), "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I.next(), columnspan=3)

f = F[FlowBA]
I = util.Counter()
util.add_entry(f, I.next(), "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=meta.filetype_ms))

f = F[FlowTS]
I = util.Counter()
util.add_entry(f, I.next(), "Target List:", V["tg"], "Select", util.askfiles(V["tg"], V["out"], filetypes=meta.filetype_tg))
util.add_combobox(f, I.next(), "List Format:", V["fmt_tg"], list(meta.fmts_tg.keys()))

f = F[FlowTA]
I = util.Counter()
util.add_entry(f, I.next(), "Target List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=meta.filetype_tg))
util.add_combobox(f, I.next(), "List Format:", V["fmt_tg"], list(meta.fmts_tg.keys()))
util.add_entry(f, I.next(), "Targeted MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "Targeted MS PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm))
util.add_entry(f, I.next(), "Original MS Data:", V["ms_"], "Select", util.askfiles(V["ms_"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "Max. MS1 Mass Error:", V["error"], "ppm")
util.add_entry(f, I.next(), "FDR Threshold:", V["fdr"], "%")
util.add_entry(f, I.next(), "MS Sim. Thres.:", V["ms_sim_thres"])
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I.next(), sticky="EW", padx=12)
ttk.Label(f, text="Optional").grid(column=0, row=I.next())
util.add_entry(f, I.next(), "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=meta.filetype_prec))
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I.next(), sticky="EW", padx=12)

f = F[FlowTAXL]
I = util.Counter()
util.add_entry(f, I.next(), "Target List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=meta.filetype_tg))
util.add_combobox(f, I.next(), "List Format:", V["fmt_tg"], list(meta.fmts_tg.keys()))
util.add_entry(f, I.next(), "Targeted MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "Targeted MS PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm_xl))
util.add_entry(f, I.next(), "Original MS Data:", V["ms_"], "Select", util.askfiles(V["ms_"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "Default Linker:", V["linker"])
util.add_entry(f, I.next(), "Max. MS1 Mass Error:", V["error"], "ppm")
util.add_entry(f, I.next(), "FDR Threshold:", V["fdr"], "%")
util.add_entry(f, I.next(), "MS Sim. Thres.:", V["ms_sim_thres"])
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I.next(), sticky="EW", padx=12)
ttk.Label(f, text="Optional").grid(column=0, row=I.next())
util.add_entry(f, I.next(), "Candidate XL List:", V["xl"], "Select", util.askfile(V["xl"], filetypes=meta.filetype_xl))
util.add_entry(f, I.next(), "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=meta.filetype_prec))
util.add_entry(f, I.next(), "Linear PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm))
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I.next(), sticky="EW", padx=12)

f = F[FlowPC]
I = util.Counter()
util.add_entry(f, I.next(), "PSM:", V["psm"], "Select", util.askfile(V["psm"], V["out"], filetypes=meta.filetype_psm))
util.add_entry(f, I.next(), "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=meta.filetype_ms))
_, f_ion, _ = util.add_entry(f, I.next(), "Ion Type:", ttk.Frame(f, height=24))
f_ion1 = ttk.Frame(f_ion)
f_ion2 = ttk.Frame(f_ion)
f_ion1.pack(fill="x")
f_ion2.pack(fill="x")
for n, t in zip(ion_names[0:6], ion_types[0:6]): ttk.Checkbutton(f_ion1, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
for n, t in zip(ion_names[6:], ion_types[6:]): ttk.Checkbutton(f_ion2, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
util.add_entry(f, I.next(), "Fragment Mass Error:", V["error"], "ppm")

f = F[FlowPCXL]
I = util.Counter()
util.add_entry(f, I.next(), "PSM:", V["psm"], "Select", util.askfile(V["psm"], V["out"], filetypes=meta.filetype_psm_xl))
util.add_entry(f, I.next(), "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=meta.filetype_ms))
_, f_ion, _ = util.add_entry(f, I.next(), "Ion Type:", ttk.Frame(f, height=24))
f_ion1 = ttk.Frame(f_ion)
f_ion2 = ttk.Frame(f_ion)
f_ion1.pack(fill="x")
f_ion2.pack(fill="x")
for n, t in zip(ion_names[0:6], ion_types[0:6]): ttk.Checkbutton(f_ion1, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
for n, t in zip(ion_names[6:], ion_types[6:]): ttk.Checkbutton(f_ion2, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
util.add_entry(f, I.next(), "Default Linker:", V["linker"])
util.add_entry(f, I.next(), "Fragment Mass Error:", V["error"], "ppm")

f = F[FlowSNRDXL]
I = util.Counter()
util.add_entry(f, I.next(), "Target List:", V["tg"], "Select", util.askfiles(V["tg"], V["out"], filetypes=meta.filetype_tg))
util.add_combobox(f, I.next(), "List Format:", V["fmt_tg"], list(meta.fmts_tg.keys()))
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I.next(), sticky="EW", padx=12)
ttk.Label(f, text="Data A").grid(column=0, row=I.next())
util.add_entry(f, I.next(), "MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "XL PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm_xl))
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I.next(), sticky="EW", padx=12)
ttk.Label(f, text="Data B").grid(column=0, row=I.next())
util.add_entry(f, I.next(), "MS Data:", V["ms_"], "Select", util.askfile(V["ms_"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "XL PSM:", V["psm_"], "Select", util.askfile(V["psm_"], filetypes=meta.filetype_psm_xl))
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I.next(), sticky="EW", padx=12)
util.add_entry(f, I.next(), "Default Linker:", V["linker"])
util.add_entry(f, I.next(), "FDR:", V["fdr"], "%")
_, f_ion, _ = util.add_entry(f, I.next(), "Ion Type:", ttk.Frame(f, height=24))
f_ion1 = ttk.Frame(f_ion)
f_ion2 = ttk.Frame(f_ion)
f_ion1.pack(fill="x")
f_ion2.pack(fill="x")
for n, t in zip(ion_names[0:6], ion_types[0:6]): ttk.Checkbutton(f_ion1, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
for n, t in zip(ion_names[6:], ion_types[6:]): ttk.Checkbutton(f_ion2, text=n, variable=V[f"ion_{t}"]).pack(side="left", expand=True)
util.add_entry(f, I.next(), "Fragment Mass Error:", V["error"], "ppm")

select_flow(None, False)
