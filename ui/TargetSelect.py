import os
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

FormSelect = "Regular Peptide"
FormSelectX = "Cross-Linked Peptide"
forms = [FormSelect, FormSelectX]
tds = ["T", "D"]
tdsx = ["TT", "TD", "DD"]
pts = ["Inter", "Intra"]
fmts = ["TW", "TmQE", "TmFu"]
fmts_name = ["TargetWizard", "Thermo Q Exactive", "Thermo Fusion"]
vars_spec = {
    "form": {"type": tk.StringVar, "value": FormSelect},
    "data": {"type": tk.StringVar, "value": ""},
    "psm": {"type": tk.StringVar, "value": ""},
    "out": {"type": tk.StringVar, "value": ""},
    "name": {"type": tk.StringVar, "value": "TargetWizard"},
    "error": {"type": tk.StringVar, "value": "20"},
    "fdr_min": {"type": tk.StringVar, "value": "-Inf"},
    "fdr_max": {"type": tk.StringVar, "value": "Inf"},
    "fdr_ge": {"type": tk.StringVar, "value": "≤"},
    "fdr_le": {"type": tk.StringVar, "value": "≤"},
    "batch": {"type": tk.StringVar, "value": "Inf"},
    "rtime": {"type": tk.StringVar, "value": "240"},
    "lc": {"type": tk.StringVar, "value": "Inf"},
}
for t in tds: vars_spec[f"td_{t}"] = {"type": tk.IntVar, "value": 1}
for t in tdsx: vars_spec[f"td_{t}"] = {"type": tk.IntVar, "value": 1}
for t in pts: vars_spec[f"pt_{t}"] = {"type": tk.IntVar, "value": 1}
for t in fmts: vars_spec[f"fmt_{t}"] = {"type": tk.IntVar, "value": t == "TW"}
task = util.Task("TargetSelect", vars_spec, path=meta.homedir, shared_vars_spec=meta.vars_spec, shared_vars=meta.vars)
V = task.vars

def run_select():
    task.call(V["targetselect"].get(), *(V["data"].get().split(";")),
        "--psm", V["psm"].get(),
        "--out", V["out"].get(),
        "--name", V["name"].get(),
        "--error", V["error"].get(),
        "--fdr_min", V["fdr_min"].get(),
        "--fdr_max", V["fdr_max"].get(),
        *(["--fdr_ge"] if V["fdr_ge"].get() == "≤" else []),
        *(["--fdr_le"] if V["fdr_le"].get() == "≤" else []),
        "--td", ",".join([t for t in tds if V[f"td_{t}"].get()]),
        "--batch", V["batch"].get(),
        "--rtime", V["rtime"].get(),
        "--lc", V["lc"].get(),
        "--fmt", ",".join([t for t in fmts if V[f"fmt_{t}"].get()]),
    )

def run_selectx():
    task.call(V["targetselectx"].get(), *(V["data"].get().split(";")),
        "--psm", V["psm"].get(),
        "--out", V["out"].get(),
        "--name", V["name"].get(),
        "--error", V["error"].get(),
        "--fdr_min", V["fdr_min"].get(),
        "--fdr_max", V["fdr_max"].get(),
        *(["--fdr_ge"] if V["fdr_ge"].get() == "≤" else []),
        *(["--fdr_le"] if V["fdr_le"].get() == "≤" else []),
        "--td", ",".join([t for t in tdsx if V[f"td_{t}"].get()]),
        "--pt", ",".join([t for t in pts if V[f"pt_{t}"].get()]),
        "--batch", V["batch"].get(),
        "--rtime", V["rtime"].get(),
        "--lc", V["lc"].get(),
        "--fmt", ",".join([t for t in fmts if V[f"fmt_{t}"].get()]),
    )

def run():
    if V["form"].get() == FormSelect: run_select()
    elif V["form"].get() == FormSelectX: run_selectx()
    else: print("Unknown Type:", V["form"].get())

def select_form(_=None, verbose=True):
    for v in forms:
        F[v].grid_forget()
    if verbose: print("selected type:", V["form"].get())
    F[V["form"].get()].grid(column=0, row=I_fs, columnspan=3, sticky="EW")

util.init_form(main)
I = 0
c = ttk.Combobox(main, textvariable=V["form"], values=forms, state="readonly", justify="center")
c.bind("<<ComboboxSelected>>", select_form)
c.grid(column=0, row=I, columnspan=3, sticky="EW", padx=16, pady=4)
I += 1
F = {v: util.init_form(ttk.Frame(main)) for v in forms}
I_fs = I
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

f = F[FormSelect]
I = 0
t = (("UMS file", "*.ums"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "Data:", V["data"], "Select", util.askfiles(V["data"], V["out"], filetypes=t))
I += 1
t = (("PSM", "*.spectra"), ("All", "*.*"))
util.add_entry(f, I, "PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
util.add_entry(f, I, "Task Name:", V["name"])
I += 1
util.add_entry(f, I, "Max. MS1 Mass Error:", V["error"], "ppm")
I += 1
_, tmp, _ = util.add_entry(f, I, "FDR Range:", ttk.Frame(f), "%")
ttk.Entry(tmp, textvariable=V["fdr_min"]).pack(side="left", fill="x", expand=True)
ttk.Label(tmp, text="% ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_ge"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Label(tmp, text=" FDR ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_le"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Entry(tmp, textvariable=V["fdr_max"]).pack(side="left", fill="x", expand=True)
I += 1
_, tmp, _ = util.add_entry(f, I, "Target / Decoy Type:", ttk.Frame(f))
for t in tds: ttk.Checkbutton(tmp, text=t, variable=V[f"td_{t}"]).pack(side="left", expand=True)

f = F[FormSelectX]
I = 0
t = (("UMS file", "*.ums"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "Data:", V["data"], "Select", util.askfiles(V["data"], V["out"], filetypes=t))
I += 1
t = (("PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
util.add_entry(f, I, "Task Name:", V["name"])
I += 1
util.add_entry(f, I, "Max. MS1 Mass Error:", V["error"], "ppm")
I += 1
_, tmp, _ = util.add_entry(f, I, "FDR Range:", ttk.Frame(f), "%")
ttk.Entry(tmp, textvariable=V["fdr_min"]).pack(side="left", fill="x", expand=True)
ttk.Label(tmp, text="% ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_ge"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Label(tmp, text=" FDR ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_le"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Entry(tmp, textvariable=V["fdr_max"]).pack(side="left", fill="x", expand=True)
I += 1
_, tmp, _ = util.add_entry(f, I, "Target / Decoy Type:", ttk.Frame(f))
for t in tdsx: ttk.Checkbutton(tmp, text=t, variable=V[f"td_{t}"]).pack(side="left", expand=True)
I += 1
_, tmp, _ = util.add_entry(f, I, "Inter- / Intra-Protein:", ttk.Frame(f))
for t in pts: ttk.Checkbutton(tmp, text=t, variable=V[f"pt_{t}"]).pack(side="left", expand=True)

select_form(None, False)
