import os
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

FlowR = "Regular Peptide"
FlowXL = "Cross-Linked Peptide"
flows = [FlowR, FlowXL]
tds = ["T", "D"]
tds_xl = ["TT", "TD", "DD"]
pts = ["Inter", "Intra"]
fmts = ["TW", "TmQE", "TmFu"]
fmts_name = ["TargetWizard", "Thermo Q Exactive", "Thermo Fusion"]
vars_spec = {
    "flow": {"type": tk.StringVar, "value": FlowR},
    "ms": {"type": tk.StringVar, "value": ""},
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
for t in tds_xl: vars_spec[f"td_{t}"] = {"type": tk.IntVar, "value": 1}
for t in pts: vars_spec[f"pt_{t}"] = {"type": tk.IntVar, "value": 1}
for t in fmts: vars_spec[f"fmt_{t}"] = {"type": tk.IntVar, "value": t == "TW"}
task = util.Task("TargetSelect", vars_spec, path=meta.homedir, shared_vars_spec=meta.vars_spec, shared_vars=meta.vars)
V = task.vars

def run_select():
    task.call(V["targetselect"].get(), *(V["ms"].get().split(";")),
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

def run_selectxl():
    task.call(V["targetselectxl"].get(), *(V["ms"].get().split(";")),
        "--psm", V["psm"].get(),
        "--out", V["out"].get(),
        "--name", V["name"].get(),
        "--error", V["error"].get(),
        "--fdr_min", V["fdr_min"].get(),
        "--fdr_max", V["fdr_max"].get(),
        *(["--fdr_ge"] if V["fdr_ge"].get() == "≤" else []),
        *(["--fdr_le"] if V["fdr_le"].get() == "≤" else []),
        "--td", ",".join([t for t in tds_xl if V[f"td_{t}"].get()]),
        "--pt", ",".join([t for t in pts if V[f"pt_{t}"].get()]),
        "--batch", V["batch"].get(),
        "--rtime", V["rtime"].get(),
        "--lc", V["lc"].get(),
        "--fmt", ",".join([t for t in fmts if V[f"fmt_{t}"].get()]),
    )

def run():
    if V["flow"].get() == FlowR: run_select()
    elif V["flow"].get() == FlowXL: run_selectxl()
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
util.add_entry(main, I.next(), "Batch Size:", V["batch"])
util.add_entry(main, I.next(), "RT Window:", V["rtime"], "sec")
util.add_entry(main, I.next(), "LC Length:", V["lc"], "min")
_, f, _ = util.add_entry(main, I.next(), "List Format:", ttk.Frame(main))
for t, n in zip(fmts, fmts_name): ttk.Checkbutton(f, text=n, variable=V[f"fmt_{t}"]).pack(side="left", expand=True)
util.add_entry(main, I.next(), "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I.next(), columnspan=3)

f = F[FlowR]
I = util.Counter()
util.add_entry(f, I.next(), "Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm))
util.add_entry(f, I.next(), "Task Name:", V["name"])
util.add_entry(f, I.next(), "Max. MS1 Mass Error:", V["error"], "ppm")
_, tmp, _ = util.add_entry(f, I.next(), "FDR Range:", ttk.Frame(f), "%")
ttk.Entry(tmp, textvariable=V["fdr_min"]).pack(side="left", fill="x", expand=True)
ttk.Label(tmp, text="% ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_ge"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Label(tmp, text=" FDR ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_le"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Entry(tmp, textvariable=V["fdr_max"]).pack(side="left", fill="x", expand=True)
_, tmp, _ = util.add_entry(f, I.next(), "Target / Decoy Type:", ttk.Frame(f))
for t in tds: ttk.Checkbutton(tmp, text=t, variable=V[f"td_{t}"]).pack(side="left", expand=True)

f = F[FlowXL]
I = util.Counter()
util.add_entry(f, I.next(), "Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=meta.filetype_ms))
util.add_entry(f, I.next(), "PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm_xl))
util.add_entry(f, I.next(), "Task Name:", V["name"])
util.add_entry(f, I.next(), "Max. MS1 Mass Error:", V["error"], "ppm")
_, tmp, _ = util.add_entry(f, I.next(), "FDR Range:", ttk.Frame(f), "%")
ttk.Entry(tmp, textvariable=V["fdr_min"]).pack(side="left", fill="x", expand=True)
ttk.Label(tmp, text="% ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_ge"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Label(tmp, text=" FDR ").pack(side="left")
ttk.Combobox(tmp, textvariable=V["fdr_le"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Entry(tmp, textvariable=V["fdr_max"]).pack(side="left", fill="x", expand=True)
_, tmp, _ = util.add_entry(f, I.next(), "Target / Decoy Type:", ttk.Frame(f))
for t in tds_xl: ttk.Checkbutton(tmp, text=t, variable=V[f"td_{t}"]).pack(side="left", expand=True)
_, tmp, _ = util.add_entry(f, I.next(), "Inter- / Intra-Protein:", ttk.Frame(f))
for t in pts: ttk.Checkbutton(tmp, text=t, variable=V[f"pt_{t}"]).pack(side="left", expand=True)

select_flow(None, False)
