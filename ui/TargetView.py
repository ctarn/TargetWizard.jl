import os
import random
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

RView = "Regular Linear Peptide View"
CView = "Comparative Linear Peptide View"
RXView = "Regular Cross-linked Peptide View"
CXView = "Comparative Cross-linked Peptide View"
views = [RView, CView, RXView, CXView]
target_fmts = {"Auto Detect": "auto", "TargetWizard": "TW", "Thermo Q Exactive": "TmQE", "Thermo Fusion": "TmFu"}
vars_spec = {
    "view": {"type": tk.StringVar, "value": RXView},
    "tg": {"type": tk.StringVar, "value": ""},
    "target_fmt": {"type": tk.StringVar, "value": "Auto Detect"},
    "xl": {"type": tk.StringVar, "value": ""},
    "ms": {"type": tk.StringVar, "value": ""},
    "ms_": {"type": tk.StringVar, "value": ""},
    "ft": {"type": tk.StringVar, "value": ""},
    "ft_": {"type": tk.StringVar, "value": ""},
    "psm": {"type": tk.StringVar, "value": ""},
    "psm_xl": {"type": tk.StringVar, "value": ""},
    "psm_xl_": {"type": tk.StringVar, "value": ""},
    "fdr": {"type": tk.StringVar, "value": "Inf"},
    "error": {"type": tk.StringVar, "value": "20.0"},
    "linker": {"type": tk.StringVar, "value": "DSSO"},
    "out": {"type": tk.StringVar, "value": ""},
    "targetview": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
    "cfg_pl": {"type": tk.StringVar, "value": ""},
    "cfg_pf": {"type": tk.StringVar, "value": ""},
    "url": {"type": tk.StringVar, "value": "127.0.0.1:30030"},
}
task = util.Task("TargetView", vars_spec, path=meta.homedir)
V = task.vars

def run_view():
    pass

def run_dualview():
    pass

def run_xview():
    task.call(os.path.join(V["targetview"].get(), "TargetXView"),
        V["tg"].get(), "--out", V["out"].get(),
        "--fmt", target_fmts[V["target_fmt"].get()],
        "--error", V["error"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--cfg_pf", V["cfg_pf"].get(),
        "--linker", V["linker"].get(),
        "--xl", V["xl"].get(),
        "--ms", V["ms"].get(),
        "--ft", V["ft"].get(),
        "--psm", V["psm_xl"].get(),
        "--psm_pf", V["psm"].get(),
        "--fdr", V["fdr"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_xdualview():
    task.call(os.path.join(V["targetview"].get(), "TargetXDualView"),
        V["tg"].get(), "--out", V["out"].get(),
        "--fmt", target_fmts[V["target_fmt"].get()],
        "--error", V["error"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--linker", V["linker"].get(),
        "--xl", V["xl"].get(),
        "--ms", V["ms"].get(), V["ms_"].get(),
        "--ft", V["ft"].get(), V["ft_"].get(),
        "--psm", V["psm_xl"].get(), V["psm_xl_"].get(),
        "--fdr", V["fdr"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run():
    if V["view"].get() == RView: run_view()
    elif V["view"].get() == CView: run_dualview()
    elif V["view"].get() == RXView: run_xview()
    elif V["view"].get() == CXView: run_xdualview()
    else: print("Unknown View Type:", V["view"].get())

def new_port():
    p = str(random.randint(49152, 65535))
    host = V["url"].get().split(":")[0]
    V["url"].set(host + ":" + p)

new_port()

def select_view(_=None, verbose=True):
    for v in views:
        F[v].grid_forget()
    if verbose: print("selected view type:", V["view"].get())
    F[V["view"].get()].grid(column=0, row=I_fs, columnspan=3, sticky="EW")

util.init_form(main)
I = 0
c = ttk.Combobox(main, textvariable=V["view"], values=views, state="readonly", justify="center")
c.bind("<<ComboboxSelected>>", select_view)
c.grid(column=0, row=I, columnspan=3, sticky="EW", padx=16, pady=4)
I += 1
F = {v: util.init_form(ttk.Frame(main)) for v in views}
I_fs = I
I += 1
util.add_entry(main, I, "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
I += 1
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I, columnspan=3)
I += 1
ttk.Separator(main, orient=tk.HORIZONTAL).grid(column=0, row=I, columnspan=3, sticky="EW")
ttk.Label(main, text="Advanced Configuration").grid(column=0, row=I, columnspan=3)
I += 1
util.add_entry(main, I, "TargetView:", V["targetview"], "Select", util.askfile(V["targetview"]))
I += 1
util.add_entry(main, I, "URL:", V["url"], "New Port", new_port)
I += 1


t = (("Target List", "*.csv"), ("All", "*.*"))

f = F[RView]
I = 0
ttk.Label(f, text=f"{RView} Not Available").grid(column=0, row=I, columnspan=3)
I += 1

f = F[CView]
I = 0
ttk.Label(f, text=f"{CView} Not Available").grid(column=0, row=I, columnspan=3)
I += 1

f = F[RXView]
I = 0
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfiles(V["tg"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["target_fmt"], values=list(target_fmts.keys()), state="readonly", justify="center"))
I += 1
t = (("Candidate XL List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Candidate XL List:", V["xl"], "Select", util.askfile(V["xl"], filetypes=t))
I += 1
t = (("MES file", "*.mes"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=t))
I += 1
t = (("Peptide Feature List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=t))
I += 1
t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "XL PSM:", V["psm_xl"], "Select", util.askfile(V["psm_xl"], filetypes=t))
I += 1
t = (("Linear PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Linear PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
util.add_entry(f, I, "FDR:", V["fdr"], "%")
I += 1
util.add_entry(f, I, "MS1 Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "pLink Cfg. Directory:", V["cfg_pl"], "Select", util.askdir(V["cfg_pl"]))
I += 1
util.add_entry(f, I, "pFind Cfg. Directory:", V["cfg_pf"], "Select", util.askdir(V["cfg_pf"]))
I += 1

f = F[CXView]
I = 0
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfiles(V["tg"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["target_fmt"], values=list(target_fmts.keys()), state="readonly", justify="center"))
I += 1
t = (("Candidate XL List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Candidate XL List:", V["xl"], "Select", util.askfile(V["xl"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Data A").grid(column=0, row=I)
I += 1
t = (("MES file", "*.mes"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=t))
I += 1
t = (("Peptide Feature List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=t))
I += 1
t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "XL PSM:", V["psm_xl"], "Select", util.askfile(V["psm_xl"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Data B").grid(column=0, row=I)
I += 1
t = (("MES file", "*.mes"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "MS Data:", V["ms_"], "Select", util.askfile(V["ms_"], filetypes=t))
I += 1
t = (("Peptide Feature List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Feature List:", V["ft_"], "Select", util.askfile(V["ft_"], filetypes=t))
I += 1
t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "XL PSM:", V["psm_xl_"], "Select", util.askfile(V["psm_xl_"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
I += 1
util.add_entry(f, I, "FDR:", V["fdr"], "%")
I += 1
util.add_entry(f, I, "MS1 Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "pLink Cfg. Directory:", V["cfg_pl"], "Select", util.askdir(V["cfg_pl"]))
I += 1

select_view(None, False)
