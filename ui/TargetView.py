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
XSView = "Cross-link Site View"
XESView = "Cross-link Exhaustive Search View"
views = [RView, CView, RXView, CXView, XSView, XESView]
fmt_targets = {"Auto Detect": "auto", "TargetWizard": "TW", "Thermo Q Exactive": "TmQE", "Thermo Fusion": "TmFu"}
vars_spec = {
    "view": {"type": tk.StringVar, "value": RXView},
    "out": {"type": tk.StringVar, "value": ""},
    "tg": {"type": tk.StringVar, "value": ""},
    "fmt_target": {"type": tk.StringVar, "value": "Auto Detect"},
    "ms": {"type": tk.StringVar, "value": ""},
    "ms_": {"type": tk.StringVar, "value": ""},
    "psm": {"type": tk.StringVar, "value": ""},
    "psm_xl": {"type": tk.StringVar, "value": ""},
    "psm_xl_": {"type": tk.StringVar, "value": ""},
    "xl": {"type": tk.StringVar, "value": ""},
    "ft": {"type": tk.StringVar, "value": ""},
    "ft_": {"type": tk.StringVar, "value": ""},
    "linker": {"type": tk.StringVar, "value": "DSSO"},
    "error": {"type": tk.StringVar, "value": "20.0"},
    "fdr": {"type": tk.StringVar, "value": "Inf"},
    "ms_sim_thres": {"type": tk.StringVar, "value": "0.5"},
    "url": {"type": tk.StringVar, "value": "127.0.0.1:30030"},
    "inten_thres": {"type": tk.StringVar, "value": "0.0"},
    "smooth": {"type": tk.StringVar, "value": "16"},
}
task = util.Task("TargetView", vars_spec, path=meta.homedir, shared_vars_spec=meta.vars_spec, shared_vars=meta.vars)
V = task.vars

def run_view():
    task.call(os.path.join(V["viewers"].get(), "TargetView"),
        V["tg"].get(),
        "--ms", V["ms"].get(),
        "--ms_old", *(V["ms_"].get().split(";")),
        "--psm", V["psm"].get(),
        "--out", V["out"].get(),
        "--ft", V["ft"].get(),
        "--fmt", fmt_targets[V["fmt_target"].get()],
        "--error", V["error"].get(),
        "--ms_sim_thres", V["ms_sim_thres"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pf"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_dualview():
    pass

def run_xview():
    task.call(os.path.join(V["viewers"].get(), "TargetViewX"),
        V["tg"].get(),
        "--ms", V["ms"].get(),
        "--ms_old", *(V["ms_"].get().split(";")),
        "--psm", V["psm_xl"].get(),
        "--out", V["out"].get(),
        "--xl", V["xl"].get(),
        "--ft", V["ft"].get(),
        "--psm_pf", V["psm"].get(),
        "--fmt", fmt_targets[V["fmt_target"].get()],
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--ms_sim_thres", V["ms_sim_thres"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--cfg_pf", V["cfg_pf"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_xdualview():
    task.call(os.path.join(V["viewers"].get(), "TargetDualViewX"),
        V["tg"].get(),
        "--ms", V["ms"].get(), V["ms_"].get(),
        "--psm", V["psm_xl"].get(), V["psm_xl_"].get(),
        "--out", V["out"].get(),
        "--xl", V["xl"].get(),
        "--ft", V["ft"].get(), V["ft_"].get(),
        "--fmt", fmt_targets[V["fmt_target"].get()],
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_xsview():
    task.call(os.path.join(V["viewers"].get(), "CrossLinkSiteView"),
        "--ms", V["ms"].get(),
        "--psm", V["psm_xl"].get(),
        "--out", V["out"].get(),
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--inten", V["inten_thres"].get(),
        "--smooth", V["smooth"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_xesview():
    task.call(os.path.join(V["viewers"].get(), "ExhaustiveSearchViewX"),
        V["tg"].get(),
        "--ms", *(V["ms"].get().split(";")),
        "--out", V["out"].get(),
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run():
    if V["view"].get() == RView: run_view()
    elif V["view"].get() == CView: run_dualview()
    elif V["view"].get() == RXView: run_xview()
    elif V["view"].get() == CXView: run_xdualview()
    elif V["view"].get() == XSView: run_xsview()
    elif V["view"].get() == XESView: run_xesview()
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
util.add_entry(main, I, "URL:", V["url"], "New Port", new_port)
I += 1
util.add_entry(main, I, "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
I += 1
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I, columnspan=3)

f = F[RView]
I = 0
t = (("Target List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["fmt_target"], values=list(fmt_targets.keys()), state="readonly", justify="center"))
I += 1
t = (("UMZ file", "*.umz"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "Targeted MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=t))
I += 1
t = (("PSM", "*.csv"), ("PSM", "*.spectra"), ("All", "*.*"))
util.add_entry(f, I, "Targeted MS PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
t = (("UMZ file", "*.umz"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "Original MS Data:", V["ms_"], "Select", util.askfiles(V["ms_"], filetypes=t))
I += 1
util.add_entry(f, I, "Max. MS1 Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "FDR Threshold:", V["fdr"], "%")
I += 1
util.add_entry(f, I, "MS Sim. Thres.:", V["ms_sim_thres"])
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Optional").grid(column=0, row=I)
I += 1
t = (("Peptide Feature List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)

f = F[CView]
I = 0
ttk.Label(f, text=f"{CView} Not Available").grid(column=0, row=I, columnspan=3)

f = F[RXView]
I = 0
t = (("Target List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["fmt_target"], values=list(fmt_targets.keys()), state="readonly", justify="center"))
I += 1
t = (("UMZ file", "*.umz"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "Targeted MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=t))
I += 1
t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Targeted MS PSM:", V["psm_xl"], "Select", util.askfile(V["psm_xl"], filetypes=t))
I += 1
t = (("UMZ file", "*.umz"), ("MS2 file", "*.ms2"), ("All", "*.*"))
util.add_entry(f, I, "Original MS Data:", V["ms_"], "Select", util.askfiles(V["ms_"], filetypes=t))
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "Max. MS1 Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "FDR Threshold:", V["fdr"], "%")
I += 1
util.add_entry(f, I, "MS Sim. Thres.:", V["ms_sim_thres"])
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Optional").grid(column=0, row=I)
I += 1
t = (("Candidate XL List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Candidate XL List:", V["xl"], "Select", util.askfile(V["xl"], filetypes=t))
I += 1
t = (("Peptide Feature List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=t))
I += 1
t = (("Linear PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Linear PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)

f = F[CXView]
I = 0
t = (("Target List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfiles(V["tg"], V["out"], filetypes=t))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["fmt_target"], values=list(fmt_targets.keys()), state="readonly", justify="center"))
I += 1
t = (("Candidate XL List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Candidate XL List:", V["xl"], "Select", util.askfile(V["xl"], filetypes=t))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Data A").grid(column=0, row=I)
I += 1
t = (("UMZ file", "*.umz"), ("MS2 file", "*.ms2"), ("All", "*.*"))
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
t = (("UMZ file", "*.umz"), ("MS2 file", "*.ms2"), ("All", "*.*"))
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
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "MS1 Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "FDR:", V["fdr"], "%")

f = F[XSView]
I = 0
t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "XL PSM:", V["psm_xl"], "Select", util.askfile(V["psm_xl"], V["out"], filetypes=t))
I += 1
t = (("UMZ file", "*.umz"), ("MS1 file", "*.ms1"), ("All", "*.*"))
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=t))
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "Max. Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "Intensity Thres.:", V["inten_thres"])
I += 1
util.add_entry(f, I, "Smooth Width:", V["smooth"])

f = F[XESView]
I = 0
t = (("List", "*.csv"), ("All", "*.*"))
util.add_entry(f, I, "Peptide List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=t))
I += 1
t = (("UMZ file", "*.umz"), ("MS1 file", "*.ms1"), ("All", "*.*"))
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], filetypes=t))
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "Max. Mass Error:", V["error"], "ppm")

select_view(None, False)
