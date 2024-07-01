import os
import random
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

FlowR = "Regular Linear Peptide View"
FlowXL = "Regular Cross-linked Peptide View"
FlowDXL = "Comparative Cross-linked Peptide View"
FlowXLS = "Cross-link Site View"
FlowESXL = "Cross-link Exhaustive Search View"
flows = [FlowR, FlowXL, FlowDXL, FlowXLS, FlowESXL]
vars_spec = {
    "flow": {"type": tk.StringVar, "value": FlowXL},
    "out": {"type": tk.StringVar, "value": ""},
    "tg": {"type": tk.StringVar, "value": ""},
    "fmt_tg": {"type": tk.StringVar, "value": "Auto Detect"},
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
        "--fmt", meta.fmts_tg[V["fmt_tg"].get()],
        "--error", V["error"].get(),
        "--ms_sim_thres", V["ms_sim_thres"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pf"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_xl_view():
    task.call(os.path.join(V["viewers"].get(), "TargetXLView"),
        V["tg"].get(),
        "--ms", V["ms"].get(),
        "--ms_old", *(V["ms_"].get().split(";")),
        "--psm", V["psm_xl"].get(),
        "--out", V["out"].get(),
        "--xl", V["xl"].get(),
        "--ft", V["ft"].get(),
        "--psm_pf", V["psm"].get(),
        "--fmt", meta.fmts_tg[V["fmt_tg"].get()],
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--ms_sim_thres", V["ms_sim_thres"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--cfg_pf", V["cfg_pf"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_dual_xl_view():
    task.call(os.path.join(V["viewers"].get(), "TargetDualXLView"),
        V["tg"].get(),
        "--ms", V["ms"].get(), V["ms_"].get(),
        "--psm", V["psm_xl"].get(), V["psm_xl_"].get(),
        "--out", V["out"].get(),
        "--xl", V["xl"].get(),
        "--ft", V["ft"].get(), V["ft_"].get(),
        "--fmt", meta.fmts_tg[V["fmt_tg"].get()],
        "--linker", V["linker"].get(),
        "--error", V["error"].get(),
        "--fdr", V["fdr"].get(),
        "--cfg", V["cfg_pl"].get(),
        "--host", V["url"].get().split(":")[0],
        "--port", V["url"].get().split(":")[1],
    )

def run_cross_link_site_view():
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

def run_exhaustive_search_xl_view():
    task.call(os.path.join(V["viewers"].get(), "ExhaustiveSearchXLView"),
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
    if V["flow"].get() == FlowR: run_view()
    elif V["flow"].get() == FlowXL: run_xl_view()
    elif V["flow"].get() == FlowDXL: run_dual_xl_view()
    elif V["flow"].get() == FlowXLS: run_cross_link_site_view()
    elif V["flow"].get() == FlowESXL: run_exhaustive_search_xl_view()
    else: print("Unknown Flow Type:", V["flow"].get())

def new_port():
    p = str(random.randint(49152, 65535))
    host = V["url"].get().split(":")[0]
    V["url"].set(host + ":" + p)

new_port()

def select_flow(_=None, verbose=True):
    for f in flows: F[f].grid_forget()
    if verbose: print("selected flow type:", V["flow"].get())
    F[V["flow"].get()].grid(column=0, row=I_fs, columnspan=3, sticky="EW")

util.init_form(main)
I = 0
c = ttk.Combobox(main, textvariable=V["flow"], values=flows, state="readonly", justify="center")
c.bind("<<ComboboxSelected>>", select_flow)
c.grid(column=0, row=I, columnspan=3, sticky="EW", padx=16, pady=4)
I += 1
F = {v: util.init_form(ttk.Frame(main)) for v in flows}
I_fs = I
I += 1
util.add_entry(main, I, "URL:", V["url"], "New Port", new_port)
I += 1
util.add_entry(main, I, "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
I += 1
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I, columnspan=3)

f = F[FlowR]
I = 0
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=meta.filetype_tg))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["fmt_tg"], values=list(meta.fmts_tg.keys()), state="readonly", justify="center"))
I += 1
util.add_entry(f, I, "Targeted MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=meta.filetype_ms))
I += 1
util.add_entry(f, I, "Targeted MS PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm))
I += 1
util.add_entry(f, I, "Original MS Data:", V["ms_"], "Select", util.askfiles(V["ms_"], filetypes=meta.filetype_ms))
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
util.add_entry(f, I, "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=meta.filetype_prec))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)

f = F[FlowXL]
I = 0
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=meta.filetype_tg))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["fmt_tg"], values=list(meta.fmts_tg.keys()), state="readonly", justify="center"))
I += 1
util.add_entry(f, I, "Targeted MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=meta.filetype_ms))
I += 1
util.add_entry(f, I, "Targeted MS PSM:", V["psm_xl"], "Select", util.askfile(V["psm_xl"], filetypes=meta.filetype_psm_xl))
I += 1
util.add_entry(f, I, "Original MS Data:", V["ms_"], "Select", util.askfiles(V["ms_"], filetypes=meta.filetype_ms))
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
util.add_entry(f, I, "Candidate XL List:", V["xl"], "Select", util.askfile(V["xl"], filetypes=meta.filetype_xl))
I += 1
util.add_entry(f, I, "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=meta.filetype_prec))
I += 1
util.add_entry(f, I, "Linear PSM:", V["psm"], "Select", util.askfile(V["psm"], filetypes=meta.filetype_psm))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)

f = F[FlowDXL]
I = 0
util.add_entry(f, I, "Target List:", V["tg"], "Select", util.askfiles(V["tg"], V["out"], filetypes=meta.filetype_tg))
I += 1
util.add_entry(f, I, "List Format:", ttk.Combobox(f, textvariable=V["fmt_tg"], values=list(meta.fmts_tg.keys()), state="readonly", justify="center"))
I += 1
util.add_entry(f, I, "Candidate XL List:", V["xl"], "Select", util.askfile(V["xl"], filetypes=meta.filetype_xl))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Data A").grid(column=0, row=I)
I += 1
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=meta.filetype_ms))
I += 1
util.add_entry(f, I, "Feature List:", V["ft"], "Select", util.askfile(V["ft"], filetypes=meta.filetype_prec))
I += 1
util.add_entry(f, I, "XL PSM:", V["psm_xl"], "Select", util.askfile(V["psm_xl"], filetypes=meta.filetype_psm_xl))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
ttk.Label(f, text="Data B").grid(column=0, row=I)
I += 1
util.add_entry(f, I, "MS Data:", V["ms_"], "Select", util.askfile(V["ms_"], filetypes=meta.filetype_ms))
I += 1
util.add_entry(f, I, "Feature List:", V["ft_"], "Select", util.askfile(V["ft_"], filetypes=meta.filetype_prec))
I += 1
util.add_entry(f, I, "XL PSM:", V["psm_xl_"], "Select", util.askfile(V["psm_xl_"], filetypes=meta.filetype_psm_xl))
I += 1
ttk.Separator(f, orient=tk.HORIZONTAL).grid(column=0, row=I, sticky="EW", padx=12)
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "MS1 Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "FDR:", V["fdr"], "%")

f = F[FlowXLS]
I = 0
util.add_entry(f, I, "XL PSM:", V["psm_xl"], "Select", util.askfile(V["psm_xl"], V["out"], filetypes=meta.filetype_psm_xl))
I += 1
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfile(V["ms"], filetypes=meta.filetype_ms))
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "Max. Mass Error:", V["error"], "ppm")
I += 1
util.add_entry(f, I, "Intensity Thres.:", V["inten_thres"])
I += 1
util.add_entry(f, I, "Smooth Width:", V["smooth"])

f = F[FlowESXL]
I = 0
util.add_entry(f, I, "Peptide List:", V["tg"], "Select", util.askfile(V["tg"], V["out"], filetypes=meta.filetype_tg))
I += 1
util.add_entry(f, I, "MS Data:", V["ms"], "Select", util.askfiles(V["ms"], filetypes=meta.filetype_ms))
I += 1
util.add_entry(f, I, "Default Linker:", V["linker"])
I += 1
util.add_entry(f, I, "Max. Mass Error:", V["error"], "ppm")

select_flow(None, False)
