import os
import tkinter as tk
from tkinter import ttk

import meta
import util

main = ttk.Frame()
main.pack(fill="both")

PrecModeIC = "by isolation center"
PrecModeIW = "by isolation window"
PrecModeExIW = "by extended isolation window"
PrecModes = {PrecModeIC: "center", PrecModeIW: "window", PrecModeExIW: "extended_window"}
fmts = ["csv", "tsv", "ms2", "mgf", "pf2"]
vars_spec = {
    "ms": {"type": tk.StringVar, "value": ""},
    "tg": {"type": tk.StringVar, "value": ""},
    "fmt_tg": {"type": tk.StringVar, "value": "Auto Detect"},
    "out": {"type": tk.StringVar, "value": ""},
    "mode": {"type": tk.StringVar, "value": PrecModeExIW},
    "error_rt": {"type": tk.StringVar, "value": "16"},
    "error_mz": {"type": tk.StringVar, "value": "10.0"},
}
for fmt in fmts: vars_spec[f"fmt_{fmt}"] = {"type": tk.IntVar, "value": fmt in ["csv"]}
task = util.Task("TargetBind", vars_spec, path=meta.homedir, shared_vars_spec=meta.vars_spec, shared_vars=meta.vars)
V = task.vars

def run():
    task.call(V["targetbind"].get(), *(V["ms"].get().split(";")),
        "--target", V["tg"].get(),
        "--out", V["out"].get(),
        "--mode", PrecModes[V["mode"].get()],
        "--error_rt", V["error_rt"].get(),
        "--error_mz", V["error_mz"].get(),
        "--fmt_target", meta.fmts_tg[V["fmt_tg"].get()],
        "--fmt", ",".join(filter(lambda x: V[f"fmt_{x}"].get(), fmts)),
    )

util.init_form(main)
I = util.Counter()
util.add_entry(main, I.next(), "Targeted MS Data:", V["ms"], "Select", util.askfiles(V["ms"], V["out"], filetypes=meta.filetype_ms))
util.add_entry(main, I.next(), "Target List:", V["tg"], "Select", util.askfile(V["tg"], filetypes=meta.filetype_tg))
util.add_entry(main, I.next(), "List Format:", ttk.Combobox(main, textvariable=V["fmt_tg"], values=list(meta.fmts_tg.keys()), state="readonly", justify="center"))
util.add_entry(main, I.next(), "Binding Mode:", ttk.Combobox(main, textvariable=V["mode"], values=list(PrecModes.keys()), state="readonly", justify="center"))
util.add_entry(main, I.next(), "Max. RTime Error:", V["error_rt"], "sec")
util.add_entry(main, I.next(), "Max. Mass Error:", V["error_mz"], "ppm")
_, f, _ = util.add_entry(main, I.next(), "Output Format", ttk.Frame(main))
for x in fmts: ttk.Checkbutton(f, text=x.upper(), variable=V[f"fmt_{x}"]).pack(side="left", expand=True)
util.add_entry(main, I.next(), "Output Directory:", V["out"], "Select", util.askdir(V["out"]))
task.init_ctrl(ttk.Frame(main), run).grid(column=0, row=I.next(), columnspan=3)
