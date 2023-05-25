import os
import threading
import tkinter as tk
from tkinter import ttk, filedialog

import meta
import util

handles = []
running = False
skip_rest = False

path_autosave = os.path.join(meta.homedir, "autosave_select.task")

main = ttk.Frame()
main.pack(fill="both")

if util.is_darwin:
    path_mono = "/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono"
else:
    path_mono = "mono"

tds = ["TT", "TD", "DD"]
pts = ["Inter", "Intra"]
fmts = ["TW", "TmQE", "TmFu"]
fmts_name = ["TargetWizard", "Thermo Q Exactive", "Thermo Fusion"]
vars_spec = {
    "data": {"type": tk.StringVar, "value": ""},
    "psm": {"type": tk.StringVar, "value": ""},
    "name": {"type": tk.StringVar, "value": "TargetWizard"},
    "fdr_min": {"type": tk.StringVar, "value": "-Inf"},
    "fdr_max": {"type": tk.StringVar, "value": "Inf"},
    "fdr_ge": {"type": tk.StringVar, "value": "≤"},
    "fdr_le": {"type": tk.StringVar, "value": "≤"},
    "batch": {"type": tk.StringVar, "value": "Inf"},
    "rtime": {"type": tk.StringVar, "value": "240"},
    "lc": {"type": tk.StringVar, "value": "Inf"},
    "out": {"type": tk.StringVar, "value": ""},
    "targetselect": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetSelect")},
    "thermorawread": {"type": tk.StringVar, "value": util.get_content("ThermoRawRead", "ThermoRawRead.exe", shared=True)},
    "mono": {"type": tk.StringVar, "value": path_mono},
}
for t in tds: vars_spec[f"td_{t}"] = {"type": tk.IntVar, "value": 1}
for t in pts: vars_spec[f"pt_{t}"] = {"type": tk.IntVar, "value": 1}
for t in fmts: vars_spec[f"fmt_{t}"] = {"type": tk.IntVar, "value": t == "TW"}
vars = {k: v["type"](value=v["value"]) for k, v in vars_spec.items()}
util.load_task(path_autosave, vars)

row = 0
util.init_form(main)

def do_select_data():
    filetypes = (("MS2", "*.ms2"), ("RAW", "*.raw"), ("All", "*.*"))
    files = filedialog.askopenfilenames(filetypes=filetypes)
    if len(files) == 0:
        return None
    elif len(files) > 1:
        print("multiple data selected:")
        for file in files: print(">>", file)
    vars["data"].set(";".join(files))
    if len(vars["data"].get()) > 0 and len(vars["out"].get()) == 0:
        vars["out"].set(os.path.join(os.path.dirname(files[0]), "out"))

util.add_entry(main, row, "MS Data:", vars["data"], "Select", do_select_data)
row += 1

t = (("PSM", "*.csv"), ("All", "*.*"))
util.add_entry(main, row, "PSM:", vars["psm"], "Select", util.askfile(vars["psm"], filetypes=t))
row += 1

util.add_entry(main, row, "Task Name:", vars["name"])
row += 1

_, f, _ = util.add_entry(main, row, "FDR Range:", ttk.Frame(main), "%")
ttk.Entry(f, textvariable=vars["fdr_min"]).pack(side="left", fill="x", expand=True)
ttk.Label(f, text="% ").pack(side="left")
ttk.Combobox(f, textvariable=vars["fdr_ge"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Label(f, text=" FDR ").pack(side="left")
ttk.Combobox(f, textvariable=vars["fdr_le"], values=("<", "≤"), state="readonly", justify="center", width=2).pack(side="left")
ttk.Entry(f, textvariable=vars["fdr_max"]).pack(side="left", fill="x", expand=True)
row += 1

_, f, _ = util.add_entry(main, row, "Target / Decoy Type:", ttk.Frame(main))
for t in tds: ttk.Checkbutton(f, text=t, variable=vars[f"td_{t}"]).pack(side="left", expand=True)
row += 1

_, f, _ = util.add_entry(main, row, "Inter- / Intra-Protein:", ttk.Frame(main))
for t in pts: ttk.Checkbutton(f, text=t, variable=vars[f"pt_{t}"]).pack(side="left", expand=True)
row += 1

util.add_entry(main, row, "Batch Size:", vars["batch"])
row += 1

util.add_entry(main, row, "RT Window:", vars["rtime"], "sec")
row += 1

util.add_entry(main, row, "LC Length:", vars["lc"], "min")
row += 1

_, f, _ = util.add_entry(main, row, "List Format:", ttk.Frame(main))
for t, n in zip(fmts, fmts_name): ttk.Checkbutton(f, text=n, variable=vars[f"fmt_{t}"]).pack(side="left", expand=True)
row += 1

util.add_entry(main, row, "Output Directory:", vars["out"], "Select", util.askdir(vars["out"]))
row += 1

def run_thermorawread(data, out):
    cmd = [vars["thermorawread"].get(), data, out]
    if not util.is_windows:
        cmd = [vars["mono"].get()] + cmd
    util.run_cmd(cmd, handles, skip_rest)
    return os.path.join(out, os.path.splitext(os.path.basename(data))[0] + ".ms1")

def run_targetselect(paths):
    cmd = [
        vars["targetselect"].get(),
        "--name", vars["name"].get(),
        "--psm", vars["psm"].get(),
        "--fdr_min", vars["fdr_min"].get(),
        "--fdr_max", vars["fdr_max"].get(),
        *(["--fdr_ge"] if vars["fdr_ge"].get() == "≤" else []),
        *(["--fdr_le"] if vars["fdr_le"].get() == "≤" else []),
        "--td", ",".join([t for t in tds if vars[f"td_{t}"].get()]),
        "--pt", ",".join([t for t in pts if vars[f"pt_{t}"].get()]),
        "--batch", vars["batch"].get(),
        "--rtime", vars["rtime"].get(),
        "--lc", vars["lc"].get(),
        "--fmt", ",".join([t for t in fmts if vars[f"fmt_{t}"].get()]),
        "--out", vars["out"].get(),
        *paths,
    ]
    util.run_cmd(cmd, handles, skip_rest)

def do_load():
    path = filedialog.askopenfilename(filetypes=(("Configuration", "*.task"), ("All", "*.*")))
    if len(path) > 0: util.load_task(path, vars)

def do_save():
    util.save_task(path_autosave, {k: v for k, v in vars.items() if v.get() != vars_spec[k]["value"]})
    path = vars["out"].get()
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)
        util.save_task(os.path.join(path, "TargetSelect.task"), vars)
    else:
        print("`Output Directory` is required")

def do_run():
    btn_run.config(state="disabled")
    global handles, running, skip_rest
    running = True
    skip_rest = False
    do_save()
    paths = []
    for p in vars["data"].get().split(";"):
        ext = os.path.splitext(p)[1].lower()
        if ext == ".raw":
            p = run_thermorawread(p, vars["out"].get())
        paths.append(p)
    run_targetselect(paths)
    running = False
    btn_run.config(state="normal")

def do_stop():
    global handles, running, skip_rest
    skip_rest = True
    for job in handles:
        if job.poll() is None:
            job.terminate()
    running = False
    handles.clear()
    btn_run.config(state="normal")
    print("TargetSelect stopped.")

frm_btn = ttk.Frame(main)
frm_btn.grid(column=0, row=row, columnspan=3)
ttk.Button(frm_btn, text="Load Task", command=do_load).grid(column=0, row=0, padx=16, pady=8)
ttk.Button(frm_btn, text="Save Task", command=do_save).grid(column=1, row=0, padx=16, pady=8)
btn_run = ttk.Button(frm_btn, text="Run Task", command=lambda: threading.Thread(target=do_run).start())
btn_run.grid(column=2, row=0, padx=16, pady=8)
ttk.Button(frm_btn, text="Stop Task", command=lambda: threading.Thread(target=do_stop).start()).grid(column=3, row=0, padx=16, pady=8)
row += 1

ttk.Separator(main, orient=tk.HORIZONTAL).grid(column=0, row=row, columnspan=3, sticky="EW")
ttk.Label(main, text="Advanced Configuration").grid(column=0, row=row, columnspan=3)
row += 1

util.add_entry(main, row, "TargetSelect:", vars["targetselect"], "Select", util.askfile(vars["targetselect"]))
row += 1

util.add_entry(main, row, "ThermoRawRead:", vars["thermorawread"], "Select", util.askfile(vars["thermorawread"]))
row += 1

if not util.is_windows:
    util.add_entry(main, row, "Mono Runtime:", vars["mono"], "Select", util.askfile(vars["mono"]))
    row += 1
