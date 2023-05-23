import os
import threading
import tkinter as tk
from tkinter import ttk, filedialog

import meta
import util

handles = []
running = False
skip_rest = False

path_autosave = os.path.join(meta.homedir, "autosave_report.task")

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
target_fmts = ["auto", "generic", "lumos", "hf"]
vars_spec = {
    "ms": {"type": tk.StringVar, "value": ""},
    "target": {"type": tk.StringVar, "value": ""},
    "target_fmt": {"type": tk.StringVar, "value": target_fmts[0]},
    "psm": {"type": tk.StringVar, "value": ""},
    "report": {"type": tk.StringVar, "value": BAReport},
    "out": {"type": tk.StringVar, "value": ""},
    "generators": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
    "thermorawread": {"type": tk.StringVar, "value": util.get_content("ThermoRawRead", "ThermoRawRead.exe", shared=True)},
    "mono": {"type": tk.StringVar, "value": path_mono},
}
vars = {k: v["type"](value=v["value"]) for k, v in vars_spec.items()}
util.load_task(path_autosave, vars)

row = 0
util.init_form(main)

def select_report(_=None, verbose=True):
    for r in reports:
        fs[r].grid_forget()
    if verbose: print("selected report type:", vars["report"].get())
    fs[vars["report"].get()].grid(column=0, row=row_fs, columnspan=3, sticky="EW")

c = ttk.Combobox(main, textvariable=vars["report"], values=reports, state="readonly", justify="center")
c.bind("<<ComboboxSelected>>", select_report)
c.grid(column=0, row=row, columnspan=3, sticky="EW", padx=16, pady=4)
row += 1

fs = {r: util.init_form(ttk.Frame(main)) for r in reports}
row_fs = row
row += 1

util.add_entry(main, row, "Output Directory:", vars["out"], "Select", util.askdir(vars["out"]))
row += 1

def run_thermorawread(data, out):
    cmd = [vars["thermorawread"].get(), data, out]
    if not util.is_windows:
        cmd = [vars["mono"].get()] + cmd
    util.run_cmd(cmd, handles, skip_rest)
    return os.path.join(out, os.path.splitext(os.path.basename(data))[0] + ".ms2")

def run_basicaquisitionreport():
    paths = []
    for p in vars["ms"].get().split(";"):
        ext = os.path.splitext(p)[1].lower()
        if ext == ".raw":
            p = run_thermorawread(p, vars["out"].get())
        paths.append(p)
    cmd = [os.path.join(vars["generators"].get(), "BasicAquisitionReport"), "--out", vars["out"].get(), *paths]
    util.run_cmd(cmd, handles, skip_rest)

def run_targetselectionreport():
    paths = vars["target"].get().split(";")
    cmd = [
        os.path.join(vars["generators"].get(), "TargetSelectionReport"),
        "--fmt", vars["target_fmt"].get(),
        "--out", vars["out"].get(),
        *paths,
    ]
    util.run_cmd(cmd, handles, skip_rest)

def run_targetaquisitionreport():
    pass

def run_identificationreport():
    pass

def do_load():
    path = filedialog.askopenfilename(filetypes=(("Configuration", "*.task"), ("All", "*.*")))
    if len(path) > 0: util.load_task(path, vars)

def do_save():
    util.save_task(path_autosave, {k: v for k, v in vars.items() if v.get() != vars_spec[k]["value"]})
    path = vars["out"].get()
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)
        util.save_task(os.path.join(path, "TargetReport.task"), vars)
    else:
        print("`Output Directory` is required")

def do_run():
    btn_run.config(state="disabled")
    global handles, running, skip_rest
    running = True
    skip_rest = False
    do_save()
    if vars["report"].get() == BAReport:
        run_basicaquisitionreport()
    elif vars["report"].get() == TSReport:
        run_targetselectionreport()
    elif vars["report"].get() == TAReport:
        run_targetaquisitionreport()
    elif vars["report"].get() == IDReport:
        run_identificationreport()
    else:
        print("Unknown Report Type:", vars["report"].get())
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
    print("TargetReport stopped.")

frm_btn = ttk.Frame(main)
frm_btn.grid(column=0, row=row, columnspan=3)
ttk.Button(frm_btn, text="Load Task", command=do_load).grid(column=0, row=0, padx=16, pady=8)
ttk.Button(frm_btn, text="Save Task", command=do_save).grid(column=1, row=0, padx=16, pady=8)
btn_run = ttk.Button(frm_btn, text="Run Task", command=lambda: threading.Thread(target=do_run).start())
btn_run.grid(column=2, row=0, padx=16, pady=8)
ttk.Button(frm_btn, text="Stop Task", command=lambda: threading.Thread(target=do_stop).start()).grid(column=3, row=0, padx=16, pady=8)
row += 1

ttk.Separator(main, orient=tk.HORIZONTAL).grid(column=0, row=row, columnspan=3, sticky="WE")
ttk.Label(main, text="Advanced Configuration").grid(column=0, row=row, columnspan=3)
row += 1

util.add_entry(main, row, "Generators:", vars["generators"], "Select", util.askdir(vars["generators"]))
row += 1

util.add_entry(main, row, "ThermoRawRead:", vars["thermorawread"], "Select", util.askfile(vars["thermorawread"]))
row += 1

if not util.is_windows:
    util.add_entry(main, row, "Mono Runtime:", vars["mono"], "Select", util.askfile(vars["mono"]))
    row += 1

def do_select_ms():
    filetypes = (("MS2", "*.ms2"), ("RAW", "*.raw"), ("All", "*.*"))
    files = filedialog.askopenfilenames(filetypes=filetypes)
    if len(files) == 0:
        return None
    elif len(files) > 1:
        print("multiple MS data selected:")
        for file in files: print(">>", file)
    vars["ms"].set(";".join(files))
    if len(vars["ms"].get()) > 0 and len(vars["out"].get()) == 0:
        vars["out"].set(os.path.join(os.path.dirname(files[0]), "out"))

def do_select_target():
    filetypes = (("Target List", "*.csv"), ("All", "*.*"))
    files = filedialog.askopenfilenames(filetypes=filetypes)
    if len(files) == 0:
        return None
    elif len(files) > 1:
        print("multiple target lists selected:")
        for file in files: print(">>", file)
    vars["target"].set(";".join(files))
    if len(vars["target"].get()) > 0 and len(vars["out"].get()) == 0:
        vars["out"].set(os.path.join(os.path.dirname(files[0]), "out"))

def do_select_psm():
    filetypes = (("PSM", "*.csv"), ("All", "*.*"))
    files = filedialog.askopenfilenames(filetypes=filetypes)
    if len(files) == 0:
        return None
    elif len(files) > 1:
        print("multiple PSM lists selected:")
        for file in files: print(">>", file)
    vars["psm"].set(";".join(files))
    if len(vars["psm"].get()) > 0 and len(vars["out"].get()) == 0:
        vars["out"].set(os.path.join(os.path.dirname(files[0]), "out"))

f = fs[BAReport]
row = 0
util.add_entry(f, row, "MS Data:", vars["ms"], "Select", do_select_ms)
row += 1

f = fs[TSReport]
row = 0
util.add_entry(f, row, "Target List:", vars["target"], "Select", do_select_target)
row += 1

util.add_entry(f, row, "Target Format:", ttk.Combobox(f, textvariable=vars["target_fmt"], values=target_fmts, state="readonly", justify="center"))
row += 1

f = fs[TAReport]
row = 0
ttk.Label(f, text=f"{TAReport} Not Available").grid(column=0, row=row, columnspan=3)
row += 1

util.add_entry(f, row, "Target List:", vars["target"], "Select", do_select_target)
row += 1

util.add_entry(f, row, "MS Data:", vars["ms"], "Select", do_select_ms)
row += 1

f = fs[IDReport]
row = 0
ttk.Label(f, text=f"{IDReport} Not Available").grid(column=0, row=row, columnspan=3)
row += 1

t = (("PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, row, "PSM:", vars["psm"], "Select", util.askfile(vars["psm"], filetypes=t))
row += 1

select_report(None, False)
