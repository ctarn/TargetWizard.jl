import os
import random
import threading
import tkinter as tk
from tkinter import ttk, filedialog

import meta
import util

handles = []
running = False
skip_rest = False

path_autosave = os.path.join(meta.homedir, "autosave_view.task")

main = ttk.Frame()
main.pack(fill="both")

RView = "Regular Linear Peptide View"
CView = "Comparative Linear Peptide View"
RXView = "Regular Cross-linked Peptide View"
CXView = "Comparative Cross-linked Peptide View"
views = [RView, CView, RXView, CXView]
vars_spec = {
    "view": {"type": tk.StringVar, "value": RXView},
    "data": {"type": tk.StringVar, "value": ""},
    "ft": {"type": tk.StringVar, "value": ""},
    "tg": {"type": tk.StringVar, "value": ""},
    "xl": {"type": tk.StringVar, "value": ""},
    "psm": {"type": tk.StringVar, "value": ""},
    "psm_xl": {"type": tk.StringVar, "value": ""},
    "fdr": {"type": tk.StringVar, "value": "Inf"},
    "error": {"type": tk.StringVar, "value": "20.0"},
    "linker": {"type": tk.StringVar, "value": "DSSO"},
    "out": {"type": tk.StringVar, "value": ""},
    "targetview": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
    "cfg_pl": {"type": tk.StringVar, "value": ""},
    "cfg_pf": {"type": tk.StringVar, "value": ""},
    "url": {"type": tk.StringVar, "value": "127.0.0.1:30030"},
}
vars = {k: v["type"](value=v["value"]) for k, v in vars_spec.items()}
util.load_task(path_autosave, vars)

row = 0
util.init_form(main)

def select_view(_=None, verbose=True):
    for v in views:
        F[v].grid_forget()
    if verbose: print("selected view type:", vars["view"].get())
    F[vars["view"].get()].grid(column=0, row=row_fs, columnspan=3, sticky="EW")

c = ttk.Combobox(main, textvariable=vars["view"], values=views, state="readonly", justify="center")
c.bind("<<ComboboxSelected>>", select_view)
c.grid(column=0, row=row, columnspan=3, sticky="EW", padx=16, pady=4)
row += 1

F = {v: util.init_form(ttk.Frame(main)) for v in views}
row_fs = row
row += 1

util.add_entry(main, row, "Output Directory:", vars["out"], "Select", util.askdir(vars["out"]))
row += 1

def run_view():
    pass

def run_dualview():
    pass

def run_xview():
    cmd = [
        os.path.join(vars["targetview"].get(), "TargetXView"),
        "--error", vars["error"].get(),
        "--cfg", vars["cfg_pl"].get(),
        "--cfg_pf", vars["cfg_pf"].get(),
        "--linker", vars["linker"].get(),
        "--xl", vars["xl"].get(),
        "--ft", vars["ft"].get(),
        "--tg", vars["tg"].get(),
        "--psm", vars["psm_xl"].get(),
        "--psm_pf", vars["psm"].get(),
        "--fdr", vars["fdr"].get(),
        "--host", vars["url"].get().split(":")[0],
        "--port", vars["url"].get().split(":")[1],
        vars["data"].get(),
    ]
    util.run_cmd(cmd, handles, skip_rest)

def run_xdualview():
    pass

def do_load():
    path = filedialog.askopenfilename(filetypes=(("Configuration", "*.task"), ("All", "*.*")))
    if len(path) > 0: util.load_task(path, vars)

def do_autosave():
    util.save_task(path_autosave, {k: v for k, v in vars.items() if v.get() != vars_spec[k]["value"]})

def do_save():
    do_autosave()
    path = vars["out"].get()
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)
        util.save_task(os.path.join(path, "TargetView.task"), vars)
    else:
        print("`Output Directory` is required")

def do_run():
    btn_run.config(state="disabled")
    global handles, running, skip_rest
    running = True
    skip_rest = False
    do_save()
    if vars["view"].get() == RView:
        run_view()
    elif vars["view"].get() == CView:
        run_dualview()
    elif vars["view"].get() == RXView:
        run_xview()
    elif vars["view"].get() == CXView:
        run_xdualview()
    else:
        print("Unknown View Type:", vars["view"].get())
    running = False
    btn_run.config(state="normal")

def do_stop():
    global handles, running, skip_rest
    running = False
    skip_rest = True
    for job in handles:
        if job.poll() is None:
            job.terminate()
    handles.clear()
    btn_run.config(state="normal")
    print("TargetView stopped.")

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

util.add_entry(main, row, "TargetView:", vars["targetview"], "Select", util.askfile(vars["targetview"]))
row += 1

def do_new_port():
    p = str(random.randint(49152, 65535))
    host = vars["url"].get().split(":")[0]
    vars["url"].set(host + ":" + p)

do_new_port()
ttk.Label(main, text="URL:").grid(column=0, row=row, sticky="W")
ttk.Entry(main, textvariable=vars["url"]).grid(column=1, row=row, **util.sty_entry)
ttk.Button(main, text="New Port", command=do_new_port).grid(column=2, row=row, **util.sty_button)
row += 1

def do_select_data():
    filetypes = (("MS2", "*.ms2"), ("All", "*.*"))
    files = filedialog.askopenfilenames(filetypes=filetypes)
    if len(files) == 0:
        return None
    elif len(files) > 1:
        print("multiple data selected:")
        for file in files: print(">>", file)
    vars["data"].set(";".join(files))
    if len(vars["data"].get()) > 0 and len(vars["out"].get()) == 0:
        vars["out"].set(os.path.join(os.path.dirname(files[0]), "out"))

f = F[RView]
row = 0
ttk.Label(f, text=f"{RView} Not Available").grid(column=0, row=row, columnspan=3)
row += 1

f = F[CView]
row = 0
ttk.Label(f, text=f"{CView} Not Available").grid(column=0, row=row, columnspan=3)
row += 1

f = F[RXView]
row = 0
util.add_entry(f, row, "MS Data:", vars["data"], "Select", do_select_data)
row += 1

t = (("Peptide Feature List", "*.csv"), ("All", "*.*"))
util.add_entry(f, row, "Feature List:", vars["ft"], "Select", util.askfile(vars["ft"], filetypes=t))
row += 1

t = (("Target List", "*.csv"), ("All", "*.*"))
util.add_entry(f, row, "Target List:", vars["tg"], "Select", util.askfile(vars["tg"], filetypes=t))
row += 1

t = (("Candidate XL List", "*.csv"), ("All", "*.*"))
util.add_entry(f, row, "Candidate XL List:", vars["xl"], "Select", util.askfile(vars["xl"], filetypes=t))
row += 1

t = (("XL PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, row, "XL PSM:", vars["psm_xl"], "Select", util.askfile(vars["psm_xl"], filetypes=t))
row += 1

t = (("Linear PSM", "*.csv"), ("All", "*.*"))
util.add_entry(f, row, "Linear PSM:", vars["psm"], "Select", util.askfile(vars["psm"], filetypes=t))
row += 1

util.add_entry(f, row, "FDR:", vars["fdr"], "%")
row += 1

util.add_entry(f, row, "MS1 Mass Error:", vars["error"], "ppm")
row += 1

util.add_entry(f, row, "Default Linker:", vars["linker"])
row += 1

util.add_entry(f, row, "pLink Cfg. Directory:", vars["cfg_pl"], "Select", util.askdir(vars["cfg_pl"]))
row += 1

util.add_entry(f, row, "pFind Cfg. Directory:", vars["cfg_pf"], "Select", util.askdir(vars["cfg_pf"]))
row += 1

f = F[CXView]
row = 0
ttk.Label(f, text=f"{CXView} Not Available").grid(column=0, row=row, columnspan=3)
row += 1

select_view(None, False)
