import os
from pathlib import Path
import tkinter as tk

import util

name = "TargetWizard"
version = "1.0.0"
author = "Tarn Yeong Ching"
url = f"http://{name.lower()}.ctarn.io"
server = f"http://api.ctarn.io/{name}/{version}"
copyright = f"{name} {version}\nCopyright © 2024 {author}\n{url}"
homedir = os.path.join(Path.home(), f".{name}", f"v{'.'.join(version.split('.')[0:2])}")

os.makedirs(homedir, exist_ok=True)

win = tk.Tk()
win.title(name)
win.iconphoto(True, tk.PhotoImage(file=util.get_content(f"{name}.png", shared=True)))

if util.is_darwin:
    path_mono = "/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono"
else:
    path_mono = "mono"

vars_spec = {
    "cfg_pl": {"type": tk.StringVar, "value": ""},
    "cfg_pf": {"type": tk.StringVar, "value": ""},
    "targetselect": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetSelect")},
    "targetselectxl": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetSelectXL")},
    "targetbind": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetBind")},
    "generators": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
    "viewers": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
}
vars = {k: v["type"](value=v["value"]) for k, v in vars_spec.items()}

filetype_tg = (("Target List", "*.csv"), ("All", "*.*"))
filetype_ms = (("UMZ file", "*.umz"), ("MS2 file", "*.ms2"), ("All", "*.*"))
filetype_psm = (("PSM", "*.csv"), ("PSM", "*.tsv"), ("PSM", "*.spectra"), ("All", "*.*"))
filetype_psm_xl = (("XL PSM", "*.csv"), ("All", "*.*"))
filetype_prec = (("Peptide Precursor List", "*.csv"), ("All", "*.*"))
filetype_xl = (("Candidate XL List", "*.csv"), ("All", "*.*"))

fmts_tg = {"Auto Detect": "auto", "TargetWizard": "TW", "Thermo Q Exactive": "TmQE", "Thermo Fusion": "TmFu"}
