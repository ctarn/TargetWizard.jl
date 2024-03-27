import os
from pathlib import Path
import tkinter as tk

import util

name = "TargetWizard"
version = "0.0.0"
author = "Tarn Yeong Ching"
url = f"http://{name.lower()}.ctarn.io"
server = f"http://api.ctarn.io/{name}/{version}"
copyright = f"{name} {version}\nCopyright © 2023 {author}\n{url}"
homedir = os.path.join(Path.home(), f".{name}", f"v{'.'.join(version.split('.')[0:2])}")

os.makedirs(homedir, exist_ok=True)

win = tk.Tk()
win.title(name)
win.iconphoto(True, tk.PhotoImage(file=util.get_content(f"{name}.png", shared=True)))
win.resizable(False, False)

if util.is_darwin:
    path_mono = "/Library/Frameworks/Mono.framework/Versions/Current/Commands/mono"
else:
    path_mono = "mono"

vars_spec = {
    "cfg_pl": {"type": tk.StringVar, "value": ""},
    "cfg_pf": {"type": tk.StringVar, "value": ""},
    "targetselect": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetSelect")},
    "targetselectx": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetSelectX")},
    "targetbind": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin", "TargetBind")},
    "generators": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
    "viewers": {"type": tk.StringVar, "value": util.get_content("TargetWizard", "bin")},
}
vars = {k: v["type"](value=v["value"]) for k, v in vars_spec.items()}
