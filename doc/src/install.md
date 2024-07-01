# [Installation](@id install)

## Download

Please download the software from [https://targetwizard.ctarn.io](https://targetwizard.ctarn.io).

## Install

### Linux

Please unzip the downloaded `.zip` file, and TargetWizard can be used directly without installation.

### macOS

For macOS users, we provide both `.pkg` and `.zip` files.

We would recommend to use the `.pkg` file which can be installed by simply double clicking it.
The software would be installed at `/Applications/TargetWizard.app` by default.

The `.zip` file contains the `.app` package software and can be used directly without installation.
If the macOS says:

```
“TargetWizard.app” is damaged and can’t be opened. You should move it to the Trash.
```

Please run 
```sh
sudo xattr -r -d com.apple.quarantine [path/to/TargetWizard.app]
```
in terminal to remove the quarantine attributions.

### Windows

Please unzip the downloaded `.zip` file, and TargetWizard can be used directly without installation.

The software is packaged using PyInstaller,
and can be detected as virus by mistake on Windows (see [the issue](https://github.com/pyinstaller/pyinstaller/issues/5932)).
Please restore the deleted file from `Protection History`,
and Windows Security should not stop or delete it again.
Otherwise,
please add the software to white list.
You can also package the software from source yourself.
