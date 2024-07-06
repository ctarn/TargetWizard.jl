# [Development](@id dev)

```@contents
Pages = ["dev.md"]
Depth = 2:2
```

## Install Julia
Please install Julia (version 1.9 or newer) from [https://julialang.org](https://julialang.org).

## Clone the Repos
Please clone [UniMZ.jl](https://github.com/UniMZ/UniMZ.jl), [UniMZUtil.jl](https://github.com/UniMZ/UniMZUtil.jl) and [TargetWizard.jl](https://github.com/ctarn/TargetWizard.jl) by:
```sh
git clone git@github.com:UniMZ/UniMZ.jl.git
git clone git@github.com:UniMZ/UniMZUtil.jl.git
git clone git@github.com:ctarn/TargetWizard.jl.git
```
or
```sh
git clone https://github.com/UniMZ/UniMZ.jl.git
git clone https://github.com/UniMZ/UniMZUtil.jl.git
git clone https://github.com/ctarn/TargetWizard.jl.git
```

## Instantiate Julia Enviroment
Please `cd` to the root folder of `TargetWizard.jl`:
```sh
cd TargetWizard.jl
```

Run the following command to register the dependencies:
```sh
julia --project=. util/instantiate.jl ../UniMZ.jl ../UniMZUtil.jl
```
You should adjust the paths acorrdingly if the packages are saved to another path.

## Compile the Project
Compile the project using:
```sh
julia util/complie.jl
```

The complied files would be located at `./tmp/{your platform}/`.

## Build GUI and Installer
First you should copy the app icon to the `tmp/shared` folder.
You can run:
```sh
mkdir -p tmp/shared
cp fig/TargetWizard.png tmp/shared/
```

Finally, please run the scripts based on your platform if you want to build the graphic user inerface and package the software:
```sh
sh util/build_linux.sh
# or 
sh util/build_macos.sh
# or
./util/build_windows.bat
```

Python, PyInstaller, and Tkinter are required to build the GUI.
You can also call TargetWizard from command line directly using the compiled files.

The packaged software would be located at `./tmp/release/`.
