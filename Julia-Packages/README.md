# Plot3D Julia library

# Developing Packages
Useful Tutorials 
https://julialang.github.io/Pkg.jl/latest/managing-packages/#Adding-a-local-package-1

https://julialang.github.io/Pkg.jl/latest/managing-packages/#developing

## Some useful things that I discovered:
To add your own module for testing open up Julia inside your scripts folder or whereever you plan on running your test code. Type `julia` hit [enter] then `using Pkg` and `]` then `dev [relativeOrAbsolutePath]` to add your module to the `Manifest.toml` file

Whenever you add your own package, it gets added to the file `Manifest.toml`. This file is located in the Julia environments folder. Everytime you load up your module `julia>using MyModule` it will look in the `Manifest.toml` to get the location of where the module is stored. 

## Standard File Structure 
https://discourse.julialang.org/t/proper-way-of-organizing-code-into-subpackages/52835/4

```
    my_package_repo
    ├── Module1
    │   ├── Project.toml (optional)
    │   ├── src
    │   │   ├── Module1.jl
    │   │   ├── Module1SubFile1.jl 
    │   │   └── Module1SubFile2.jl
    │   └── test
    │       └── runtests.jl
    ├── Module2
    │   └src
    │       └── Module2.jl
    ├── Module3
    │   └src
    │       └── Module3.jl
    ...
    └── ModuleN
        └src
            └── ModuleN.jl
```