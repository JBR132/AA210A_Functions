# AA210A_Functions
Compressible Flow Formulae Implemented in Julia

Material from "AA210A Course Reader: Fundamentals of Compressible Flow" by Prof. Brian J. Cantwell
https://web.stanford.edu/~cantwell/

## Usage
AA210A_Functions.jl is a Julia module, and can be loaded by either directly `include()`-ing it after opening the REPL as you would any other Julia script, like so:

    include("/...<filepath>.../AA210A_Functions.jl")

or by `using` it like a Julia module. For the latter method, the following line needs to be added to your `startup.jl` file or run in the REPL before `using` the module, since this will direct Julia to search your chosen directory for local modules rather than its internal module library.

    push!(LOAD_PATH,"/Users/JBR132/Documents/_Stanford/AA210A Comp Flow")

Once this has been done, you just need to call:

    using AA210A_Functions

and Julia should recognize and load the module. Once the module has been loaded, you can type `?AA210A_Functions` to see the module documentation, which includes a list of functions.

## Improvements
You can submit an issue if you want me to fix something or add another function, and I may or may not oblige.
