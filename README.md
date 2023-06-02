# ISPhvdc.jl

This package provides the code to populate the 2000 bus NEM system with the PV, wind and demand traces coming from the ISP 2022.

## Usage:
You can clone the package and add it to your julia environment using 

```julia
] develop https://github.com/hakanergun/ISPhvdc.jl.git
```

You can add all required packages, such as CbaOPF, PowerModelsACDC, PowerModels etc. using

```julia
] add CbaOPF
```

Make sure that your Julia registry is up to date. To detect and download the latest version of the packages use

```julia
] update
```
## Examples:

There is an example script located under "scripts/NEM_2050.jl" where you can manipulate the input section for

- Selecting the ISP scenario and climate year
- Selecting if data should be dowloaded from AEMO (to avoid putting ~500 MB of data on the repository)  
    - Data only needs to be downloaded in the first run. .gitignore takes care it is not pushed to the repository.
- Selecting the (half) hours to run the calculations
- Selecting the OPF formulation, e.g. AC OPF, DC OPF or LPAC OPF 
- Selecting if parallel branches should be merged to a single branch
- Assigning solvers

The rest of the functions you don't necesserily need to manipulate. Once the OPFs are performed, there are some sample code snippets to analyse the results in more detail.

# To Do:

- Further validate results
- Create functions to run OPF in batch
- Write functions for better processing of results
- Apply clustering to operating hours to limit computation time for AC OPF
- Extend with WEM model if data becomes available