module ISPhvdc

import CSV
import DataFrames
const _DF = DataFrames
import JuMP
import Memento
import PowerModels
const _PM = PowerModels
import InfrastructureModels
const _IM = InfrastructureModels
import PowerModelsACDC
const _PMACDC = PowerModelsACDC
import CbaOPF
import XLSX
import ExcelFiles
const _EF = ExcelFiles
import StringDistances
const _SD = StringDistances


# Create our module level logger (this will get precompiled)
const _LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(PowerModels)`
# NOTE: If this line is not included then the precompiled `_PM._LOGGER` won't be registered at runtime.
__init__() = Memento.register(_LOGGER)

include("core/get_input_data.jl")
include("core/hosting_capacity.jl")
include("core/opf_timeseries.jl")


end # module ISPhvdc
