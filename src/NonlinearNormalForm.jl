module NonlinearNormalForm

import Base: ∘,
             *,
             literal_pow,
             +,
             -,
             ^,
             show,
             convert,
             inv,
             ==

using LinearAlgebra,
      Printf,
      Reexport,
      DelimitedFiles

@reexport using GTPSA

import LinearAlgebra: norm, dot

import GTPSA: numtype, 
              Desc, 
              RTPSA, 
              CTPSA,
              compose!,
              jacobian

export TaylorMap, Quaternion, Probe, TPSAMap, DAMap, TPSAMap, checksymp, qmul!,
        normalize!, dot, to_SO3, read_fpp_map, test


include("quaternion.jl")
include("probe.jl")
include("map.jl")
include("show.jl")
include("methods.jl")
include("utils.jl")

end
