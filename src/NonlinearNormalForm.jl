module NonlinearNormalForm

import Base: ∘,
             *,
             literal_pow,
             +,
             -,
             /,
             \,
             ^,
             show,
             convert,
             inv,
             zero,
             one,
             ==

using LinearAlgebra,
      Printf,
      Reexport,
      DelimitedFiles

@reexport using GTPSA

import LinearAlgebra: norm, dot

import GTPSA: Desc, 
              RTPSA, 
              CTPSA,
              compose!,
              jacobian, 
              cutord,
              cutord!

export TaylorMap, Quaternion, Probe, TPSAMap, DAMap, TPSAMap, checksymp, mul!,
        normalize!, dot, to_SO3, read_fpp_map, cut, cut!, gofix, test, compose!


include("quaternion.jl")
include("probe.jl")
include("map.jl")
include("show.jl")
include("methods.jl")
include("utils.jl")

end
