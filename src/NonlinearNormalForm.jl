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
             complex,
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
              cutord!,
              clear!,
              getdesc,
              numvars,
              numparams,
              numnn

export TaylorMap, Quaternion, Probe, TPSAMap, DAMap, TPSAMap, checksymp, mul!,
        normalize!, dot, to_SO3, read_fpp_map, cut, cut!, gofix, gofix!, normal, test, compose!, I


include("quaternion.jl")
include("probe.jl")
include("map.jl")
include("work.jl")
include("compose_it.jl")
include("operators.jl")
include("show.jl")
include("methods.jl")
include("normal.jl")
include("utils.jl")

end
