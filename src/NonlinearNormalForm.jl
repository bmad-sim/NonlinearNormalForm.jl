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
      SkewLinearAlgebra,
      #StaticArrays,
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
              jacobiant, jacobiant!,
              cutord,
              cutord!,
              clear!,
              getdesc,
              numvars,
              numparams,
              numnn,
              numtype,
              lowtype

export TaylorMap, Quaternion, Probe, TPSAMap, DAMap, TPSAMap, checksymp, mul!,
        normalize!, dot, to_SO3, read_fpp_map, cutord, cutord!, gofix, gofix!, 
        normal, compose!, I, jacobian, jacobiant, linear_a, moveback_unstable!, normalize_evecs!,
        mat_eigen, mat_eigen!, locate_modes!, S, check_eigen, check_evecs_norm, inv!, testallocs!


include("quaternion.jl")
include("probe.jl")
include("map.jl")
include("work.jl")
include("compose_it.jl")
include("operators.jl")
include("show.jl")
include("methods.jl")
include("normal.jl")
include("matrix.jl")
include("symplectic_s.jl")
include("utils.jl")

end
