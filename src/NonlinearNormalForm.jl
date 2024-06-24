module NonlinearNormalForm

import Base: ∘,
             *,
             +,
             -,
             /,
             \,
             ^,
             literal_pow,
             inv, 
             zero,
             zeros,
             one,
             complex,
             log,
             exp,
             ==,
             copy!,
             convert,
             show,
             rand,
             promote_rule,
             eltype

import LinearAlgebra: norm,
                      dot,
                      mul!

using LinearAlgebra,
      SkewLinearAlgebra,
      Printf,
      Reexport,
      DelimitedFiles
      
using ReferenceFrameRotations: Quaternion

@reexport using GTPSA

# We want these guys in our workspace:
import GTPSA: Desc, 
              RTPSA, 
              CTPSA,
              getdesc,
              numvars,
              numparams,
              numnn, 
              numtype,
              lowtype,

              add!,
              sub!,
              div!,

              jacobian,
              jacobiant,
              clear!,
              cutord,
              cutord!,
              getord,
              getord!,
              compose!             


export        TaylorMap, 
              Quaternion, 
              Probe,      
              TPSAMap, 
              DAMap, 
              VectorField,
              
              norm,
              dot,
              mul!,
              inv!,
              to_SO3,
      
              compose,
              inv!,
      
              checksymp,
              jacobian,
              jacobiant,
              I,
              S,
      
              read_fpp_map,
      
              mat_eigen,
              mat_eigen!,
              normalize_eigenmode!,
              locate_modes!,
              moveback_unstable!,
      
              normal!,
              normal,
              gofix!,
              gofix,
              linear_a!,
              linear_a,
              from_phasor!,
              from_phasor,
              to_phasor!,
              to_phasor,

              log!,
              pb,



              inv_with_log,
              equilibrium_moments











#include("utils/quaternion.jl")
include("types.jl")
include("utils/matrix.jl")
include("utils/symplectic_s.jl")
include("utils/gtpsa.jl")

include("probe.jl")
include("map/ctors.jl")
include("map/compose.jl")
include("map/compose_it.jl")

include("map/methods.jl")
include("map/operators.jl")

include("vectorfield/ctors.jl")
include("vectorfield/methods.jl")
include("vectorfield/operators.jl")
include("map/inv.jl")
include("work.jl")
include("normal.jl")

include("utils/misc.jl")
include("show.jl")






end
