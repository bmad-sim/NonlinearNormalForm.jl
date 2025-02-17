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
             one,
             complex,
             real,
             imag,
             log,
             exp,
             ==,
             copy!,
             copy,
             convert,
             show,
             rand,
             promote_rule

import LinearAlgebra: norm

import TPSAInterface as TI
import TPSAInterface: AbstractTPSADef, getdef, nvars, nparams, ndiffs, maxord, nmonos
using LinearAlgebra,
      SkewLinearAlgebra,
      Printf,
      StaticArrays,
      DelimitedFiles
      
using ReferenceFrameRotations: Quaternion

export        TaylorMap, 
              Quaternion,    
              TPSAMap, 
              DAMap, 
              VectorField,
              
              norm,
              to_SO3,
      
              compose,
      
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
              equilibrium_moments,
              factorize




const COAST = eps(Float64)

macro _DEFAULT_X0(NV)
  return :(MVector{$(esc(NV))})
end

macro _DEFAULT_X(NN)
  return :(MVector{$(esc(NN))})
end

macro _DEFAULT_S(NV)
  return :(MMatrix{$(esc(NV)),$(esc(NV))})
end

#include("utils/quaternion.jl")

include("map.jl")
include("vectorfield.jl")
include("quaternion.jl")
include("utils.jl")
include("matrix.jl")
include("set.jl")
include("sanity.jl")
include("operators.jl")
include("methods.jl")


#include("staticarrays.jl")
#=
include("utils/matrix.jl")
include("utils/symplectic_s.jl")
include("sanity.jl")

include("methods.jl")
include("operators.jl")


include("map/compose.jl")
include("map/inv.jl")
include("map/map_methods.jl")
include("map/map_operators.jl")


include("vectorfield/ctors.jl")
include("vectorfield/vf_methods.jl")

include("work.jl")
include("normal.jl")

include("utils/misc.jl")
=#
#include("show.jl")






end
