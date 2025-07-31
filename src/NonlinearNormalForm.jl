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

import LinearAlgebra: norm, factorize

import TPSAInterface as TI
import TPSAInterface: AbstractTPSAInit, getinit, ndiffs, maxord, nmonos
using LinearAlgebra,
      SkewLinearAlgebra,
      Printf,
      StaticArrays,
      DelimitedFiles
      
using ReferenceFrameRotations: Quaternion, vect

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
              c_map,
              ci_map,
              c_jacobian,
              ci_jacobian,

              log!,
              pb,



              inv_with_log,
              equilibrium_moments,
              factorize,
              fast_canonize,
              compute_de_moivre,
              compute_sagan_rubin,

              LinearDAMap


# After experimenting I have found MVector
# to be the fastest for the entire NNF workflow,
# even though the cost of constructing a map 
# initially (before nv and np are known) is higher

macro _DEFAULT_X0(NV)
  return :(MVector{$(esc(NV))})
end

macro _DEFAULT_X(NN)
  return :(MVector{$(esc(NN))})
end

macro _DEFAULT_S(NV)
  return :(MMatrix{$(esc(NV)),$(esc(NV))})
end

const DEFAULT_NVARS = 6

#include("utils/quaternion.jl")

include("map.jl")
include("vectorfield.jl")
include("quaternion.jl")
include("utils.jl")
#include("show.jl")
include("matrix.jl")
include("set.jl")
include("sanity.jl")
include("operators.jl")
include("methods.jl")
include("normal.jl")
include("linear.jl")
include("de_moivre.jl")
include("sagan_rubin.jl")


end
