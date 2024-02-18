#= 
Simple differentiable (bi)quaternion bc Julia 
doesn't seem to have anything...
=#
struct Quaternion{T <: Number}
  q::NTuple{4, T}
end
