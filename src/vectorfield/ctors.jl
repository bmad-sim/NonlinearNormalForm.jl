struct VectorField{T<:Union{TPS,ComplexTPS}, U<:Union{Quaternion{T},Nothing}}
  x::Vector{T}  
  Q::U           
end


"""
    VectorField(h::T; Q::U=nothing, spin::Union{Bool,Nothing}=nothing) where {T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing}}

Constructs a `VectorField` from the passed Hamiltonian `h`. Explicity, 
for `h`, constructs a vector field `F` such that
  
`F.x = [∂h/∂p₁, -∂h/∂q₁, ...]`
"""
function VectorField(h::T; Q::U=nothing, spin::Union{Bool,Nothing}=nothing) where {T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing}}

end