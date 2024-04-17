struct VectorField{T<:Union{TPS,ComplexTPS}, U<:Union{Quaternion{T},Nothing}}
  x::Vector{T}  
  Q::U           
end

