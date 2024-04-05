# Generic symplectic skew symmetric S matrix (size inferred from 
# other matrix) using SkewLinearAlgebra JMatrix
struct SymplecticS end

const S = SymplecticS()

(S::SymplecticS)(n::Integer) = JMatrix{Int8,+1}(n)

for t = (:DAMap, :TPSAMap)
@eval begin

end
end
