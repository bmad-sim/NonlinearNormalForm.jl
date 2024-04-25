"""
Generic symplectic skew symmetric S matrix (size inferred from 
other matrix in operations) using SkewLinearAlgebra's `JMatrix`
"""
struct SymplecticS end

"""
Generic symplectic skew symmetric S matrix (size inferred from 
other matrix in operations) using SkewLinearAlgebra's `JMatrix`
"""
const S = SymplecticS()

(S::SymplecticS)(n::Integer) = JMatrix{Int8,+1}(n)

for op = (:+, :-, :*, :/)
@eval begin
Base.$op(S::SymplecticS,M) = Base.$op(JMatrix{Int8,+1}(size(M,1)), M)
Base.$op(M,S::SymplecticS) = Base.$op(M, JMatrix{Int8,+1}(size(M,2)))
end
end

"""
    checksymp(M)

Returns `tranpose(M)*S*M - S`, where `S` is the skew-symmetric matrix 
`S = blkdiag([0 1; -1 0], ...)`. If `M` is symplectic, then the result should be a matrix 
containing all zeros. The non-symplectic parts of the matrix can be identified 
by those nonzero elements in the result.
"""
function checksymp(M)
  res = transpose(M)*S*M-S
  return res
end