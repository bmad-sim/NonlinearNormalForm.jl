"""
    mat_eigen!(evecs::Matrix{<:Number}, evals::Vector{<:Number}, mat::Matrix{<:Number}; keep_mat=true, sort=true)

Given a square matrix `mat` with an even number of rows/cols, calculates the eigenvectors 
and eigenvalues. For complex conjugate pairs, the eigenvectors are normalized so that 
`conj(vⱼ)*S*vⱼ = +im` for odd `j`, and `-im` for even `j`. 
"""
function mat_eigen!(mat::Matrix{<:Number}; sort=true) 
  F = eigen!(mat)
  if sort && !locate_planes!(F, sort=true)
    @warn "Plane sorting of eigenvectors failed; eigenvectors in arbitrary order..."
  end
  normalize_evecs!(F.vectors)
  return F
end

function mat_eigen(mat::Matrix{<:Number}; sort=true)
  F = eigen!(mat)
  if sort && !locate_planes!(F, sort=true)
    @warn "Plane sorting of eigenvectors failed; eigenvectors in arbitrary order..."
  end
  normalize_evecs!(F.vectors)
  return F
end


function pairup_eigen!(F::Eigen)
  vecs = F.vectors
  vals = F.vals
  nv = size(vals,1)
  npl = 

end


"""
    normalize_evecs!(evecs::Matrix{<:Complex})

Assuming the eigenvectors are in complex-conjugate pais, the eigenvectors are 
normalized so that `conj(vⱼ)*S*vⱼ = +im` for odd `j`, and `-im` for even `j` where 
`S` is the skew-symmetric matrix e.g. [0 1; -1 0] for 2D.
"""
function normalize_evecs!(evecs::Matrix{<:Complex})
  npl = Int(size(evecs,1)/2) # Number planes = Number harmonic variables / 2

  for i=1:npl # For each eigenvector pair
    fnorm = 0
    for j=1:nhpl # For each harmonic plane
      fnorm += 2*imag(conj(evecs[2*j-1,2*i-1])*evecs[2*j,2*i-1])
    end
    println(fnorm)
    @views evecs[:,2*i-1] = evecs[:,2*i-1]/sqrt(abs(fnorm))
    @views evecs[:,2*i] = evecs[:,2*i]/sqrt(abs(fnorm))
  end
end

"""
    locate_planes!(F::Union{Matrix,Eigen}; sort::Bool=true, planes::Union{Vector{Int},Nothing}=nothing)
  
For each plane (every pair of rows in `F` or `F.vectors`), determines the eigenvector 
(each column in `F` or `F.vectors`) with the largest `norm` in that plane. The eigenvectors 
must have dimensionality of 2n, n ∈ ℤ.

If `sort` is true, then the eigenvectors are sorted in-place. The eigenvalues will also be 
sorted if `F isa Eigen`. If `sort` is false, then the `planes` vector is filled with a mapping 
of the current eigenvector position to its plane (index of `planes` is the current eigenvector 
position, and the element at that index in `planes` is the  plane that it corresponds to). For 
`sort=true`, `planes == 1:n` after calling this function.

If the plane locating is successful, `true` is returned. For unsuccessful plane locating, 
`false` is returned, `planes` will contain junk, and the eigenvectors may be partially sorted.

Assumes the eigenvectors (and eigenvalues if `F isa Eigen`) in `F` are already in pairs (next to 
each other) starting at the first column.

### Input:
- `F`      -- 2n x 2n matrix of eigenvectors already in pairs, or an Eigen object 
- `sort`   -- (Optional) kwarg to specify whether or not to sort the eigenvectors, default is true

### Output:
- `ret`    -- Returns true if the plane locating is successful, false if otherwise
- `planes` -- (Optional) kwarg for eigvec -> plane mapping, `length(vecs) == Int(size(vecs,1)/2)`
"""
function locate_planes!(F::Union{Matrix,Eigen}; sort::Bool=true, planes::Union{Vector{Int},Nothing}=nothing)
  # eigenvectors already in pairs
  # we just need to check which has the biggest component
  # nv = number variables
  if F isa Eigen # Also sort eigenvals
    evecs = F.vectors
    evals = F.values
    nv = size(F.values,1)
  else 
    evecs = F
    nv = size(evecs,1)
  end
  
  npl = Int(nv/2)
  if isnothing(planes)
    planes = Vector{Int}(undef, npl)
  else
    length(planes) == npl || "Incorrect size of planes vector given matrix dimensions"
  end
  planes .= -1

  # For each plane
  for i=1:npl
    maxamp = @views norm(evecs[2*i-1:2*i,1])
    maxidx = 1
    # For each eigenvector
    for j=2:npl
      amp = @views norm(evecs[2*i-1:2*i,2*j-1])
      if amp > maxamp
        maxamp = amp
        maxidx = j
      end
    end
    if planes[maxidx] < 0
      if sort
        # current plane = i
        # eigenvector that should be in that plane = maxidx
        if i != maxidx # then must sort
          if F isa Eigen # also sort eigvals
            evals[2*i-1], evals[2*maxidx-1] = evals[2*maxidx-1], evals[2*i-1]
            evals[2*i], evals[2*maxidx] = evals[2*maxidx], evals[2*i]
          end
          for j=1:nv
            evecs[j,2*i-1], evecs[j,2*maxidx-1] = evecs[j, 2*maxidx-1], evecs[j,2*i-1]
            evecs[j,2*i], evecs[j,2*maxidx] = evecs[j, 2*maxidx], evecs[j,2*i]
          end
        end
        planes[i] = i # sorted
      else
        planes[maxidx] = i
      end
    else
      return false
    end
  end

  return true
end


"""
    checksymp(M::Matrix{T}) where T<:Number

Returns `tranpose(M)*S*M - S`, where `S` is the skew-symmetric matrix 
`S = [0 I; -I 0]`. If `M` is symplectic, then the result should be a matrix 
containing all zeros. The non-symplectic parts of the matrix can be identified 
by those nonzero elements in the result.
"""
function checksymp(M::Matrix{T}) where T<:Number
  s = size(M)
  nv = first(s)
  nv == last(s) || error("Non-square matrix!")
  iseven(nv) || error("Matrix contains odd number of rows/columns!")
  S = zeros(nv,nv)
  for i=1:2:nv
    S[i:i+1,i:i+1] = [0 1; -1 0];
  end
  res = transpose(M)*S*M-S
  return res
end
