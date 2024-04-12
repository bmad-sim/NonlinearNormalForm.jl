"""
    mat_eigen!(mat; sort=true)

Same as `mat_eigen`, but mutates `mat` for speed.
"""
function mat_eigen!(mat; sort=true)
  F = eigen!(mat)
  low_mat_eigen!(F,sort)
  return F
end

"""
    mat_eigen(mat; sort=true)

Given a square matrix `mat` with an even number of rows/cols, calculates the eigenvectors 
and eigenvalues. For complex conjugate pairs, the eigenvectors are normalized so that 
`conj(vⱼ)*S*vⱼ = +im` for odd `j`, and `-im` for even `j`. Unstable eigenvectors/values (with 
eigenvalues having a zero imaginary component) are moved to the end and normalized to unity.
"""
function mat_eigen(mat; sort=true)
  F = eigen(mat)
  low_mat_eigen!(F,sort)
  return F
end

function low_mat_eigen!(F, sort)
  # Move unstable modes to the end:
  num_unstable = moveback_unstable!(F)

  if sort
    # Attempt to locate modes, but do not sort yet
    # Must first normalize stable after locating modes before sorting
    nv = length(F.values)
    n_modes = Int(nv/2) # number of modes
    modes = Vector{Int}(undef, n_modes)
    success = locate_modes!(F.vectors, F.values, sort=false, modes=modes)

    # Now normalize the stable modes:
    @views for i=1:Int((nv-num_unstable)/2)
      normalize_eigenmode!(F.vectors[:,2*i-1:2*i], F.values[2*i-1:2*i])
    end


    if success
      # Now finally sort by modes:
      F.vectors[:,2*modes.-1] = F.vectors[:,1:2:nv]
      F.vectors[:,2*modes] = F.vectors[:,2:2:nv]
      F.values[2*modes.-1] = F.values[1:2:nv]
      F.values[2*modes] = F.values[2:2:nv]
    else
      @warn "Mode sorting of eigenvectors failed; eigenvectors in arbitrary order..."
    end
  else
    # Normalize the stable modes:
    nv = length(F.values)
    @views for i=1:Int((nv-num_unstable)/2)
      normalize_eigenmode!(F.vectors[:,2*i-1:2*i], F.values[2*i-1:2*i])
    end
  end

  return F
end

"""
    normalize_eigenmode!(evec_pair, eval_pair)

Normalizes the complex-conjugate eigenvector pair so that `vⱼ'*S*vⱼ = +im` 
for `j=1`, and `-im` for `j=2`, and ensures that `vⱼ'*vⱼ` has positive imaginary 
part for `j=1`, ensuring that the `j=1` eigenvector/eigenvalue is associated with 
the tune and `j=2` is associated with the negative of the tune
"""
function normalize_eigenmode!(evec_pair, eval_pair)
  Base.require_one_based_indexing(evec_pair, eval_pair)
  nv = size(evec_pair,1)
  n_modes = Int(nv/2)

  size(evec_pair,2) == 2 || error("Eigenvector pair not provided")
  length(eval_pair) == 2 || error("Eigenvalue pair not provided")

  fnorm = 0
  for i=1:n_modes # For each mode
    fnorm += 2*imag(conj(evec_pair[2*i-1,1])*evec_pair[2*i,1])
  end
  sgn = sign(fnorm)
  if sgn == -1
    # Flip eval_pair
    eval_pair[1] = conj(eval_pair[1])
    eval_pair[2] = conj(eval_pair[2])
    #eval_pair[1], eval_pair[2] = eval_pair[2], eval_pair[1]
    fnorm=1/sqrt(-fnorm)
  else
    fnorm=1/sqrt(fnorm)
  end
  for i=1:nv
    evec_pair[i,1] = (sgn*real(evec_pair[i,1])+im*imag(evec_pair[i,1]))*fnorm
    evec_pair[i,2] = (sgn*real(evec_pair[i,2])+im*imag(evec_pair[i,2]))*fnorm
  end
  return
end


"""
    moveback_unstable!(F::Eigen) -> Int

This function moves back eigenvectors with eigenvalues having a zero imaginary component 
to the end of the `values` and `vectors` arrays in the Eigen struct, and returns the number 
of unstable eigenvectors. 

Note that if more than 1 mode is unstable, the pair of eigenvectors corresponding to a mode 
will not necessarily be next to each other at the end of the matrix.
"""
function moveback_unstable!(F::Eigen)
  evecs = F.vectors
  evals = F.values
  Base.require_one_based_indexing(evecs,evals)

  cnt = 0
  nv = size(evals,1)

  for i=1:nv
    if imag(evals[i]) == 0
      cnt += 1
    else
      evals[i-cnt], evals[i] = evals[i], evals[i-cnt]
      for k=1:nv
        evecs[k,i-cnt], evecs[k,i] = evecs[k,i], evecs[k,i-cnt]
      end
    end
  end

  return cnt
end

"""
    locate_modes!(evecs, evals=nothing; sort=true, modes=nothing)
  
For each mode (every pair of rows in `evecs`), determines the eigenvector (each 
column in `evecs`) with the largest `norm` in that mode. The eigenvectors must 
have dimensionality of 2n, n ∈ ℤ and be next to each other pairs.

If `sort` is true, then the pairs of eigenvectors are sorted in-place. The eigenvalues will also be 
sorted if `evals` is provided. If `sort` is false, then the `modes` vector is filled with a mapping 
of the current eigenvector position to its mode (index of `modes` is the current eigenvector 
position, and the element at that index in `modes` is the  mode that it corresponds to). For 
`sort=true`, `modes == 1:n` after calling this function.

If the mode locating is successful, `true` is returned. For unsuccessful mode locating, 
`false` is returned, `modes` will contain junk, and the eigenvectors/values may be partially sorted.

Assumes the eigenvectors (and eigenvalues if provided) are already in pairs (next to 
each other) starting at the first column.

### Input:
- `evecs`  -- 2n x 2n matrix of eigenvectors already in pairs
- `evals`  -- Vector of length 2n containing the eigenvalues corresponding to the eigenvectors in `evecs`
- `sort`   -- (Optional) kwarg to specify whether or not to sort the eigenvectors, default is true

### Output:
- `ret`    -- Returns true if the mode locating is successful, false if otherwise
- `modes` -- (Optional) kwarg for eigvec -> mode mapping, `length(evecs) == Int(size(evecs,1)/2)`
"""
function locate_modes!(evecs, evals=nothing; sort=true, modes=nothing)
  if isnothing(evals)
    Base.require_one_based_indexing(evecs)
  else
    Base.require_one_based_indexing(evecs,evals)
  end

  # eigenvectors already in pairs
  # we just need to check which has the biggest component
  # nv = number variables
  nv = size(evecs,1)
  n_modes = Int(nv/2)

  if isnothing(modes)
    modes = Vector{Int}(undef, n_modes)
  else
    Base.require_one_based_indexing(modes)
    length(modes) == n_modes || "Incorrect size of modes vector given matrix dimensions"
  end
  modes .= -1

  # For each mode
  for i=1:n_modes
    maxamp = @views norm(evecs[2*i-1:2*i,1])
    maxidx = 1
    # For each eigenvector
    for j=2:n_modes
      amp = @views norm(evecs[2*i-1:2*i,2*j-1])
      if amp > maxamp
        maxamp = amp
        maxidx = j
      end
    end
    if modes[maxidx] < 0
      if sort
        # current mode = i
        # eigenvector that should be in that mode = maxidx
        if i != maxidx # then must sort
          if !isnothing(evals) # also sort eigvals
            evals[2*i-1], evals[2*maxidx-1] = evals[2*maxidx-1], evals[2*i-1]
            evals[2*i], evals[2*maxidx] = evals[2*maxidx], evals[2*i]
          end
          for j=1:nv
            evecs[j,2*i-1], evecs[j,2*maxidx-1] = evecs[j, 2*maxidx-1], evecs[j,2*i-1]
            evecs[j,2*i], evecs[j,2*maxidx] = evecs[j, 2*maxidx], evecs[j,2*i]
          end
        end
        modes[i] = i # sorted
      else
        modes[maxidx] = i
      end
    else
      return false
    end
  end

  return true
end
