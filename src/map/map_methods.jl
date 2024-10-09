#=

Non-arithmetic functions acting on maps.

=#

# --- jacobian/jacobiant --- 
jacobian(m::TaylorMap;include_params=false) = GTPSA.jacobian(view(m.x, 1:numvars(m)),include_params=include_params)
jacobiant(m::TaylorMap;include_params=false) = GTPSA.jacobiant(view(m.x, 1:numvars(m)), include_params=include_params)

# --- checksymp ---
checksymp(m::TaylorMap) = checksymp(GTPSA.jacobian(m))

# --- setmatrix! ---
function setmatrix!(m::TaylorMap, M::AbstractMatrix)
  Base.require_one_based_indexing(M)
  nv = numvars(m)

  nv >= size(M,1) || error("Number of rows in matrix > number of variables in GTPSA!")

  for i=1:size(M,1)
    for j=1:size(M,2)
      m.x[i][j] = M[i,j]
    end
  end
end