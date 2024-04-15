var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = NonlinearNormalForm","category":"page"},{"location":"#NonlinearNormalForm","page":"Home","title":"NonlinearNormalForm","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NonlinearNormalForm.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NonlinearNormalForm]","category":"page"},{"location":"#NonlinearNormalForm.DAMap","page":"Home","title":"NonlinearNormalForm.DAMap","text":"DAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}\n\nTaylorMap that composes and inverses as a DAMap (with the scalar part ignored). See TaylorMap for more information.\n\n\n\n\n\n","category":"type"},{"location":"#NonlinearNormalForm.DAMap-Union{Tuple{Any}, Tuple{V}, Tuple{U}, Tuple{S}} where {S, U<:Union{Nothing, Quaternion{<:Union{ComplexTPS, TPS}}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.DAMap","text":"DAMap(M; x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nM must represent a matrix with linear indexing.\n\nConstructs a DAMap with the passed matrix of scalars M as the linear part of the TaylorMap, and optionally the entrance  coordinates x0, Quaternion for spin Q, and stochastic matrix E as keyword arguments. The helper keyword  arguments spin and radiation may be set to true to construct a DAMap with an identity quaternion/stochastic  matrix, or false for no spin/radiation. Note that setting spin/radiation to any Bool value without Q or E  specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA  Descriptor.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.DAMap-Union{Tuple{Probe{S, T, U, V}}, Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}} where {S, T<:Union{ComplexTPS, TPS}, U<:Union{Nothing, Quaternion{T}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.DAMap","text":"DAMap(p::Probe{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nCreates a DAMap from the Probe, which must contain TPSs. \n\nIf use is not specified, then the same GTPSA Descriptor as p will be used. If use is  specified (could be another Descriptor, TaylorMap, or a Probe containing TPSs), then the  p promoted to a DAMap will have the same Descriptor as in use. The number of variables  and parameters must agree, however the orders may be different.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.DAMap-Union{Tuple{TaylorMap{S, T, U, V}}, Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}} where {S, T<:Union{ComplexTPS, TPS}, U<:Union{Nothing, Quaternion{T}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.DAMap","text":"DAMapTaylorMap{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nCreates a new copy of the passed TaylorMap as a DAMap. \n\nIf use is not specified, then the same GTPSA Descriptor as m will be used. If use is  specified (could be another Descriptor, TaylorMap, or a Probe containing TPSs), then the  copy of m as a new DAMap will have the same Descriptor as in use. The number of variables  and parameters must agree, however the orders may be different.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.DAMap-Union{Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}, Tuple{UndefInitializer, TaylorMap{S, T, U, V}}} where {S, T, U, V}","page":"Home","title":"NonlinearNormalForm.DAMap","text":"DAMap(u::UndefInitializer, m::TaylorMap{S,T,U,V}) where {S,T,U,V}\n\nCreates an undefined DAMap based on m (same Descriptor and with radiation/spin on/off.)\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.DAMap-Union{Tuple{}, Tuple{Vector{T}}, Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}} where {S, T<:Union{ComplexTPS, TPS}, U<:Union{Nothing, Quaternion{T}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.DAMap","text":"DAMap(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nConstructs a DAMap with the passed vector of TPS/ComplexTPS as the orbital ray, and optionally the entrance  coordinates x0, Quaternion for spin Q, and stochastic matrix E as keyword arguments. The helper keyword  arguments spin and radiation may be set to true to construct a DAMap with an identity quaternion/stochastic  matrix, or false for no spin/radiation. Note that setting spin/radiation to any Bool value without Q or E  specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA  Descriptor. The use kwarg may also be used to change the Descriptor of the TPSs, provided the number of variables  and parameters agree (orders may be different).\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.TPSAMap","page":"Home","title":"NonlinearNormalForm.TPSAMap","text":"TPSAMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}} <: TaylorMap{S,T,U,V}\n\nTaylorMap that composes and inverses as a TPSAMap (with the scalar part included). See TaylorMap for more information.\n\n\n\n\n\n","category":"type"},{"location":"#NonlinearNormalForm.TPSAMap-Union{Tuple{Any}, Tuple{V}, Tuple{U}, Tuple{S}} where {S, U<:Union{Nothing, Quaternion{<:Union{ComplexTPS, TPS}}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.TPSAMap","text":"TPSAMap(M; x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nM must represent a matrix with linear indexing.\n\nConstructs a TPSAMap with the passed matrix of scalars M as the linear part of the TaylorMap, and optionally the entrance  coordinates x0, Quaternion for spin Q, and stochastic matrix E as keyword arguments. The helper keyword  arguments spin and radiation may be set to true to construct a TPSAMap with an identity quaternion/stochastic  matrix, or false for no spin/radiation. Note that setting spin/radiation to any Bool value without Q or E  specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA  Descriptor.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.TPSAMap-Union{Tuple{Probe{S, T, U, V}}, Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}} where {S, T<:Union{ComplexTPS, TPS}, U<:Union{Nothing, Quaternion{T}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.TPSAMap","text":"TPSAMap(p::Probe{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nCreates a TPSAMap from the Probe, which must contain TPSs. \n\nIf use is not specified, then the same GTPSA Descriptor as p will be used. If use is  specified (could be another Descriptor, TaylorMap, or a Probe containing TPSs), then the  p promoted to a TPSAMap will have the same Descriptor as in use. The number of variables  and parameters must agree, however the orders may be different.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.TPSAMap-Union{Tuple{TaylorMap{S, T, U, V}}, Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}} where {S, T<:Union{ComplexTPS, TPS}, U<:Union{Nothing, Quaternion{T}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.TPSAMap","text":"TPSAMapTaylorMap{S,T,U,V}; use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nCreates a new copy of the passed TaylorMap as a TPSAMap. \n\nIf use is not specified, then the same GTPSA Descriptor as m will be used. If use is  specified (could be another Descriptor, TaylorMap, or a Probe containing TPSs), then the  copy of m as a new TPSAMap will have the same Descriptor as in use. The number of variables  and parameters must agree, however the orders may be different.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.TPSAMap-Union{Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}, Tuple{UndefInitializer, TaylorMap{S, T, U, V}}} where {S, T, U, V}","page":"Home","title":"NonlinearNormalForm.TPSAMap","text":"TPSAMap(u::UndefInitializer, m::TaylorMap{S,T,U,V}) where {S,T,U,V}\n\nCreates an undefined TPSAMap based on m (same Descriptor and with radiation/spin on/off.)\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.TPSAMap-Union{Tuple{}, Tuple{Vector{T}}, Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}} where {S, T<:Union{ComplexTPS, TPS}, U<:Union{Nothing, Quaternion{T}}, V<:Union{Nothing, Matrix}}","page":"Home","title":"NonlinearNormalForm.TPSAMap","text":"TPSAMap(x::Vector{T}=zeros(TPS, numvars(GTPSA.desc_current)); x0::Vector{S}=zeros(numtype(first(x)), numvars(x)), Q::U=nothing, E::V=nothing,  spin::Union{Bool,Nothing}=nothing, radiation::Union{Bool,Nothing}=nothing, use::Union{Descriptor,TaylorMap,Probe{S,T,U,V},Nothing}=nothing) where {S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix,Nothing}}\n\nConstructs a TPSAMap with the passed vector of TPS/ComplexTPS as the orbital ray, and optionally the entrance  coordinates x0, Quaternion for spin Q, and stochastic matrix E as keyword arguments. The helper keyword  arguments spin and radiation may be set to true to construct a TPSAMap with an identity quaternion/stochastic  matrix, or false for no spin/radiation. Note that setting spin/radiation to any Bool value without Q or E  specified is type-unstable. This constructor also checks for consistency in the length of the orbital ray and GTPSA  Descriptor. The use kwarg may also be used to change the Descriptor of the TPSs, provided the number of variables  and parameters agree (orders may be different).\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.TaylorMap","page":"Home","title":"NonlinearNormalForm.TaylorMap","text":"TaylorMap{S,T<:Union{TPS,ComplexTPS},U<:Union{Quaternion{T},Nothing},V<:Union{Matrix{S},Nothing}}\n\nAbstract type for TPSAMap and DAMap used for normal form analysis. \n\nAll TaylorMaps contain x0 and x as the entrance coordinates and transfer map  as a truncated power series respectively. If spin is included, a field Q containing  a Quaternion as a truncated power series is included, else Q is nothing. If  radiation is included, a field E contains a matrix of the envelope for stochastic  radiation, else E is nothing.\n\nFields\n\nx0 – Entrance coordinates of the map, Taylor expansion point\nx  – Orbital ray as a truncated power series, expansion around x0 + scalar part equal to EXIT coordinates of map\nQ  – Quaternion as a truncated power series if spin is included, else nothing\nE  – Matrix of the envelope for stochastic radiation if included, else nothing\n\n\n\n\n\n","category":"type"},{"location":"#Base.inv-Tuple{Vector{<:Union{ComplexTPS, TPS}}}","page":"Home","title":"Base.inv","text":"inv(ma::Vector{<:Union{TPS,ComplexTPS}})\n\nInverts the map ma such that ma ∘ inv(ma) = 1 in the variables.\n\nExample\n\n\njulia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];\n\njulia> time = 0.01; k = 2; m = 0.01;\n\njulia> h = p^2/(2m) + 1/2*k*x^2;\n\njulia> hf = getvectorfield(h);\n\njulia> map = exppb(-time*hf, [x, p]);\n\njulia> map ∘ inv(map)\n2-element Vector{TPS}:\n  Out  Coefficient                Order   Exponent\n-------------------------------------------------\n   1:   1.0000000000000000e+00      1      1   0\n-------------------------------------------------\n   2:   1.0000000000000002e+00      1      0   1\n\n\n\n\n\n","category":"method"},{"location":"#Base.inv-Union{Tuple{TaylorMap{S, T, U, V}}, Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}} where {S, T, U, V}","page":"Home","title":"Base.inv","text":"inv(m1::TaylorMap{S,T,U,V}; dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_inv_work_low(m1)) where {S,T,U,V}\n\nInverts the TaylorMap.\n\nKeyword Arguments\n\ndospin – Specify whether to invert the quaternion as well, default is true\nwork_low – Temporary vector to hold the low-level C pointers. Default is output from prep_inv_work_low\n\n\n\n\n\n","category":"method"},{"location":"#GTPSA.compose!-Tuple{DAMap, DAMap, DAMap}","page":"Home","title":"GTPSA.compose!","text":"compose!(m::DAMap, m2::DAMap, m1::DAMap; dospin::Bool=true, keep_scalar::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_comp_work_low(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=prep_comp_work_prom(m,m2,m1))\n\nIn-place DAMap composition which calculates m = m2 ∘ m1, ignoring the scalar part of m1.\n\nAliasing m === m2 is allowed, but NOT m === m1. m2 === m1 is fine. The destination map m should be  properly set up (with correct types promoted if necessary), and have m.x[1:nv] (and m.Q.q if spin)  containing allocated TPSs.\n\nIf keep_scalar is set to true, then the scalar part of m1 is retained. In this case, a temporary  vector work_ref must be used to store the scalar part before calling the low-level compose_it!,  then added back into m1 after composition. work_ref can be optionally provided, or will be created  internally in this case. Default is true.\n\nIf dospin is true, then the quaternion part of the maps will be composed as well. Default is true.\n\nSee the documentation for compose_it! for information on work_low and work_prom.\n\nKeyword Arguments\n\nkeep_scalar – Specify whether to keep the scalar part in m1 or throw it away. If true, a temporary vector storing the scalar part must be used. Default is true\nwork_ref – If keep_scalar is true, the temporary vector can be provided in this keyword argument. If nothing is provided, the temporary will be created internally. Default is nothing\ndospin – Specify whether or not to include the quaternions in the concatenation. Default is true\nwork_low – Temporary vector to hold the low-level C pointers. See the compose_it! documentation for more details. Default is output from prep_comp_work_low(m)\nwork_prom – Temporary vector of allocated ComplexTPSs when there is implicit promotion. See the compose_it! documentation for more details. Default is output from prep_comp_work_prom(m, m2, m1)\n\n\n\n\n\n","category":"method"},{"location":"#GTPSA.compose!-Tuple{TPSAMap, TPSAMap, TPSAMap}","page":"Home","title":"GTPSA.compose!","text":"compose!(m::TPSAMap, m2::TPSAMap, m1::TPSAMap; dospin::Bool=true, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_comp_work_low(m), work_prom::Union{Nothing,Tuple{Vararg{Vector{<:ComplexTPS}}}}=prep_comp_work_prom(m,m2,m1))\n\nIn-place TPSAMapcomposition which calculatesm = m2 ∘ m1, including the scalar part ofm1`.\n\nAliasing m === m2 is allowed, but NOT m === m1. m2 === m1 is fine. The destination map m should be  properly set up (with correct types promoted if necessary), and have m.x[1:nv] (and m.Q.q if spin)  containing allocated TPSs.\n\nIf dospin is true, then the quaternion part of the maps will be composed as well. Default is true.\n\nSee the documentation for compose_it! for information on work_low and work_prom.\n\nKeyword Arguments\n\ndospin – Specify whether or not to include the quaternions in the concatenation. Default is true\nwork_low – Temporary vector to hold the low-level C pointers. See the compose_it! documentation for more details. Default is output from prep_comp_work_low(m)\nwork_prom – Temporary vector of allocated ComplexTPSs when there is implicit promotion. See the compose_it! documentation for more details. Default is output from prep_comp_work_prom(m, m2, m1)\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.checksymp-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Home","title":"NonlinearNormalForm.checksymp","text":"checksymp(M)\n\nReturns tranpose(M)*S*M - S, where S is the skew-symmetric matrix  S = blkdiag([0 1; -1 0], ...). If M is symplectic, then the result should be a matrix  containing all zeros. The non-symplectic parts of the matrix can be identified  by those nonzero elements in the result.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.compose_it!-Tuple{DAMap, DAMap, DAMap}","page":"Home","title":"NonlinearNormalForm.compose_it!","text":"compose_it!(m, m2, m1; dospin=true, work_low=nothing, work_prom=nothing)\n\nLow level composition function, m = m2 ∘ m1. Aliasing m with m2 is allowed, but not m with m1. Assumes the destination map is properly set up (with correct types promoted if necessary), and  that m.x[1:nv] (and m.Q.q if spin) contain allocated TPSs. The parameters part of m.x (m.x[nv+1:nn])  does not need to contain allocated TPSs.\n\nFor all compositions, 5 temporary vectors must be generated that contain Ptr{RTPSA} or Ptr{CTPSA} for each TPS in the map (depending on output type), to pass to the low-level C composition function in GTPSA.  Three vectors are for the orbital part (m.x, m2.x, m1.x correspondingly referred in the code as  outx_low, m2x_low, and m1x_low) and two are for the spin part (m.Q.q and m2.Q.q correspondingly  referred in the code as outQ_low and m2Q_low).  These 5 temporaries can be optionally passed as a tuple  in work_low, and must satisfy the following requirements:\n\nwork_low[1] = outx_low   # Length >= nv\nwork_low[2] = m2x_low    # Length >= nv\nwork_low[3] = m1x_low    # Length >= nv+np\nwork_low[4] = outQ_low   # Length >= 4, could be = work_low[1] or work_low[2] if nv >= 4\nwork_low[5] = m2Q_low    # Length >= 4, could be = work_low[1] or work_low[2] if nv >= 4\n\nFurthermore, for the spin part both of outx_low and m2x_low could be reused for outQ_low and m2Q_low  if nv >= 4 , however m1x_low MUST NOT BE REUSED!\n\nIf promotion is occuring, then one of the maps is ComplexTPS and the other TPS, with output map ComplexTPS.  Note that the spin part is required to agree with the orbital part in terms of type by definition of the TaylorMap  struct. work_prom can optionally be passed as a tuple containing the temporary ComplexTPSs if promotion is occuring:\n\nIf eltype(m.x) != eltype(m1.x) (then m1 must be promoted): work_prom[1] = m1x_prom  # Length >= nv+np, Vector{ComplexTPS}\n\nelse if eltype(m.x) != eltype(m2.x) (then m2 must be promoted):\n\nwork_prom[1] = m2x_prom  # Length >= nv, Vector{ComplexTPS}\nwork_prom[2] = m2Q_prom  # Length >= 4, Vector{ComplexTPS}\n\nNote that the ComplexTPSs in the vector(s) must be allocated and have the same Descriptor.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.compose_it!-Tuple{TPSAMap, TPSAMap, TPSAMap}","page":"Home","title":"NonlinearNormalForm.compose_it!","text":"compose_it!(m, m2, m1; dospin=true, work_low=nothing, work_prom=nothing)\n\nLow level composition function, m = m2 ∘ m1. Aliasing m with m2 is allowed, but not m with m1. Assumes the destination map is properly set up (with correct types promoted if necessary), and  that m.x[1:nv] (and m.Q.q if spin) contain allocated TPSs. The parameters part of m.x (m.x[nv+1:nn])  does not need to contain allocated TPSs.\n\nFor all compositions, 5 temporary vectors must be generated that contain Ptr{RTPSA} or Ptr{CTPSA} for each TPS in the map (depending on output type), to pass to the low-level C composition function in GTPSA.  Three vectors are for the orbital part (m.x, m2.x, m1.x correspondingly referred in the code as  outx_low, m2x_low, and m1x_low) and two are for the spin part (m.Q.q and m2.Q.q correspondingly  referred in the code as outQ_low and m2Q_low).  These 5 temporaries can be optionally passed as a tuple  in work_low, and must satisfy the following requirements:\n\nwork_low[1] = outx_low   # Length >= nv\nwork_low[2] = m2x_low    # Length >= nv\nwork_low[3] = m1x_low    # Length >= nv+np\nwork_low[4] = outQ_low   # Length >= 4, could be = work_low[1] or work_low[2] if nv >= 4\nwork_low[5] = m2Q_low    # Length >= 4, could be = work_low[1] or work_low[2] if nv >= 4\n\nFurthermore, for the spin part both of outx_low and m2x_low could be reused for outQ_low and m2Q_low  if nv >= 4 , however m1x_low MUST NOT BE REUSED!\n\nIf promotion is occuring, then one of the maps is ComplexTPS and the other TPS, with output map ComplexTPS.  Note that the spin part is required to agree with the orbital part in terms of type by definition of the TaylorMap  struct. work_prom can optionally be passed as a tuple containing the temporary ComplexTPSs if promotion is occuring:\n\nIf eltype(m.x) != eltype(m1.x) (then m1 must be promoted): work_prom[1] = m1x_prom  # Length >= nv+np, Vector{ComplexTPS}\n\nelse if eltype(m.x) != eltype(m2.x) (then m2 must be promoted):\n\nwork_prom[1] = m2x_prom  # Length >= nv, Vector{ComplexTPS}\nwork_prom[2] = m2Q_prom  # Length >= 4, Vector{ComplexTPS}\n\nNote that the ComplexTPSs in the vector(s) must be allocated and have the same Descriptor.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.exppb","page":"Home","title":"NonlinearNormalForm.exppb","text":"exppb(F::Vector{<:Union{TPS,ComplexTPS}}, m::Vector{<:Union{TPS,ComplexTPS}}=vars(first(F)))\n\nCalculates exp(F⋅∇)m = m + F⋅∇m + (F⋅∇)²m/2! + .... If m is not provided, it is assumed  to be the identity. \n\nExample\n\njulia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];\n\njulia> time = 0.01; k = 2; m = 0.01;\n\njulia> h = p^2/(2m) + 1/2*k*x^2;\n\njulia> hf = getvectorfield(h);\n\njulia> map = exppb(-time*hf, [x, p])\n2-element Vector{TPS}:\n  Out  Coefficient                Order   Exponent\n-------------------------------------------------\n   1:   9.9001665555952290e-01      1      1   0\n   1:   9.9666999841313930e-01      1      0   1\n-------------------------------------------------\n   2:  -1.9933399968262787e-02      1      1   0\n   2:   9.9001665555952378e-01      1      0   1\n\n\n\n\n\n","category":"function"},{"location":"#NonlinearNormalForm.fgrad-Tuple{Vector{<:Union{ComplexTPS, TPS}}, Union{ComplexTPS, TPS}}","page":"Home","title":"NonlinearNormalForm.fgrad","text":"fgrad(F::Vector{<:Union{TPS,ComplexTPS}}, g::Union{TPS,ComplexTPS})\n\nCalculates F⋅∇g.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.gethamiltonian-Tuple{Vector{<:Union{ComplexTPS, TPS}}}","page":"Home","title":"NonlinearNormalForm.gethamiltonian","text":"gethamiltonian(F::Vector{<:Union{TPS,ComplexTPS}})\n\nAssuming the variables in the TPSA are canonically-conjugate, and ordered so that the canonically- conjugate variables are consecutive (q₁, p₁, q₂, p₂, ...), this function calculates the Hamiltonian  from a vector field F that can be obtained from a Hamiltonian (e.g. by getvectorfield). Explicitly,  ∫ F₁ dp₁ - ∫ F₂ dq₁ + ... + ∫ F₂ₙ₋₁ dpₙ - ∫ F₂ₙ dqₙ\n\nExample\n\njulia> d = Descriptor(2,10); x = vars();\n\njulia> h = (x[1]^2 + x[2]^2)/2\nTPS:\n Coefficient                Order   Exponent\n  5.0000000000000000e-01      2      2   0\n  5.0000000000000000e-01      2      0   2\n\n\njulia> F = getvectorfield(h)\n2-element Vector{TPS}:\n  Out  Coefficient                Order   Exponent\n-------------------------------------------------\n   1:  -1.0000000000000000e+00      1      0   1\n-------------------------------------------------\n   2:   1.0000000000000000e+00      1      1   0\n\n\njulia> gethamiltonian(F)\nTPS:\n Coefficient                Order   Exponent\n  5.0000000000000000e-01      2      2   0\n  5.0000000000000000e-01      2      0   2\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.getvectorfield-Tuple{Union{ComplexTPS, TPS}}","page":"Home","title":"NonlinearNormalForm.getvectorfield","text":"getvectorfield(h::Union{TPS,ComplexTPS})::Vector{<:typeof(h)}\n\nAssuming the variables in the TPSA are canonically-conjugate, and ordered so that the canonically- conjugate variables are consecutive (q₁, p₁, q₂, p₂, ...), calculates the vector field (Hamilton's  equations) from the passed Hamiltonian, defined as [∂h/∂p₁, -∂h/∂q₁, ...]\n\nExample\n\njulia> d = Descriptor(2,10); x = vars();\n\njulia> h = (x[1]^2 + x[2]^2)/2\nTPS:\n Coefficient                Order   Exponent\n  5.0000000000000000e-01      2      2   0\n  5.0000000000000000e-01      2      0   2\n\n\njulia> getvectorfield(h)\n2-element Vector{TPS}:\n  Out  Coefficient                Order   Exponent\n-------------------------------------------------\n   1:  -1.0000000000000000e+00      1      0   1\n-------------------------------------------------\n   2:   1.0000000000000000e+00      1      1   0\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.inv!-Union{Tuple{V}, Tuple{U}, Tuple{T}, Tuple{S}, Tuple{TaylorMap{S, T, U, V}, TaylorMap{S, T, U, V}}} where {S, T, U, V}","page":"Home","title":"NonlinearNormalForm.inv!","text":"inv!(m::TaylorMap{S,T,U,V}, m1::TaylorMap{S,T,U,V}; dospin::Bool=true, work_ref::Union{Nothing,Vector{<:Union{Float64,ComplexF64}}}=nothing, work_low::Tuple{Vararg{Vector{<:Union{Ptr{RTPSA},Ptr{CTPSA}}}}}=prep_inv_work_low(m1)) where {S,T,U,V}\n\nIn-place inversion of the TaylorMap setting m = inv(m1). Aliasing m === m1 is allowed, however  in this case a temporary vector must be used to store the scalar part of m1 prior to inversion so  that the entrance/exit coordinates of the map can be properly handled.\n\nKeyword Arguments\n\ndospin – Specify whether to invert the quaternion as well, default is true\nwork_ref – If m === m1, then a temporary vector must be used to store the scalar part. If not provided and m === m1, this temporary will be created internally. Default is nothing\nwork_low – Temporary vector to hold the low-level C pointers. Default is output from prep_inv_work_low\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.lb-Tuple{Vector{<:Union{ComplexTPS, TPS}}, Vector{<:Union{ComplexTPS, TPS}}}","page":"Home","title":"NonlinearNormalForm.lb","text":"lb(A::Vector{<:Union{TPS,ComplexTPS}}, F::Vector{<:Union{TPS,ComplexTPS}})\n\nComputes the Lie bracket of the vector functions A and F, defined over N variables as  Σᵢᴺ Aᵢ (∂F/∂xᵢ) - Fᵢ (∂A/∂xᵢ)\n\nExample\n\njulia> d = Descriptor(2,10); x = vars();\n\njulia> A = [-x[2], x[1]]\n2-element Vector{TPS}:\n  Out  Coefficient                Order   Exponent\n-------------------------------------------------\n   1:  -1.0000000000000000e+00      1      0   1\n-------------------------------------------------\n   2:   1.0000000000000000e+00      1      1   0\n\n\njulia> F = [-x[1]^2, 2*x[1]*x[2]]\n2-element Vector{TPS}:\n  Out  Coefficient                Order   Exponent\n-------------------------------------------------\n   1:  -1.0000000000000000e+00      2      2   0\n-------------------------------------------------\n   2:   2.0000000000000000e+00      2      1   1\n\n\njulia> lb(A,F)\n2-element Vector{TPS}:\n  Out  Coefficient                Order   Exponent\n-------------------------------------------------\n   1:   4.0000000000000000e+00      2      1   1\n-------------------------------------------------\n   2:   3.0000000000000000e+00      2      2   0\n   2:  -2.0000000000000000e+00      2      0   2\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.locate_modes!","page":"Home","title":"NonlinearNormalForm.locate_modes!","text":"locate_modes!(evecs, evals=nothing; sort=true, modes=nothing)\n\nFor each mode (every pair of rows in evecs), determines the eigenvector (each  column in evecs) with the largest norm in that mode. The eigenvectors must  have dimensionality of 2n, n ∈ ℤ and be next to each other pairs. \n\nIf sort is true, then the pairs of eigenvectors are sorted in-place. The eigenvalues will also be  sorted if evals is provided. If sort is false, then the modes vector is filled with a mapping  of the current eigenvector position to its mode (index of modes is the current eigenvector  position, and the element at that index in modes is the mode that it corresponds to). For  sort=true, modes == 1:n after calling this function.\n\nIf the mode locating is successful, true is returned. For unsuccessful mode locating,  false is returned, modes will contain junk, and the eigenvectors/values may be partially sorted.\n\nAssumes the eigenvectors (and eigenvalues if provided) are already in pairs (next to  each other) starting at the first column.\n\nInput:\n\nevecs  – 2n x 2n matrix of eigenvectors already in pairs\nevals  – Vector of length 2n containing the eigenvalues corresponding to the eigenvectors in evecs\nsort   – (Optional) kwarg to specify whether or not to sort the eigenvectors, default is true\n\nOutput:\n\nret    – Returns true if the mode locating is successful, false if otherwise\nmodes – (Optional) kwarg for eigvec -> mode mapping, length(evecs) == Int(size(evecs,1)/2)\n\n\n\n\n\n","category":"function"},{"location":"#NonlinearNormalForm.logpb","page":"Home","title":"NonlinearNormalForm.logpb","text":"logpb(mf::Vector{<:Union{TPS,ComplexTPS}}, mi::Vector{<:Union{TPS,ComplexTPS}}=vars(first(F)))\n\nGiven a final map mf and initial map mi, this function calculates the vector field F such that mf=exp(F⋅∇)mi. If mi is not provided, it is assumed to be the identity.\n\njulia> d = Descriptor(2,10); x = vars()[1]; p = vars()[2];\n\njulia> time = 0.01; k = 2; m = 0.01;\n\njulia> h = p^2/(2m) + 1/2*k*x^2;\n\njulia> hf = getvectorfield(h);\n\njulia> map = exppb(-time*hf);\n\njulia> logpb(map) == -time*hf\ntrue\n\n\n\n\n\n","category":"function"},{"location":"#NonlinearNormalForm.mat_eigen!-Tuple{Any}","page":"Home","title":"NonlinearNormalForm.mat_eigen!","text":"mat_eigen!(mat; sort=true, phase_modes=true)\n\nSame as mat_eigen, but mutates mat for speed. See the documentation for mat_eigen  for more details.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.mat_eigen-Tuple{Any}","page":"Home","title":"NonlinearNormalForm.mat_eigen","text":"mat_eigen(mat; sort=true, phase_modes=true)\n\nGiven a matrix mat with an even number of rows/cols, calculates the eigenvectors  and eigenvalues. For stable eigenmodes (complex conjugate eigenvector/eigenvalue pairs),  the eigenvectors are normalized so that vⱼ'*S*vⱼ = +im for odd j, and -im for even j.\n\nIf sort is true, then each eigenvector/eigenvalue pair will be sorted according to the mode  it best identifies with. A warning will be printed if the mode sorting fails. Mode sorting  will automatically fail if more than 1 mode is unstable. Default is true.\n\nIf phase_modes is true, then each stable mode will be multiplied by a phase factor to make  the first component of that mode (e.g. the x, y, or z component, NOT the px, py, or pz component)  be fully real. This both makes the eigenvector \"pretty\", and in this case of weak coupling, ensures  that vⱼ'*vⱼ₊₁ for j=(1, 2, 3) has a positive imaginary part for odd j so that the eigenvectors/ eigenvalues for odd j are associated with the positive tune and even j the negative tune. This  phase factor is harmless/useless to include in a highly-coupled matrix,\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.moveback_unstable!-Tuple{LinearAlgebra.Eigen}","page":"Home","title":"NonlinearNormalForm.moveback_unstable!","text":"moveback_unstable!(F::Eigen) -> Int\n\nMoves back eigenvectors with eigenvalues having a zero imaginary component to the end of  the values and vectors arrays in the Eigen struct, and returns the number of unstable  eigenvectors. \n\nNote that if more than 1 mode is unstable, the pair of eigenvectors corresponding to a mode  will not necessarily be next to each other at the end of the matrix.\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.normalize_eigenmode!","page":"Home","title":"NonlinearNormalForm.normalize_eigenmode!","text":"normalize_eigenmode!(evec_pair, eval_pair, phase_mode::Integer=-1)\n\nNormalizes the complex-conjugate eigenvector pair so that vⱼ'*S*vⱼ = +im  for j=1, and -im for j=2. This may involve swapping the complex-conjugate  eigenvectors/eigenvalues.\n\nIf phase_mode is set to a particular mode (1, 2, 3), then a phase will be multiplied by the  eigenvectors to make them \"pretty\"; the phase will make the first component of that mode (e.g.  the x, y, or z component) be fully real. This gives the eigenvectors a simple form in the  weakly-coupled case, with the product v₁'*v₂ having a positive imaginary part, and ensures the  odd numbered eigenvector/eigenvalues are associated with the tune and the even numbered ones have  the negative of the tune. In cases with a lot of coupling, there is no simple way to define the  positive or negative tune, and multiplying by this phase factor is harmless. See the \"Tunes From  One-Turn Matrix Eigen Analysis\" section in the Bmad manual for more details. \n\n\n\n\n\n","category":"function"},{"location":"#NonlinearNormalForm.pb-Tuple{Union{ComplexTPS, TPS}, Union{ComplexTPS, TPS}}","page":"Home","title":"NonlinearNormalForm.pb","text":"pb(f::Union{TPS, ComplexTPS}, g::Union{TPS, ComplexTPS})\n\nAssuming the variables in the TPSA are canonically-conjugate, and ordered so that the canonically- conjugate variables are consecutive (q₁, p₁, q₂, p₂, ...), computes the Poisson bracket  of the scalar functions f and g. The Poisson bracket of two functions {f, g} is defined as  Σᵢ (∂f/∂qᵢ)(∂g/∂pᵢ) - (∂g/∂qᵢ)(∂f/∂pᵢ).\n\nExamples\n\njulia> d = Descriptor(4,10);\n\njulia> x = vars(d);\n\njulia> f = (x[1]^2 + x[2]^2)/2 + (x[3]^2 + x[4]^2)/2;\n\njulia> pb(f,x[1])\nTPS:\n  Coefficient              Order     Exponent\n  -1.0000000000000000e+00    1        0    1    0    0\n\n\njulia> pb(f,x[2])\nTPS:\n  Coefficient              Order     Exponent\n   1.0000000000000000e+00    1        1    0    0    0\n\n\njulia> pb(f,x[3])\nTPS:\n  Coefficient              Order     Exponent\n  -1.0000000000000000e+00    1        0    0    0    1\n\n\njulia> pb(f,x[4])\nTPS:\n  Coefficient              Order     Exponent\n   1.0000000000000000e+00    1        0    0    1    0\n\n\n\n\n\n","category":"method"},{"location":"#NonlinearNormalForm.ptinv-Tuple{Vector{<:Union{ComplexTPS, TPS}}, Vector{<:Integer}}","page":"Home","title":"NonlinearNormalForm.ptinv","text":"ptinv(ma::Vector{<:Union{TPS,ComplexTPS}}, vars::Vector{<:Integer})\n\nPartially-inverts the map ma, inverting only the variables specified by index in vars.\n\n\n\n\n\n","category":"method"}]
}
