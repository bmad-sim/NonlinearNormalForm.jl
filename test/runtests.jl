using NonlinearNormalForm
using TPSAInterface
using TPSAInterface: TPSAInterface as TI
using Test, GTPSA

include("readfpp.jl")

@testset "Composition and inversion" begin
    d = Descriptor(1,2,1,2) 
    x1 = @vars(d)[1]
    m1 = DAMap(v=[1+2*x1+2*x1^2], v0=[4])
    m2 = DAMap(v=[1+2*x1+2*x1^2], v0=[3])

    mt1 = TPSAMap(m1)
    mt2 = TPSAMap(m2)
    
    tol = 1e-10

    @test norm(m2∘m1 - DAMap(v=[1+4*x1+12*x1^2], v0=[4])) < tol
    @test norm(mt2∘mt1 - TPSAMap(v=[5-12*x1-4*x1^2], v0=[4])) < tol
    @test norm(m2^2 - m2∘m2) < tol
    @test norm(mt2^2 - mt2∘mt2) < tol
    @test norm(m2^3 - m2∘m2∘m2) < tol
    @test norm(mt2^3 - mt2∘mt2∘mt2) < tol
    @test norm(m2^3 - DAMap(v=[1+8*x1+56*x1^2], v0=[3])) < tol
    @test norm(mt2^3 - TPSAMap(v=[13+8*x1-8*x1^2], v0=[3])) < tol

    @test norm(m1^-2 - inv(m1^2)) < tol
    @test norm(m1^-3 - DAMap(v=[4+0.125*x1-0.109375*x1^2],v0=[1])) < tol
    @test norm(m1^-3 - inv(m1^3)) < tol
    @test norm(mt1^-2 - inv(mt1^2)) < tol
    @test norm(mt1^-3 - inv(mt1^3)) < tol
    @test norm(mt1^-3 - TPSAMap(v=[4+0.9615384615384616E-02*x1-0.4978379608557124E-04*x1^2],v0=[-35])) < tol 
end

@testset "Comparison with FPP" begin

    tol = 1e-12

    # Normal form -------------------------------------
    # 1D all pseudo-harmonic oscillators order 3
    m = read_fpp_map("order2/test.map",spin=false)
    r_fpp = read_fpp_map("order2/R.map",spin=false)
    c = c_map(m)
    a = normal(m)
    @test norm(inv(c)∘inv(a)∘m∘a∘c - r_fpp) < tol

    # 1D all pseudo-harmonic oscillators order 10
    m = read_fpp_map("order10/test.map",spin=false)
    r_fpp = read_fpp_map("order10/R.map",spin=false)
    c = c_map(m)
    a = normal(m)
    @test norm(inv(c)∘inv(a)∘m∘a∘c - r_fpp) < tol

    # 2D all pseudo-harmonic oscillators order 6
    m = read_fpp_map("order6var4/test.map",spin=false)
    r_fpp = read_fpp_map("order6var4/R.map",spin=false)
    c = c_map(m)
    a = normal(m)
    @test norm(inv(c)∘inv(a)∘m∘a∘c - r_fpp) < 2e-8

    # 3D coasting last plane order 3
    m = read_fpp_map("coast/test.map",spin=false,coast=true)
    r_fpp = read_fpp_map("coast/R.map",spin=false,coast=true)
    c = c_map(m)
    a = normal(m)
    @test norm(inv(c)∘inv(a)∘m∘a∘c - r_fpp) < tol

    # Equilibrium moments -----------------------------
    Σ_fpp = [ 0.3799021729309453E-03   0.6751391708793735E-04   0.5210559309521585E-04   0.5356022290151431E-04   0.7855831912608747E-04   0.2787256761813361E-04;
              0.6751391708793728E-04   0.1198436499688306E-03   0.2355446920731669E-05   0.2584036833779584E-04   0.4495800939032602E-04  -0.5279019399759065E-05;
              0.5210559309521584E-04   0.2355446920731667E-05   0.9446395990088929E-04   0.6887819755086758E-04   0.1212047267833650E-03   0.1627981838787403E-03;
              0.5356022290151429E-04   0.2584036833779585E-04   0.6887819755086762E-04   0.2881897937276110E-03   0.8460461412432578E-04  -0.8485345202908079E-05;
              0.7855831912608744E-04   0.4495800939032603E-04   0.1212047267833649E-03   0.8460461412432579E-04   0.1116532341434929E-03   0.9345624000572918E-04;
              0.2787256761813363E-04  -0.5279019399759087E-05   0.1627981838787403E-03  -0.8485345202908038E-05   0.9345624000572919E-04   0.2503415143184575E-03]
    m = read_fpp_map("radiation/test.map",spin=false)
    r_fpp = read_fpp_map("radiation/R.map",spin=false)
    c = c_map(m)
    a = normal(m)
    @test norm(inv(c)*inv(a)*m*a*c - r_fpp) < 1e-9
    Σ = equilibrium_moments(m,a)
    @test norm(Σ - Σ_fpp) < tol

    # Factorize -----------------------------
    # 3D all pseudo harmonic oscillators
    a = read_fpp_map("factorize1/a.map")
    a0_fpp = read_fpp_map("factorize1/a0.map")
    a1_fpp = read_fpp_map("factorize1/a1.map")
    a2_fpp = read_fpp_map("factorize1/a2.map")
    
    as, a0, a1, a2 = factorize(a)
    @test norm(a0-a0_fpp) < tol
    @test norm(a1-a1_fpp) < tol
    @test norm(a2-a2_fpp) < tol

    # 3D coasting beam
    a = read_fpp_map("factorize2/a.map",coast=true,spin=true)
    as_fpp = read_fpp_map("factorize2/as.map",coast=true,spin=true)
    a0_fpp = read_fpp_map("factorize2/a0.map",coast=true,spin=true)
    a1_fpp = read_fpp_map("factorize2/a1.map",coast=true,spin=true)
    a2_fpp = read_fpp_map("factorize2/a2.map",coast=true,spin=true)
    
    as, a0, a1, a2 = factorize(a)
    @test norm(a0-a0_fpp) < tol
    @test norm(a1-a1_fpp) < tol
    @test norm(a2-a2_fpp) < tol
    @test norm(as-as_fpp) < tol

    # Spin resonance  (Q_y - Q_s resonance)
    # WHEN LEAVING RESONANCES IN THE MAP, the PHASE of the normal 
    # form will affect things with the spin. Therefore we need to 
    # make sure a and a_fpp are in the same phase
    m = read_fpp_map("spin1/test.map")
    a_fpp = include("spin1/a.jl")
    r_fpp = read_fpp_map("spin1/R.map")
    c = c_map(m)
    a = normal(m,res=[0; 1], spin_res=[-1])
    @test norm(TI.norm_tps.((inv(c)*inv(a)*m*a*c).v .- r_fpp.v)) < 4e-9

    R_rot = jacobian(a, NonlinearNormalForm.HVARS)\jacobian(a_fpp, NonlinearNormalForm.HVARS)    
    r_rot = one(a_fpp)
    NonlinearNormalForm.setray!(r_rot.v, v_matrix=R_rot)
    a = a ∘ r_rot
    @test norm(a - a_fpp) < 2e-7

end

