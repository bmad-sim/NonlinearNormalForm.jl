using NonlinearNormalForm
using Test

@testset "Composition and inversion" begin
    d = Descriptor(1,2)
    x1 = vars()[1]
    m1 = DAMap(x=[1+2*x1+2*x1^2], x0=[4])
    m2 = DAMap(x=[1+2*x1+2*x1^2], x0=[3])

    mt1 = TPSAMap(m1)
    mt2 = TPSAMap(m2)
    
    tol = 1e-10

    @test norm(m2∘m1 - DAMap(x=[1+4*x1+12*x1^2], x0=[4])) < tol
    @test norm(mt2∘mt1 - TPSAMap(x=[5-12*x1-4*x1^2], x0=[4])) < tol
    @test norm(m2^3 - m2∘m2∘m2) < tol
    @test norm(mt2^3 - mt2∘mt2∘mt2) < tol
    @test norm(m2^3 - DAMap(x=[1+8*x1+56*x1^2], x0=[3])) < tol
    @test norm(mt2^3 - TPSAMap(x=[13+8*x1-8*x1^2], x0=[3])) < tol

    @test norm(m1^-3 - DAMap(x=[4+0.125*x1-0.109375*x1^2],x0=[1])) < tol
    @test norm(m1^-3 - inv(m1^3)) < tol
    @test norm(mt1^-3 - inv(mt1^3)) < tol
    @test norm(mt1^-3 - TPSAMap(x=[4+0.9615384615384616E-02*x1-0.4978379608557124E-04*x1^2],x0=[-35])) < tol 
end

@testset "Comparison with FPP" begin

    tol = 1e-12


    # Normal form -------------------------------------
    # 2D all pseudo-harmonic oscillators order 3
    m = read_fpp_map("order2/test.map",spin=false)
    R_fpp = read_fpp_map("order2/R.map",spin=false)
    c = to_phasor(m)
    a = normal(m).a
    @test norm(inv(c)∘inv(a)∘m∘a∘c - R_fpp) < tol

    # 2D all pseudo-harmonic oscillators order 10
    m = read_fpp_map("order10/test.map",spin=false)
    R_fpp = read_fpp_map("order10/R.map",spin=false)
    c = to_phasor(m)
    a = normal(m).a
    @test norm(inv(c)∘inv(a)∘m∘a∘c - R_fpp) < tol

    # 4D all pseudo-harmonic oscillators order 6
    m = read_fpp_map("order6var4/test.map",spin=false)
    R_fpp = read_fpp_map("order6var4/R.map",spin=false)
    c = to_phasor(m)
    a = normal(m).a
    @test norm(inv(c)∘inv(a)∘m∘a∘c - R_fpp) < 1e-8

    # 6D coasting last plane order 3
    m = read_fpp_map("coast/test.map",idpt=true,spin=false)
    R_fpp = read_fpp_map("coast/R.map",idpt=true,spin=false)
    c = to_phasor(m)
    a = normal(m).a
    @test norm(inv(c)∘inv(a)∘m∘a∘c - R_fpp) < tol

    # Equilibrium moments -----------------------------
    Σ_fpp = [ 0.3799021729309453E-03   0.6751391708793735E-04   0.5210559309521585E-04   0.5356022290151431E-04   0.7855831912608747E-04   0.2787256761813361E-04;
         0.6751391708793728E-04   0.1198436499688306E-03   0.2355446920731669E-05   0.2584036833779584E-04   0.4495800939032602E-04  -0.5279019399759065E-05;
         0.5210559309521584E-04   0.2355446920731667E-05   0.9446395990088929E-04   0.6887819755086758E-04   0.1212047267833650E-03   0.1627981838787403E-03;
         0.5356022290151429E-04   0.2584036833779585E-04   0.6887819755086762E-04   0.2881897937276110E-03   0.8460461412432578E-04  -0.8485345202908079E-05;
         0.7855831912608744E-04   0.4495800939032603E-04   0.1212047267833649E-03   0.8460461412432579E-04   0.1116532341434929E-03   0.9345624000572918E-04;
         0.2787256761813363E-04  -0.5279019399759087E-05   0.1627981838787403E-03  -0.8485345202908038E-05   0.9345624000572919E-04   0.2503415143184575E-03]
    m = read_fpp_map("radiation/test.map",spin=false)
    R_fpp = read_fpp_map("radiation/R.map",spin=false)
    c = to_phasor(m)
    a = normal(m).a
    @test norm(inv(c)*inv(a)*m*a*c - R_fpp) < tol
    Σ = equilibrium_moments(m,a)
    @test norm(Σ - Σ_fpp) < tol



end

