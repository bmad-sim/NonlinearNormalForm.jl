using NonlinearNormalForm
using Test

@testset "NonlinearNormalForm.jl" begin
    d = Descriptor(1,2)
    x1 = vars()[1]
    m1 = DAMap(x=[1+2*x1+2*x1^2], x0=[4.])
    m2 = DAMap(x=[1+2*x1+2*x1^2], x0=[3.])

    mt1 = TPSAMap(m1)
    mt2 = TPSAMap(m2)
    
    tol = 1e-10

    @test norm(m2∘m1 - DAMap(x=[1+4*x1+12*x1^2], x0=[4.])) < tol
    @test norm(mt2∘mt1 - TPSAMap(x=[5-12*x1-4*x1^2], x0=[4.])) < tol
    @test norm(m2^3 - m2∘m2∘m2) < tol
    @test norm(mt2^3 - mt2∘mt2∘mt2) < tol
    @test norm(m2^3 - DAMap(x=[1+8*x1+56*x1^2], x0=[3.])) < tol
    @test norm(mt2^3 - TPSAMap(x=[13+8*x1-8*x1^2], x0=[3.])) < tol

    # with temporary inverter"
    #@test norm(m1^-3 - DAMap(x=[4+0.125*x1-0.109375*x1^2],x0=[1])) < tol
    #@test norm(m1^-3 - inv(m1^3)) < tol
    #@test norm(mt1^-3 - inv(mt1^3)) < tol
    #@test norm(mt1^-3 - TPSAMap(x=[4+0.9615384615384616E-02*x1-0.4978379608557124E-04*x1^2],x0=[-35])) < tol 
end