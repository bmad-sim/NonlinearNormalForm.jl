using NonlinearNormalForm
using Test

@testset "NonlinearNormalForm.jl" begin
    d = Descriptor(1,2)
    x1 = vars()[1]
    m1 = DAMap([1+2*x1+2*x1^2], x0=[4])
    m2 = DAMap([1+2*x1+2*x1^2], x0=[3])

    mt1 = TPSAMap(m1)
    mt2 = TPSAMap(m2)

    @test m2∘m1 == DAMap([1+4*x1+12*x1^2], x0=[4])
    @test mt2∘mt1 == TPSAMap([5-12*x1-4*x1^2], x0=[4])
end
