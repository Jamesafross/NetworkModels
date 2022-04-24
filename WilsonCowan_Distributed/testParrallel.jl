 using LinearAlgebra,Distributed,SharedArrays,Parameters,DifferentialEquations

mutable struct testP
    b::Matrix{Float64}
end

b = [1.0 2.0 ; 1.0 4.0]
c = zeros(2,2)
c .= b


function test(b)
  tP = testP(b.+2.0)
  tP.b .= tP.b .+4.0

  return tP.b
end

b[:,:] = test(b)