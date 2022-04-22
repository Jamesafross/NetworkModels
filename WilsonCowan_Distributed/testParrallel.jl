@everywhere using LinearAlgebra,Distributed,SharedArrays,Parameters,DifferentialEquations

if nprocs() < 2
addprocs(2)
end

@everywhere @with_kw struct params{R}
    p1::R = 0.1
    p2::R = 0.2
end

Pa = params()

@everywhere function addParams(P)
    @unpack p1,p2 = P
    return p1+p2
end

a = SharedArray(zeros(5,5,10))

B = [1 2 3 4 5 6 7 8 9 10]

@sync @distributed for i in 1:10

  println("working on i=$i")
  for j = 1:10
  a[:,:,i] += 0.1*randn(5,5)
  end
end
