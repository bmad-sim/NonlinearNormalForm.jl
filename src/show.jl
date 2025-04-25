function show(io::IO, m::TaylorMap)
  println(io, typeof(m),":")
  #lines_used=Ref{Int}(2)
  println(io, "Reference Orbit ", typeof(m.v0),":")
  println(io, m.v0)
  #lines_used += 1m 
  #=
  for i in 1:length(m.v0)
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    @printf(io, "%-3s  ", "$(i):"); println(io, m.v0[i])
    lines_used[] += 1
  end
  =#
  #!get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
  println(io)
  #lines_used[] += 1
  ndpt = coastidx(m)
  if ndpt != -1
    #!get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    println("Last plane is coasting: variable #", ndpt, " is constant")
    #lines_used[] += 1
  end
  #!get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
  println(io, "Orbital Ray ", typeof(m.v),":")
  #lines_used[] += 1
  println(io, m.v)
#=
  if eltype(m.v) <: TPS
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    def = true
    for i in eachindex(m.v)
      if !isassigned(m.v, i)
        def = false
      end
    end
    if def
      GTPSA.show_map!(io, m.v, lines_used, true, 1:numvars(m))
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || return
    else
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      println(io)
      lines_used[] += 1
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      println(io, "\tAtleast one $(eltype(m.v)) is undefined!")
      lines_used[]+=1
    end
  else
    for i in 1:length(m.v0)
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      @printf(io, "%-3s  ", "$(i):"); println(io, m.v[i])
      lines_used[] += 1
    end
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 ||  (println(io, "\t⋮"); return)
  end
=#

  if !isnothing(m.q)
    println(io)
    #lines_used[]+= 1
    #!get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    println(io,typeof(m.q),":")
    #lines_used[] += 1

    println(io, m.q)
    #!get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    #GTPSA.show_map!(io, collect(m.q), lines_used, false, [" q0:"," q1:"," q2:"," q3:"])
  end
  
end

#=
function show(io::IO, q::Quaternion{<:TPS})
  println(io, "$(typeof(q)):")
  GTPSA.show_map!(io, collect(q), Ref{Int}(1), false, [" q0:"," q1:"," q2:"," q3:"])
end
=#
#show(io::IO, ::MIME"text/plain", q::Quaternion{<:TPS}) = (println(io, "$(typeof(q)):"); GTPSA.show_map!(io, collect(q), Ref{Int}(1), false, [" q0:"," q1:"," q2:"," q3:"]))

