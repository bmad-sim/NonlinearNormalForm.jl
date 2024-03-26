function show(io::IO, m::Union{Probe,TaylorMap})
  println(io, typeof(m),":")
  lines_used=Ref{Int}(2)

  if GTPSA.show_header
    println(io, "-----------------------")
    desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(first(m.x).tpsa).d))
    lines_used[] += 2 + GTPSA.show_GTPSA_info(io, desc)
    println(io, "-----------------------")
  end
  println(io, "Reference Orbit ", typeof(m.x0),":")
  for i =1:length(m.x0)
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    @printf(io, "%-3s  ", "$(i):"); println(io, m.x0[i])
    lines_used[] += 1
  end
  !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
  println(io)
  lines_used[] += 1
  !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
  println(io, "Orbital Ray ", typeof(m.x),":")
  lines_used[] += 1

  if eltype(m.x) <: Union{TPS,ComplexTPS}
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    def = true
    for i in eachindex(m.x)
      if !isassigned(m.x, i)
        def = false
      end
    end
    if def
      GTPSA.show_map!(io, m.x, lines_used, true, 1:numvars(m))
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || return
    else
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      println(io)
      lines_used[] += 1
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      println(io, "\tAtleast one $(eltype(m.x)) is undefined!")
      lines_used[]+=1
    end
  else
    for i =1:length(m.x0)
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      @printf(io, "%-3s  ", "$(i):"); println(io, m.x[i])
      lines_used[] += 1
    end
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 ||  (println(io, "\t⋮"); return)
  end

  if typeof(m.Q) != Nothing
    println(io)
    lines_used[]+= 1
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    println(io,typeof(m.Q),":")
    lines_used[] += 1

    
    if eltype(m.Q.q) <: Union{TPS,ComplexTPS}
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      def = true
      for i in eachindex(m.Q.q)
        if !isassigned(m.Q.q, i)
          def = false
        end
      end
      if def
        GTPSA.show_map!(io, collect(m.Q.q), lines_used, false, [" q0:"," q1:"," q2:"," q3:"])
      else
        !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
        println(io)
        lines_used[] += 1
        !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
        println(io, "\tAtleast one $(eltype(m.Q.q)) is undefined!")
        lines_used[]+=1
      end
    else
      for i=1:4
        !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
        @printf(io, "%-3s  ", "$(i):"); println(io, m.Q.q[i])
        lines_used[] += 1
      end
    end
  end
  
end

#=
function show(io::IO, Q::Quaternion)

end=#