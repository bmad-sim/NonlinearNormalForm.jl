function show(io::IO, m::Union{Probe,TaylorMap})
  println(io, typeof(m),":")
  lines_used=Ref{Int}(2)

  if GTPSA.show_header
    println(io, "-----------------------")
    desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(first(m.x).tpsa).d))
    GTPSA.show_GTPSA_info(io, desc)
    println(io, "-----------------------")
    lines_used[] += 2 + (desc.nv > 0 ? 2 : 0) + (desc.np > 0 ? 2 : 0)
  end
  println(io, "Reference Orbit ", typeof(m.x0),":")
  for i in LinearIndices(m.x)
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
    GTPSA.show_map!(io, m.x, lines_used, true)
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || return
  else
    for i in LinearIndices(m.x)
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      @printf(io, "%-3s  ", "$(i):"); println(io, m.x[i])
      lines_used[] += 1
    end
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 ||  (println(io, "\t⋮"); return)
  end

  println(io)
  lines_used[]+= 1
  !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
  println(io,typeof(m.q),":")
  lines_used[] += 1

  
  if eltype(m.q.q) <: Union{TPS,ComplexTPS}
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    GTPSA.show_map!(io, [m.q.q...], lines_used, false, [" q0:"," q1:"," q2:"," q3:"])
  else
    for i=1:4
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      @printf(io, "%-3s  ", "$(i):"); println(io, m.q.q[i])
      lines_used[] += 1
    end
  end
  
end

function show(io::IO, Q::Quaternion)

end