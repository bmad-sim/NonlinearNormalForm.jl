function show(io::IO, m::Union{Probe,TaylorMap})
  println(io, typeof(m),":")
  lines_used=Ref{Int}(2)
  if eltype(m.x) <: TPS
    diffdescsx = false
    for i in eachindex(m.x)
      if !diffdescsx && getdesc(first(m.x)) != getdesc(m.x[i])
        println(io, "WARNING: Atleast one $(eltype(m.x)) in the orbital ray has a different Descriptor!")
        diffdescsx = true
        lines_used[] += 1
      end
    end
    if eltype(m.Q) <: TPS
      diffdescsq = false
      for qi in m.Q.q
        if !diffdescsq && getdesc(first(m.Q.q)) != getdesc(qi)
          println(io, "WARNING: Atleast one $(eltype(m.Q.q.q)) in the quaternion has a different Descriptor!")
          diffdescsq = true
          lines_used[] += 1
        end
      end
      diffdescsxq = false
      if getdesc(first(m.x)) != getdesc(first(m.Q.q))
        println(io, "WARNING: First element in orbital ray has different Descriptor than first element in quaternion!")
        diffdescsxq = true
        lines_used[] += 1
      end
    else
      diffdescsq = false
      diffdescsxq = false
    end
    if GTPSA.show_header
      if diffdescsq || diffdescsx || diffdescsxq
        println(io, "Cannot show GTPSA header (non-unique Descriptor).")
        lines_used[] += 1
      else
        println(io, "-----------------------")
        desc = getdesc(first(m.x))
        lines_used[] += 2 + GTPSA.show_GTPSA_info(io, desc)
        println(io, "-----------------------")
      end
    end
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
  if !isnothing(m.idpt)
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    println("Last plane is coasting: variable #", numvars(m)-1+m.idpt, " is constant")
    lines_used[] += 1
  end
  !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
  println(io, "Orbital Ray ", typeof(m.x),":")
  lines_used[] += 1

  if eltype(m.x) <: TPS
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

  if !isnothing(m.Q)
    println(io)
    lines_used[]+= 1
    !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
    println(io,typeof(m.Q),":")
    lines_used[] += 1

    
    if eltype(m.Q) <: TPS
      !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
      GTPSA.show_map!(io, collect(m.Q.q), lines_used, false, [" q0:"," q1:"," q2:"," q3:"])
    else
      i=1
      for qi in m.Q.q
        !get(io, :limit, false) || lines_used[] < displaysize(io)[1]-5 || (println(io, "\t⋮"); return)
        @printf(io, "%-3s  ", "q$(i):"); println(io, qi)
        lines_used[] += 1
        i+=1
      end
    end
  end
  
end
