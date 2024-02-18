function show(io::IO, m::TaylorMap)
  println(io, typeof(m),":")
  extralines=0
  if GTPSA.show_header
    println(io, "-----------------------")
    desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(first(m.v).tpsa).d))
    show_GTPSA_info(io, desc)
    println(io, "-----------------------")
    extralines = 2 + (desc.nv > 0 ? 2 : 0) + (desc.np > 0 ? 2 : 0)
  end
  println(io, "x0 = ", m.x0)
  GTPSA.show_map(io, m.v, 1+extralines, true)
end

function show(io::IO, m::Probe{TPS})
  println(io, typeof(m),":")
  extralines=0
  if GTPSA.show_header
    println(io, "-----------------------")
    desc = unsafe_load(Base.unsafe_convert(Ptr{Desc}, unsafe_load(first(m.v).tpsa).d))
    show_GTPSA_info(io, desc)
    println(io, "-----------------------")
    extralines = 2 + (desc.nv > 0 ? 2 : 0) + (desc.np > 0 ? 2 : 0)
  end
  #println(io, "x0 = ", m.x0)
  GTPSA.show_map(io, m.x0 + m.v, extralines, true)
end