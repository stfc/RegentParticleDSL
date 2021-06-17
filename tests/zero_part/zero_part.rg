import "regent"

fspace part{
  a : int32,
  b : uint32,
  c: int64,
  d : uint64,
  e : bool,
  f : float,
  g : double,
  h : int1d,
  i : int2d, 
  j : int3d,
  k : int64[24],
  l : bool[3],
  m : int3d[2][26]
}
require("src/particles/init_part")

local zero_part_task = generate_zero_part_func()


task main()

  var parts = region(ispace(int1d, 1), part);
  zero_part_task(parts)
  regentlib.assert(parts[0].a == 0, "Failed to allocate a")
  regentlib.assert(parts[0].b == 0, "Failed to allocate b")
  regentlib.assert(parts[0].c == 0, "Failed to allocate c")
  regentlib.assert(parts[0].d == 0, "Failed to allocate d")
  regentlib.assert(parts[0].e == false, "Failed to allocate e")
  regentlib.assert(parts[0].f == 0.0, "Failed to allocate f")
  regentlib.assert(parts[0].g == 0.0, "Failed to allocate g")
  regentlib.assert(parts[0].h == int1d(0), "Failed to allocate h")
  regentlib.assert(parts[0].i == int2d({0,0}), "Failed to allocate i")
  regentlib.assert(parts[0].j == int3d({0,0,0}), "Failed to allocate j")
  for i = 0, 24 do
     regentlib.assert(parts[0].k[i] == 0, "Failed to zero k")
  end
  for i = 0,3 do
    regentlib.assert(parts[0].l[i] == false, "Failed to zero l")
  end
  for i = 0, 2 do
    for j = 0, 26 do
        regentlib.assert(parts[0].m[j][i] == int3d({0,0,0}), "Failed to zero m")
    end
  end
end

regentlib.start(main)


