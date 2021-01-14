import "regent"

--h = 1.3 * dx
--kernelchoice=7
--h2 = h*h
--h3 = h*h*h
--sup_size_2   = 9.0d0*h2
--sup_size     = 3.0d0*h
--uno_sup_size = 1.0d0/(3.0d0*h)
--ad_7=3.d0/(359.d0*pi*h3)
--ad_7h=ad_7/h

local sqrt = regentlib.sqrt(double)

terra fac(qq : double) : double

  var factemp : double
  var q3 : double
  var q2 : double
  var q1 : double
  var q34 : double
  var q24 : double
  var q14 : double
  
  q3 = 3.0 - qq
  q2 = 2.0 - qq
  q1 = 1.0 - qq
  
  q34 = q3 * q3
  q34 = q34 * q34
  
  q24 = q2 * q2
  q24 = q24 * q24
  
  q14 = q1 * q1
  q14 = q14 * q14
  
  if (qq > 0.0 and qq < 1.0) then
    factemp = ad_7h * (-5.0*q34 + 30.0*q24 - 75.0*q14)
  elseif (qq >= 1.0 and qq < 2.0) then
    factemp = ad_7h * (-5.0 * q34 + 30.0 * q24)
  elseif (qq >= 2.0 and qq < 3.0) then
    factemp = ad_7h * (-5.0 * q34)
  elseif (qq >= 3.0) then
    factemp = 0.0
  end

factemp = factemp / (qq*h)

return factemp
end

--Asymmetric kernel
function divergence(part1, part2, r2)

local kernel = rquote
  var rad = sqrt(r2)
  var qq = rad / h
  var facs = fac(qq)
  
  var drx = part1.core_part_space.pos_x - part2.core_part_space.pos_x
  var dry = part1.core_part_space.pos_y - part2.core_part_space.pos_y
  var drz = part1.core_part_space.pos_z - part2.core_part_space.pos_z

  var frx = facs * drx
  var fry = facs * dry
  var frz = facs * drz

  var temp = drx*frx + dry*fry + drz*frz
  part1.divergence = part1.divergence - part2.volume*temp
  part1.interactions = part1.interactions + 1
end
return kernel
end
