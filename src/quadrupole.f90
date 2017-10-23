program main

  use PolyMapMod
  use TPSAMod
  Implicit None

  type (dctps)   :: x,px,y,py,m
  
  call dctps_Initialize(4,6)
  
  call assign(x,  0.0d0, 1)
  call assign(px, 0.0d0, 2)
  call assign(y,  0.0d0, 3)
  call assign(py, 0.0d0, 4)

  m = px
  
  call m%output()
  
end program main
