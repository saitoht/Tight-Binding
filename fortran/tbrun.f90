!!! run the Tight-Binding code !!!
program tbrun
  use m_read_prms
  use m_TB
  implicit none
  
  ! get_parameters
  call get_prms()
  ! band calculation
  call init_band()
end program tbrun
