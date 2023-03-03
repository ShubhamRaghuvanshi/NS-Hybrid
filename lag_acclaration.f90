	subroutine lag_acclaration 

	use mod_serial_fluid
  implicit none

	Ap = (Vp-Vp0)/delta	
	Vp0 = Vp
		
	end subroutine lag_acclaration 
