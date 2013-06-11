program real
	use init                ! init_model
	use boundary_conditions ! bc_*    routines, constant in an ideal simulation? 
	use data_structures     ! *_type  types
	use microphysics        ! mp_init
	use output              ! write_domain
	use time_step			! step
	
	implicit none
	type(options_type) :: options
	type(domain_type)  :: domain
	type(bc_type)      :: boundary
	integer::i
	
	call init_model("real_options.namelist",options,domain,boundary)
! 	initialize microphysics code (e.g. compute look up tables in Thompson et al)
	call mp_init(options%physics%microphysics) !this could easily be moved to init_model...
! 	read initial conditions from the boundary file
	call bc_init(domain,boundary,options)
	
! 	note that a timestep here is an IO timestep O(1hr), not a physics timestep O(20s)
	do i=1,options%ntimesteps
! 		update boundary conditions (dXdt variables)
		call bc_update(domain,boundary,options)
! 		this is the meat of the model, run all the physics for the current time step looping over internal timesteps
		call step(domain,options,boundary)
! 		finally write the output for this timestep
		call write_domain(domain,options,i)
	end do
	
end program real
