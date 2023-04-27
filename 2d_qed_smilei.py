begin:constant
 t_sim = 60.0 * 1.0e-6 / c
 omega_r = 2.0 * pi * c / 1.0e-6
end:constant

begin:control
 nx = 1920
 ny = 256
 t_end = t_sim
 x_min = 0.0
 x_max = 30.0e-6
 y_min = -2.0e-6
 y_max = 2.0e-6
 dt_multiplier = 0.95
 field_order = 2
 maxwell_solver = yee
 use_random_seed = T
 particle_tstart = (29.0 / 60) * t_sim
 dlb_threshold = 0.6
end:control

begin:boundaries
  bc_x_min = simple_laser
  bc_x_max = simple_laser
  bc_y_max = periodic
  bc_y_min = periodic
end:boundaries

begin:constant
  a0 = 100.0 # 1., 10., 270.,
  w0 = 1000 / (2 * pi)    
  I_peak_Wcm2 = (a0/0.86)^2 * 1.0e18
  las_lambda = 1.0e-6      
  foc_dist = 15.0e-6        
  width = t_sim / (6* 2 * (loge(2))^(1/4))
end:constant

begin:constant
 las_k = 2.0 * pi / las_lambda    
 ray_rang = pi * w0^2 / las_lambda                   
 w_boundary = w0 * sqrt(1.0 + (foc_dist/ray_rang)^2)  
 I_boundary = I_peak_Wcm2 * (w0 / w_boundary)^2       
 rad_curve = foc_dist * (1.0 + (ray_rang/foc_dist)^2) 
 gouy = atan(-foc_dist/rad_curve)                     
end:constant

begin:laser
    boundary = x_min
    intensity_w_cm2 = I_boundary
    lambda = las_lambda
    phase = las_k * y^2 / (2.0 * rad_curve) - gouy
    profile = gauss(y, 0, w_boundary)
    pol_angle = 0
    t_profile = supergauss(time, (50.0/(60.0*2)) * t_sim, width, 4)
    t_start = 0.0
    t_end = (50.0/60.0) * t_sim
end:laser

begin:qed
 use_qed = T
 qed_start_time = 0.0
 produce_photons = T
 photon_energy_min = 0.0
 photon_dynamics = F
 produce_pairs = F
 use_radiation_reaction = T
end:qed

begin:constant
 gamma = 1000.0 / 0.511
 drift = gamma * me * c * sqrt(1 - (1.0/((gamma)^2)) )
 n_0 = 1.0e-5
 density_0 = n_0 * critical(omega_r)
 l_x = 30.0
 x_left = 0.97 * l_x * 1e-6
 x_right = 0.99 * l_x * 1e-6
 l_y = 2.0
 y_abs = (2/50) * l_y * 1e-6
end:constant

begin:species
 name = electron
 npart_per_cell = 32
 temp = 0.0
 drift_x = -drift
 number_density = if( (x gt x_left) and (x lt x_right) 
                    and (abs(y) lt y_abs), density_0, 0.0)
 identify:electron 
end:species

begin:species
 name = photon
 nparticles = 0
 identify:photon
end:species

begin:output
 nstep_snapshot = 10
 poynt_flux = always
 charge_density = always
 grid = always
 ey = always
 ex = always
 ez = always
 particles = always
 particle_energy = always + species
 particle_weight = always + species
 px = always + species
 py = always + species
 pz = always + species
 total_energy_sum = always + species
end:output