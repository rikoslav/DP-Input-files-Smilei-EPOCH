from math import pi, sqrt

c = 299792458
lambdar = 1e-6
wr = 2*pi*c/lambdar
l0 = 2.0*pi             
t0 = l0                      
Lx = 30*l0                   
Ly = 4*l0                     
n0 = 1e-5                   
Tsim = 60.*t0               
resx = 64.                   
resy = 64.
dx = l0/resx                           
dy = l0/resy
dt  = 0.95 * 1./np.sqrt(1./dx**2 + 1./dy**2)         
start = 0                               
fwhm = 10*t0                            
duration = 50*t0                       
center = duration*0.5                  
order = 4                             
gamma = 4000./0.511                    
v = sqrt(1 - 1./gamma**2)         

def n0_electron(x,y):
        if ((0.97*Lx < x < 0.99*Lx)and(0.48*Ly < y < 0.52*Ly)):
                return n0
        else:
                return 0.

def n0_positron(x,y):
                return 0.

def n0_photon(x,y):
                return 0.

Main(
    geometry = "2Dcartesian",
    interpolation_order = 4 ,
    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],
    number_of_patches = [4,4],
    timestep = dt,
    simulation_time = Tsim,
    EM_boundary_conditions = [
        ["silver-muller", "silver-muller"],
        ["periodic", "periodic"]
    ],
    reference_angular_frequency_SI = wr,
    random_seed = 0
)

LaserGaussian2D(
    box_side         = "xmin",
    a0              = 100.,
    omega           = 1.,
    focus           = [0.5*Lx, 0.5*Ly],
    waist           = 1E9,
    incidence_angle = 0.,
    polarization_phi = 0.,
    ellipticity     = 0,
    time_envelope  = tgaussian(start=start,duration=duration,
                               fwhm=fwhm,
                               center=center,
                               order=order)
)

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 32,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = [-v, 0.0, 0.0],
    temperature = [0.],
    pusher = "boris", #zmena z vay
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 0,
    time_frozen = 29*t0,
    boundary_conditions = [["remove","remove"],
                            ["periodic","periodic"]]
)

Species(
    name = "positron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    mass = 1.0,
    charge = 1.0,
    charge_density = n0_positron,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.],
    pusher = "boris",
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 0,
    time_frozen = 29*t0,
    boundary_conditions = [["remove","remove"],
    ["periodic","periodic"]]
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    mass = 0,
    charge = 0.,
    number_density = n0_photon,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "boris",
    multiphoton_Breit_Wheeler = ["electron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = [["remove","remove"],
    ["periodic","periodic"]]
)

RadiationReaction(
    minimum_chi_continuous = 1e-4,
    minimum_chi_discontinuous = 1e-4,
)

MultiphotonBreitWheeler(
)

DiagScalar(
    every = 10,
)

DiagFields(
    every = 500,
    fields = ['Ex','Ey','Ez','By','Bz']
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 500]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["positron"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 500]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["photon"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 500]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 500]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 500,
    time_average = 1,
    species = ["positron"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 500]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 500,
    time_average = 1,
    species = ["photon"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 500]
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight_chi",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 1000],
    ["y", 0., Ly, 500]
    ]
)

DiagParticleBinning(
     deposited_quantity = "weight_chi",
     every = 500,
     time_average = 1,
     species = ["positron"],
     axes = [
         ["x", 0., Lx, 1000],
     ["y", 0., Ly, 500]
     ]
)


DiagParticleBinning(
     deposited_quantity = "weight_chi",
     every = 500,
     time_average = 1,
     species = ["photon"],
     axes = [
         ["x", 0., Lx, 1000],
     ["y", 0., Ly, 500]
     ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["gamma", 10., 3000, 200,'logscale']
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["positron"],
    axes = [
        ["gamma", 10., 3000, 200,'logscale']
    ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["photon"],
    axes = [
        ["gamma", 2., 3000, 200,'logscale']
    ]
)