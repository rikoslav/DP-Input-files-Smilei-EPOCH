import math
import datetime

c = 299792458
electron_mass = 9.10938356e-31
electron_charge = 1.60217662e-19
lambdar = 1e-6                         
wr = 2*math.pi*c/lambdar               
l0 = 2.0*math.pi                       
t0 = l0                                 
Lx = 30*l0                             
n0 = 1e-5                               
Tsim = 50.*t0                           
resx = 128.                             
dx = l0/resx                          
dt  = 0.95 * dx                        
start = 0                            
fwhm = 10*t0                          
duration = 50*t0                      
center = duration*0.5                  
order = 4                              
gamma = 4000./0.511                     
v = math.sqrt(1 - 1./gamma**2)          
pusher = "boris" 
datetime = datetime.datetime.now()
random_seed = datetime.microsecond

# Density profile for inital location of the particles
def n0_(x):
        if (Lx - 10*dx < x < Lx - dx):
                return n0
        else:
                return 0.

def n0_positron(x):
                return 0.

def n0_photon(x):
                return 0.

Main(
    geometry = "1Dcartesian",
    interpolation_order = 4 ,
    cell_length = [dx],
    grid_length  = [Lx],
    number_of_patches = [16],
    timestep = dt,
    simulation_time = Tsim,
    EM_boundary_conditions = [['silver-muller']],
    reference_angular_frequency_SI = wr,
    random_seed = random_seed
)

LaserPlanar1D(
    box_side         = "xmin",
    a0              = 100.,
    omega           = 1.,
    polarization_phi = 0.,
    ellipticity     = 1,
    time_envelope  = tgaussian(start=start,duration=duration,
                               fwhm=fwhm,
                               center=center,
                               order=order)
)

Species(
    name = "electron_MC",
    position_initialization = "random", 
    momentum_initialization = "cold",
    particles_per_cell = 64,
    mass = 1.0,
    charge = -1.0,
    number_density = n0_,
    mean_velocity = [-v, 0.0, 0.0],
    temperature = [0.],
    pusher = pusher,
    radiation_model = "Monte-Carlo",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 2,
    time_frozen = 29*t0,
    boundary_conditions = [
        ["remove", "remove"],
    ],
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
    radiation_photon_gamma_threshold = 2,
    time_frozen = 29*t0,
    boundary_conditions = [["remove"]],
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
    multiphoton_Breit_Wheeler = ["electron_MC","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = [["remove"]],
)

RadiationReaction(
    minimum_chi_discontinuous = 1e-5,
)

MultiphotonBreitWheeler(
)

# Scalar Diagnostics
DiagScalar(
    every = 10
)

    # Weight spatial-distribution
DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["electron_MC"],
    axes = [
        ["x", 0., Lx, 1000],
    ]
)


    # Weight x chi spatial-distribution
DiagParticleBinning(
    deposited_quantity = "weight_chi",
    every = 500,
    time_average = 1,
    species = ["electron_MC"],
    axes = [
        ["x", 0., Lx, 1000],
    ]
)

   # Chi-distribution
DiagParticleBinning(
    deposited_quantity = "weight",
    every = [5000,6500,100],
    time_average = 1,
    species = ["electron_MC"],
    axes = [
        ["chi", 1e-3, 1., 256,"logscale"],
    ]
)

    # Gamma-distribution
DiagParticleBinning(
    deposited_quantity = "weight",
    every = [5000,6500,100],
    time_average = 1,
    species = ["electron_MC"],
    axes = [
        ["gamma", 1., 1.1*gamma, 256,"logscale"],
    ]
)