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
gamma = 1000./0.511                    
v = math.sqrt(1 - 1./gamma**2)         
pusher = "boris" 
radiation_list = ["Landau-Lifshitz","corrected-Landau-Lifshitz",
                    "Niel","Monte-Carlo"] 
species_name_list = ["LL","CLL","Niel","MC"]            
datetime = datetime.datetime.now()
random_seed = datetime.microsecond

# Density profile for inital location of the particles
def n0_(x):
        if (Lx - 10*dx < x < Lx - dx):
                return n0
        else:
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

# Laser definition
LaserPlanar1D(
    box_side         = "xmin",
    a0              = 100., # = 1., 10., 270.,
    omega           = 1.,
    polarization_phi = 0.,
    ellipticity     = 1,
    time_envelope  = tgaussian(start=start,duration=duration,
                               fwhm=fwhm,
                               center=center,
                               order=order)
)

# Loop to create all the species
for i,radiation in enumerate(radiation_list):
    Species(
        name = "electron_" + species_name_list[i],
        position_initialization = "random", 
        momentum_initialization = "cold",
        particles_per_cell = 64,
        mass = 1.0,
        charge = -1.0,
        number_density = n0_,
        mean_velocity = [-v, 0.0, 0.0],
        temperature = [0.],
        pusher = pusher,
        radiation_model = radiation,
        time_frozen = 29*t0,
        boundary_conditions = [
            ["remove", "remove"],
        ],
    )

# Radiation parameters
RadiationReaction(
    minimum_chi_discontinuous = 1e-5,
)

# Scalar Diagnostics
DiagScalar(
    every = 10
)

# Loop to create all the species particle binning diagnostics
# One species per radiation implementations
for i,radiation in enumerate(radiation_list):

    # Weight spatial-distribution
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = 500,
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["x", 0., Lx, 1000],
        ]
    )


for i,radiation in enumerate(radiation_list):
    # Weight x chi spatial-distribution
    DiagParticleBinning(
        deposited_quantity = "weight_chi",
        every = 500,
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["x", 0., Lx, 1000],
        ]
    )


for i,radiation in enumerate(radiation_list):
    # Chi-distribution
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = [5000,6500,100],
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["chi", 1e-3, 1., 256,"logscale"],
        ]
    )

for i,radiation in enumerate(radiation_list):
    # Gamma-distribution
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = [5000,6500,100],
        time_average = 1,
        species = ["electron_" + species_name_list[i]],
        axes = [
            ["gamma", 1., 1.1*gamma, 256,"logscale"],
        ]
    )
