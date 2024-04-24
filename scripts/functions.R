##### FUNCTIONS #####

#### Particle Characteristics Equations ####
#surface area equation for elongated spheres (fragments)
SAfnx = function(a, # length
                 b, # width
                 c){ # height
  SA = 4*pi*(((a*b)^1.6 + (a*c)^1.6 + (b*c)^1.6) / 3)^(1/1.6)
  return(SA)}

#cylinder equation for SA. S = 2pi*r*h + 2pi*r^2, where r = width/2 and h = length 
SAfnx_fiber = function(width, length){
  radius = width / 2
  SA = 2*pi*radius*length + 2*pi*radius^2
  return(SA)
}

# equation for volume
volumefnx_poly = function(width, length){
  height = width #0.67 * width
  volume = (4/3) * pi * (length/2) * (width/2) * (height/2) 
  return(volume)}

#Volume equation for elongated sphere (fragments)
volumefnx = function(R, L){
  volume = 0.111667 * pi * R^2 * L^3 #assumes height = 0.67 * Width, and Width:Length ratio is 'R' (0.77 average in marine surface water)
  return(volume)}

#equation for fibers (cylinder) V = pi*r^2*h (where r = particle width/2 and h = particle length). Assume 15 um if width not reported (kooi et al 2021)
volumefnx_fiber = function(width, length){
  radius = width/2
  volume = pi * (radius) ^ 2 * length
  return(volume)
}

massfnx_poly = function(width, length, p){
  height = width #0.67 * width
  volume = (4/3) * pi * (length/2) * (width/2) * (height/2)  
  mass = p * #density (g/cm^3)
    volume * # volume (um^3): assumes height = 0.67 * Width, and Width:Length ratio is 'R' (compartment-specific)
    1/1e12 * 1e6 #correction factor
  return(mass)}

#### Ecologically Relevant Metric Functions (used in reactives with user-input params) ####

###function to derive correction factor (CF) from Koelmans et al (equation 2)
CFfnx = function(a, #default alpha from Koelmans et al (2020)
                 x2D, #set detault values to convert ranges to (1-5,000 um) #5mm is upper defuault 
                 x1D, #1 um is lower default size
                 x2M, x1M){
  CF = (x2D^(1-a)-x1D^(1-a))/(x2M^(1-a)-x1M^(1-a)) 
  return(CF)}

#### equations for mu_x_poly (note that there are three depending on certain alphas for limits of equation)
##### if alpha does not equal 2 #####
mux.polyfnx = function(a.x, 
                       x_UL, 
                       x_LL){
  mux.poly = ((1-a.x)/(2-a.x)) * ((x_UL^(2-a.x) - x_LL^(2-a.x))/(x_UL^(1-a.x) - x_LL^(1-a.x)))
  return(mux.poly)}

##### If alpha does equal 2 #####
mux.polyfnx.2 = function(x_UL,x_LL){
  mux.poly = (log(x_UL/x_LL))/(x_LL^(-1) - x_UL^-1)
  return(mux.poly)}

#max ingestible specific surface area
SSA.inversefnx = function(sa, #surface area, calcaulted elsewhere
                          m){ #mass, calculated elsewhere
  SSA.inverse = m / sa
  return(SSA.inverse)}

#data tidying functions from Ana

############## Levels summary ##################
summarize_and_print <- function(data, column_name)
{
  result <- data %>%
    group_by({{ column_name }}) %>%
    summarise(n_datapoints = n()) %>%
    arrange(as.numeric(as.character({{ column_name }}))) %>%
    print(n = 1000)
  return(result)
}


############## Change particle length ##################
update_particle_length <- function(data, doi, length, polymer, shape, new_value) {
  data$Particle.Length..μm.[
    data$DOI == doi &
      data$Particle.Length..μm. == length &
      data$Polymer == polymer &
      data$Shape == shape
  ] = new_value
  
  return(data)
}


###### check what is missing #########
generate_structure_checks <- function(data) {
  structure.checks <- data.frame(
    na.counts = sapply(data, function(x) sum(is.na(x))),
    na.percent = round(sapply(data, function(x) sum(is.na(x)) / nrow(data) * 100), digits = 1),
    n.levels = sapply(data, function(x) length(unique(x)))
  )
  return(structure.checks)
  
}




