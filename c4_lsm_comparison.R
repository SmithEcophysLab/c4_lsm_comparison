# c4_lsm_comparison.R
## script to compare c4 leaf photosynthesis responses in lsms
## all models use some version of collatz
## steps:
### 1. temperature response function for each lsm
### 2. collatz model function
### 3. lsm-specific simulations in response to light, temperature, CO2
### 4. make figures

### 0. pressure function
calc_patm = function(z) {
  
  kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
  kTo = 298.15   # base temperature, K (Prentice, unpublished)
  kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
  kG = 9.80665   # gravitational acceleration, m/s**2 (Allen, 1973)
  kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
  kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  
  patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
  
  patm
  
}

### 1. temperature response function for each lsm

ftemp_clm <- function(tleaf = 298.15){
  
  q10_func = 2 ^ (0.1 * (tleaf - 298.15))
  temp_lim_1 = 1 + exp(0.3 * (tleaf - 313.15))
  temp_lim_2 = 1 + exp(0.2 * (288.15 - tleaf))
  
  ftemp = q10_func / (temp_lim_1 * temp_lim_2)
  
  ftemp
  
}

ftemp_jules <- function(tleaf = 298.15){
  
  q10_func = 2 ^ (0.1 * (tleaf - 298.15))
  temp_lim_1 = 1 + exp(0.3 * (tleaf - 318.15))
  temp_lim_2 = 1 + exp(0.2 * (286.15 - tleaf))
  
  ftemp = q10_func / (temp_lim_1 * temp_lim_2)
  
  ftemp
  
}

#### 2. collatz model function

collatz1992_func <- function(tleaf = 298.15,
                             ca = 400,
                             cica = 0.4,
                             par = 1000,
                             z = 0,
                             vcmax25 = 10,
                             model = 'clm',
                             kp = 200000,
                             use_ftemp = FALSE,
                             scale_kpvcmax = FALSE,
                             phi = 0.04){
  
  if(use_ftemp == FALSE){
    
    vcmax = vcmax25
    
  }else{
    
    if(model == 'clm'){
      
      vcmax = vcmax25 * ftemp_clm(tleaf)
      
    }else if(model == 'jules'){
      
      vcmax = vcmax25 * ftemp_jules(tleaf)
      
    }
    
  }
  
  if(scale_kpvcmax == FALSE){
    
    kp = kp
    
  }else{
    
    kp = 200000 * vcmax
    
  }
  
  if(model == 'pmodel'){
    
    phi = -0.008 + (0.00375 * (tleaf-273.15)) - (0.000058 * (tleaf-273.15) * (tleaf-273.15))
    
  }else{
    
    phi = phi
    
  }
  
  initial_slope = kp / calc_patm(z) # initial slope of co2 response curve
  
  ci = ca * cica * 1e-6 * calc_patm(z)
  
  a_p = initial_slope * ci
  a_c = vcmax
  a_j = phi * par
  
  assimilation = pmin(a_p, a_c, a_j)
  
  return(data.frame(assimilation, 
         a_p, 
         a_c, 
         a_j, 
         tleaf, 
         par, 
         ca, 
         cica, 
         ci, 
         phi,
         vcmax,
         vcmax25,
         initial_slope, 
         kp, 
         z, 
         model))
  
}


### 3. lsm-specific simulations in response to light, temperature, CO2
#### light
bethy_light <- collatz1992_func(par = seq(50, 2000, 50), phi = 0.04, model = 'bethy')
clm_light <- collatz1992_func(par = seq(50, 2000, 50), phi = 0.05, model = 'clm')
jules_light <- collatz1992_func(par = seq(50, 2000, 50), phi = 0.04, model = 'jules')
pmodel_light <- collatz1992_func(par = seq(50, 2000, 50), model = 'pmodel')

#### temperature
bethy_temperature <- collatz1992_func(tleaf = seq(278.15, 313.15, 5), 
                                phi = 0.04, model = 'bethy')
clm_temperature <- collatz1992_func(tleaf = seq(278.15, 313.15, 5), 
                              phi = 0.05, model = 'clm', use_ftemp = TRUE)
jules_temperature <- collatz1992_func(tleaf = seq(278.15, 313.15, 5), 
                                phi = 0.04, model = 'jules', use_ftemp = TRUE)
pmodel_temperature <- collatz1992_func(tleaf = seq(278.15, 313.15, 5), 
                                 model = 'pmodel')

### 4. make figures



