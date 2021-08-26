# This source code is related to the Handbook of Quenchants and Quenching Technology (ASM International)
# The mathematical model of Heat Transfer is represented in the chapter
# Heat Transfer during Quenching Processes (Solving method for heat conduction equation)
#
# The Heat Transfer model of a cylindrical solid object is considered
# Due to the symmetry of the object , the temopro-spatial distribution of temperature
# at the half of the cross section of the cylinder is calculated. Therefore, we are using a 1D axis-symmetrical model
#
# The transient heat transfer during quenching can mathematically be described by an appropriate
# form of a Fourier’s heat conduction equation (Eq.2). For the 1D axis-symmetrical model, the (Eq.8)
# was taken into account. The solution of (Eq.8) is obtained by the Finite Difference Method (FDM).
#
# The variables and constants ar the followings:
#  nn               The number of the local coordinate nodes
#  radius           The radius of the cylindrical object (s)
#  htc              The value of Heat Transfer Coefficient (Wm-2k-1)
#  t_start          The starting temperature (C))
#  t_quenchant      The temperature of the quenchant (C)
#  time_            The time elapsed (s)
#  dt               The time step (s)
#  dx               The space interval between the local coordinates (m)
#  k                The heat conductivity  (Wm-1K-1)
#  Cp               The specific heat (JKg-1K-1)
#  rho              The density in (Kgm-3)
#  alpha            The thermal diffusivity (Jm-3K-1)
#
#  temperature_vec  The vector of recent temperature distribution at the half cross section of the cylinder
#  exchange_vec     The vector for former temperature distribution at the half cross section of the cylinder
#
#
# Notes:
# This Python code require the module "Numpy".
# Please install this module
#
#

import numpy as np

nn: int = 10
radius: float = 0.01

htc: float = 100
t_start: float = 850
t_quenchant: float = 20

dx: float = radius / nn
dt: float = 0.01
time_: float = 0



# Materials Data --------
k: float = 27
cp: float = 485
rho: float = 7850
alpha = k / (cp * rho)
# Materials Data --------


temperature_vec = np.arange(nn + 1)
temperature_vec = temperature_vec.astype("float")

exchange_vec = np.arange(nn + 1)
exchange_vec = exchange_vec.astype("float")


# Initialization of start temperature---------
for i in range(0, nn+1):
    temperature_vec[i] = t_start
    exchange_vec[i] = t_start
# Initialization of start temperature---------

# Initialization of data file---------
result_file = open('1D_HT_results.txt','w')
str_coordinates = ' '
for ii in range(nn+1):
    str_coordinates = str_coordinates + ' ' + str(round(ii*dx,5))
str_coordinates = str_coordinates + '\n'
result_file.write(str_coordinates)
# Initialization of data file---------



for j in range(0, 6000):
    for i in range(1, nn):
        temperature_vec[i] = exchange_vec[i] + dt * alpha * (1 / (dx*dx) * (exchange_vec[i-1] + exchange_vec[i+1] - 2 * exchange_vec[i]) + 1 / (i*dx) * 1 / (2*dx) * (exchange_vec[i+1] - exchange_vec[i-1]))
    temperature_vec[0] = exchange_vec[0] + dt * alpha * (1 / (dx * dx) * 2 * (exchange_vec[1] - exchange_vec[0]))
    temperature_vec[nn] = exchange_vec[nn] + dt * alpha * (1 / (dx * dx) * 2 * (exchange_vec[nn - 1] - exchange_vec[nn] - dx / k* (htc * (exchange_vec[nn] - t_quenchant))) + 1 / (nn * dx) * (-1 / k) * (htc * (exchange_vec[nn] - t_quenchant)))

    time_ = time_ + dt
    exchange_vec = temperature_vec
    # print(time_, ' ', temperature_vec

    str_temperature = ''
    for ii in range(nn):
        str_temp_round = round(temperature_vec[ii],2)
        str_temperature = str_temperature + ' ' + str(str_temp_round)
    str_temperature = str_temperature + '\n'
    result_file.write(str_temperature)




# Closing data file---------
result_file.close()
# Closing data file---------
