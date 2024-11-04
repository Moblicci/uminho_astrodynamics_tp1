import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

###########################
##### EARTH CONSTANTS #####
###########################
μ_earth = 398600 #km^3/s^2
radius_earth = 6378 #km

###########################
##### EARTH CONSTANTS #####
###########################
μ_moon = 4903 #km^3/s^2
ω_moon = 2.6617e-6 #rad/s -> angular velocity => ω = (2 * π)/T where T = 27.32 * 24 * 3600 = 2360582 seconds

###########################
##### OTHER CONSTANTS #####
###########################
distance_earth_moon = 3.844e5 #Km

initial_ship_altitude = 300 #km
initial_ship_radius = initial_ship_altitude + radius_earth
delta_v = float(input("Digite o valor de delta V (km/s): "))
initial_moon_angle = float(input("Digite o ângulo inicial da lua em relação a terra (em graus): "))

initial_ship_velocity = np.sqrt(μ_earth/initial_ship_radius) + delta_v
initial_Θ_moon = np.radians(initial_moon_angle)

initial_position_ship = np.array([0.0, initial_ship_radius])
initial_velocity_ship = np.array([initial_ship_velocity, 0.0])
initial_conditions = [initial_position_ship[0], initial_position_ship[1], initial_velocity_ship[0], initial_velocity_ship[1]]

current_ship_velocity = initial_ship_velocity

#######################################################
################### MOON FUNCTIONS ####################
#######################################################
def get_moon_position_at_time(time, initial_Θ =initial_Θ_moon):
    Θ_moon = initial_Θ + ω_moon * time # Θ_moon
    x_moon = distance_earth_moon * np.cos(Θ_moon)
    y_moon = distance_earth_moon * np.sin(Θ_moon)
    return np.array([x_moon, y_moon])

#######################################################
####### FUNCTIONS TAKING THE SHIP AS REFERENCE #######
#######################################################
def get_scalar_distance(position_vector): #from ship to the body
    position_x, position_y = position_vector
    return np.sqrt(position_x**2 + position_y**2)

def get_gravitational_acceleration_vector(μ, body_position_vector):
    scalar_distance = np.linalg.norm(body_position_vector)
    if scalar_distance == 0:
        return np.array([0.0, 0.0])
    return -μ * body_position_vector / scalar_distance**3 # a = -μ * r / |r|^3

#######################################################
# PARAMETERS PREPARATION FOR EDO INTEGRATION FUNCTION #
#######################################################
def orbital_state_derivates(time, ship_orbital_space_vector):
    ship_position_x, ship_position_y, ship_velocity_x, ship_velocity_y = ship_orbital_space_vector
    earth_fixed_position = np.array([0.0, 0.0])

    r_earth_vector = np.array([ship_position_x, ship_position_y]) - earth_fixed_position
    gravitational_acceleration_vector_earth = get_gravitational_acceleration_vector(μ_earth, r_earth_vector)

    current_moon_position = get_moon_position_at_time(time)
    r_moon_vector = np.array([ship_position_x, ship_position_y]) - current_moon_position
    distance_to_moon = get_scalar_distance(r_moon_vector)

    if distance_to_moon < 200000: 
        gravitational_acceleration_vector_moon = get_gravitational_acceleration_vector(μ_moon, r_moon_vector)
    else:
        gravitational_acceleration_vector_moon = np.array([0.0, 0.0])

    gravitational_acceleration_x, gravitational_acceleration_y = (gravitational_acceleration_vector_earth + gravitational_acceleration_vector_moon)

    #returns position derivates (velocity) and velocity derivates (acceleration) to be used by the EDO solver
    return [ship_velocity_x, ship_velocity_y, gravitational_acceleration_x, gravitational_acceleration_y] # F = m * a

t_span = (0, 850000)
t_eval = np.linspace(0, 850000, 50000)

solution = solve_ivp(orbital_state_derivates, t_span, initial_conditions, method='DOP853', t_eval=t_eval, rtol=1e-12, atol=1e-14)

#######################
##### PLOT ############
#######################

###################################
### Show Trajectory apogee info ###
###################################
distances_from_earth = np.sqrt(solution.y[0]**2 + solution.y[1]**2)
apogee_index = np.argmax(distances_from_earth)
apogee_position = np.array([solution.y[0][apogee_index], solution.y[1][apogee_index]])
apogee_velocity = [solution.y[2][apogee_index], solution.y[3][apogee_index]]
apogee_distance = distances_from_earth[apogee_index]
apogee_gravitational_force = get_gravitational_acceleration_vector(μ_earth, apogee_position)
kinetic_energy = 0.5 * (apogee_velocity[0]**2 + apogee_velocity[1]**2)
potential_energy = -μ_earth / apogee_distance
specific_energy = kinetic_energy + potential_energy
print("### Variáveis no Apogeu ###")
print(f"Posição (x, y): {apogee_position}")
print(f"Distância da Terra: {apogee_distance:.2f} km")
print(f"Velocidade (vx, vy): {apogee_velocity}")
print(f"Força Gravitacional: {apogee_gravitational_force}")
print(f"Velocidade Total: {np.sqrt(apogee_velocity[0]**2 + apogee_velocity[1]**2):.3f} km/s")
print(f"Energia Cinética Específica: {kinetic_energy:.3f} km^2/s^2")
print(f"Energia Potencial Específica: {potential_energy:.3f} km^2/s^2")
print(f"Energia Específica Total: {specific_energy:.3f} km^2/s^2")

############
### PLOT ###
############
# ship_velocity_x = solution.y[2]
# ship_velocity_y = solution.y[3]
# current_ship_velocity = np.round(np.sqrt(ship_velocity_x**2 + ship_velocity_y**2), 3)
# print(current_ship_velocity)
# plt.figure(figsize=(12, 6))

# plt.subplot(1, 2, 1)
# plt.plot(solution.y[0], solution.y[1], label="Ship trajectory")
# plt.plot(0, 0, 'o', label="Earth", markersize=10)
# plt.plot(get_moon_position_at_time(t_eval)[0], get_moon_position_at_time(t_eval)[1], 'o', label="Moon (trajectory)", markersize=5, alpha=0.3)
# plt.xlabel("Distance in x (km)")
# plt.ylabel("Distance in y (km)")
# plt.legend()
# plt.axis("equal")

# print("Current Ship Velocity:", current_ship_velocity)

# plt.subplot(1, 2, 2)
# plt.plot(solution.t, current_ship_velocity, label="Ship velocity", color="orange")
# plt.xlabel("Time (s)")
# plt.ylabel("Velocity (km/s)")
# plt.legend()

# plt.tight_layout()
# plt.show()

#################
### Animation ###
#################
ship_position_x = solution.y[0]
ship_position_y = solution.y[1]
ship_velocity_x = solution.y[2]
ship_velocity_y = solution.y[3]
current_ship_velocity = np.round(np.sqrt(ship_velocity_x**2 + ship_velocity_y**2), 3)
print(current_ship_velocity)

moon_positions = np.array([get_moon_position_at_time(t) for t in solution.t])
moon_position_x = moon_positions[:, 0]
moon_position_y = moon_positions[:, 1]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.set_xlim(min(ship_position_x.min(), moon_position_x.min()) - 10000, 
             max(ship_position_x.max(), moon_position_x.max()) + 10000)
ax1.set_ylim(min(ship_position_y.min(), moon_position_y.min()) - 10000, 
             max(ship_position_y.max(), moon_position_y.max()) + 10000)
ax1.set_xlabel("Distância X (km)")
ax1.set_ylabel("Distância Y (km)")
ax1.set_title("Trajetória da Nave e da Lua")
nave_traj, = ax1.plot([], [], 'b-', label="Trajetória da Nave")
lua_traj, = ax1.plot([], [], 'gray', label="Trajetória da Lua")
terra = ax1.plot(0, 0, 'go', label="Terra", markersize=10)[0]
nave_ponto, = ax1.plot([], [], 'ro', markersize=8, label="Posição da Nave")
lua_ponto, = ax1.plot([], [], 'yo', markersize=8, label="Posição da Lua")

ax2.set_xlim(t_span[0], t_span[1])
ax2.set_ylim(min(current_ship_velocity), max(current_ship_velocity))
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Velocity (km/s)")
ax2.set_title("Velocidade da Nave")
velocidade_line, = ax2.plot([], [], color="orange", label="Velocidade da Nave")

def init():
    nave_traj.set_data([], [])
    lua_traj.set_data([], [])
    nave_ponto.set_data([], [])
    lua_ponto.set_data([], [])
    velocidade_line.set_data([], [])
    return nave_traj, lua_traj, nave_ponto, lua_ponto, velocidade_line

def atualizar(frame):
    nave_traj.set_data(ship_position_x[:frame], ship_position_y[:frame])
    lua_traj.set_data(moon_position_x[:frame], moon_position_y[:frame])

    nave_ponto.set_data([ship_position_x[frame]], [ship_position_y[frame]])
    lua_ponto.set_data([moon_position_x[frame]], [moon_position_y[frame]])

    velocidade_line.set_data(solution.t[:frame], current_ship_velocity[:frame])

    return nave_traj, lua_traj, nave_ponto, lua_ponto, velocidade_line

frame_step = 50
frames = range(0, len(solution.t), frame_step)

ani = FuncAnimation(fig, atualizar, frames=frames, init_func=init, interval=30, blit=True)

plt.tight_layout()
plt.legend()
plt.grid(True)
plt.show()