To run the program:
- python3 main.py
- Enter the desired velocity change in km/s (e.g.: 3.2)
- Enter the desired Moon theta (The Moon start position angle (degrees) in relation of the earth e.g.: 36)

* To run a simulation with a fixed moon, comment the '* time' in the get_moon_position_at_time function, as below:
```
def get_moon_position_at_time(time, initial_Θ =initial_Θ_moon):
    Θ_moon = initial_Θ + ω_moon #* time # Θ_moon
    x_moon = distance_earth_moon * np.cos(Θ_moon)
    y_moon = distance_earth_moon * np.sin(Θ_moon)
    return np.array([x_moon, y_moon])
```

* To create a ship trajectory in the same direction of the Moon orbit:
```
initial_position_ship = np.array([0.0, -initial_ship_radius])
```

* To create a ship trajectory in the opposite direction of the Moon orbit:
```
initial_position_ship = np.array([0.0, initial_ship_radius])
```

* To increase or decrease the simulation time, update t_span and second parameter of t_eval, setting the desired total simulation time in seconds:
```
t_span = (0, 850000)
t_eval = np.linspace(0, 850000, 50000)
```

* To increase or decrease the simulation quality, update t_eval third parameter, setting the desired number of calculation points:
```
t_eval = np.linspace(0, 850000, 50000)
```

* Results graph may be showed in fixed graphics or animations. Just comment and uncomment PLOT and ANIMATION sections in the end of the code
