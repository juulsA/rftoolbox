# rftoolbox

The *rftoobox* helps you to read *Touchstone* files (*.sXp, X = number of ports) and to perform the necessary calculations to evaluate the measurement results.

## The *Signal* class

When a time-dependent vector is returned by a function, the *rftoolbox* uses a *Signal* consisting of:
- a time vector (`time`)
- a signal or amplitude vector (`amplitude`)
- the time between the samples (`dt`) and the sampling frequency (`fs`)
- and the signal length (`length`)

## The *Network* class

### Reading a touchstone file

The *Network* class takes a filename as input argument, reads the data and - dependent on the format - converts it into complex numbers for better processing.

```
import rftoolbox as rftbx

network = rftbx.Network("example.s2p")
```

### *Example:* Plot $S_{21}$ - Parameter

The *Network* class allows direct access to the s-parameters through `S(port0, port1)`. The frequency or the frequencies at which the s-parameters should be returned can be passed as an optional argument: `S(port0,port1, (1,3,100)`). Otherwise all frequency points will be returned.

```
import numpy as np
import matplotlib.pyplot as plt
import rftoolbox as rftbx

network = rftbx.Network("example.s2p")

#  Plot S21
fig, ax = plt.subplots(figsize=(15, 10))
ax.plot(network.f/1e9, 20*np.log10(np.abs(network.S(2,1))))

# add labels
ax.set_title('transfer function')
ax.set_xlabel('f / GHz')
ax.set_ylabel('H / dB')
ax.set_xlim(0,20)
ax.set_ylim(-80,0)

ax.grid(True, which='major', color='#C0C0C0', linestyle='-')
ax.grid(True, which='minor', color='#C5C9C7', linestyle='--')

ax.minorticks_on()
```

![transfer function](/img/transfer_function.png)

### *Example:* Calculate impulse response

The impulse response of the s-parameter can be calculated in the same manner the s-parameters are accessible.

```
import matplotlib.pyplot as plt
import rftoolbox as rftbx

network = rftbx.Network("example.s2p")

# Impulse response of S21
impulse_response = network.impulse_response(2,1)

fig, ax = plt.subplots(figsize=(15, 10))

samples = int( 10e-9 / impulse_response.dt )
ax.plot(impulse_response.time[0:samples]/1e-9, impulse_response.amplitude[0:samples])

# add labels
ax.set_title('impulse response')
ax.set_xlabel('t / ns')
ax.set_ylabel('1 / s')
ax.set_xlim(0,10)
ax.set_ylim(-0.01,0.12)

ax.grid(True, which='major', color='#C0C0C0', linestyle='-')
ax.grid(True, which='minor', color='#C5C9C7', linestyle='--')

ax.minorticks_on()
```

![impulse response](/img/impulse_response.png)

### Example: Simulate data transmission and eye diagram

Although the s-parameters contain all the information needed to determine the quality a network, it is easier to evaluate a lossy channel based on the shape of the received signal and the resulting eye diagram. Therefore, you can stimulate the port of interest by any signal or use the built-in function to generate a random digital signal. The signal must be of type `Signal`. 

--- 

The `create_digital_signal_vector` takes the baud rate as the first input argument and the number of symbols as the second. The third argument is the optional peak-to-peak output voltage and 1.0 by default. The fourth argument describes the number of amplitude steps, if the signal is an *ASK* modulated signal. This argument is also optional and 2 by default. 

--- 

The `stimulate_port` function takes the two ports as the first two inputs. The third input is the signal used to stimulate the port. Since the sampling rate of the signal does not always match to the sampling rate of the impulse response, the impulse response is resampled based on a third optional argument - the minimum oversampling. The value is 16 by default. The minimum oversampling is increased by one until the remainder of the oversampling times the sampling rate of the signal divided by the frequency step of the impulse response is zero (`oversampling * vector.fs % self.df`) or the oversampling reaches the limit of `max_oversampling` or 128 by default. Then the oversampling with the smallest remainder is used.

```
import numpy as np
import matplotlib.pyplot as plt
import rftoolbox as rftbx

# Set transmission parameters
data_rate = 10.3125e9
output_swing = 2.4
nb_symbols = 128

network = rftbx.Network("example.s2p")

# Simulate channel
vector = rftbx.create_digital_signal_vector(data_rate, nb_symbols, output_swing, 2)
simulated_data_rate, h, input, output = network.stimulate_port(2, 1, vector, 64)

# Ploting
dt = h.dt
samples = int( 10e-9 / dt )

fig1, ax = plt.subplots(2, figsize=(15,10))

ax[0].plot(input.time[0:samples] / 1e-9, input.amplitude[0:samples])
ax[0].set_title('input vector @ %1.4f Gbit' % (simulated_data_rate / 1e9))
ax[0].set_xlabel('time / ns')
ax[0].set_ylabel('voltage')
ax[0].set_xlim(0,10)

ax[1].plot(output.time[0:samples] / 1e-9, output.amplitude[0:samples])
ax[1].set_title('output vector @ %1.4f Gbit' % (simulated_data_rate / 1e9))
ax[1].set_xlabel('time / ns')
ax[1].set_ylabel('voltage')
ax[1].set_xlim(0,10)

# Set grid
ax[0].grid(True, which='major', color='#C0C0C0', linestyle='-')
ax[0].grid(True, which='minor', color='#C5C9C7', linestyle='--')
ax[1].grid(True, which='major', color='#C0C0C0', linestyle='-')
ax[1].grid(True, which='minor', color='#C5C9C7', linestyle='--')

ax[0].minorticks_on()
ax[1].minorticks_on()

fig1.tight_layout(pad=2.0)
plt.show()
```

![transmission](/img/transmission.png)

--- 

To create an eye diagram, the following line needs to be added to the previous code example:

```
rftbx.create_eye_diagram(h, input, output, 128, 5, 13)
```

The `create_eye_diagram` function takes the impulse response *h*, the input vector and output vector as arguments. The fourth argument is the resolution of the y-axis and by default 64. The last two arguments specify the number of points on the x- and y-scale.

![eye diagram](/img/eye_diagram.png)

#### Example of PAM4 signal

![eye diagram](/img/eye_diagram_pam4.png)

### Example: Plot port impedance

To locate discontinuities in the impedance, it is useful to plot the impedance over time.
This can be done with the function `calculate_port_impedance`. The first two arguments are the ports of interest and the third argument is the reference impedance of your measurement equipment. Although the reference impedance is included in the option line of the *Touchstone* file, the impedance must be passed manually because in case of a differential measurement, the reference impedance is twice times the reference impedance and does not match the reference impedance of the equipment. 

```
import matplotlib.pyplot as plt
import rftoolbox as rftbx

network = rftbx.Network("example.s2p")

impedance = network.calculate_port_impedance(1, 1, 50)

fig, ax = plt.subplots(figsize=(15,10))

# Plot frequency
samples = int(2e-9 / impedance.dt)

ax.plot(impedance.time[0:samples], impedance.amplitude[0:samples])
ax.set_title('Impedance')
ax.set_xlabel('time / ns')
ax.set_ylabel('Z / Ohms')
ax.set_xlim(0, 2e-9)

ax.grid(True, which='major', color='#C0C0C0', linestyle='-')
ax.grid(True, which='minor', color='#C5C9C7', linestyle='--')

ax.minorticks_on()
```

![port impedance](/img/impedance.png)

### Example: Impedance transformation

Sometimes you need to measure a 2-port-network with a specific input and output impedance, but in most cases the measurement equipment only supports one reference impedance, typically 50 ohms. Therefore, you can perform your measurement as usual and then use the `impedance_transformation` function.

The function takes two input arguments: the source impedance and the load impedance as the second.

```
import numpy as np
import matplotlib.pyplot as plt
import rftoolbox as rftbx

network = rftbx.Network("example_2.s2p")

transformed = network.impedance_transformaton(75,25)

fig, ax = plt.subplots(2, 2, figsize=(15,10))

ax[0,0].plot(network.f/1e9, 20*np.log10(np.abs(network.s(1,1))))
ax[0,0].plot(network.f/1e9, 20*np.log10(np.abs(transformed[:,0,0])))

ax[0,1].plot(network.f/1e9, 20*np.log10(np.abs(network.s(2,1))))
ax[0,1].plot(network.f/1e9, 20*np.log10(np.abs(transformed[:,1,0])))

ax[1,0].plot(network.f/1e9, 20*np.log10(np.abs(network.s(1,2))))
ax[1,0].plot(network.f/1e9, 20*np.log10(np.abs(transformed[:,0,1])))

ax[1,1].plot(network.f/1e9, 20*np.log10(np.abs(network.s(2,2))))
ax[1,1].plot(network.f/1e9, 20*np.log10(np.abs(transformed[:,1,1])))

ax[0,0].set_xlim(2, 17)
ax[0,1].set_xlim(2, 17)
ax[1,0].set_xlim(2, 17)
ax[1,1].set_xlim(2, 17)


ax[0,0].set_title('$S_{11}$')
ax[0,1].set_title('$S_{12}$')
ax[1,0].set_title('$S_{21}$')
ax[1,1].set_title('$S_{22}$')

ax[0,0].set_xlabel('f / GHz')
ax[0,0].set_ylabel('S / dB')
ax[0,0].legend(['$Z_{out}$ = $Z_{in}$ = 50$\,$ohms', '$Z_{out}$ = 75$\,$ohms, $Z_{in}$ = 25$\,$ohms'])

ax[0,1].set_xlabel('f / GHz')
ax[0,1].set_ylabel('S / dB')
ax[0,1].legend(['$Z_{out}$ = $Z_{in}$ = 50$\,$ohms', '$Z_{out}$ = 75$\,$ohms, $Z_{in}$ = 25$\,$ohms'])

ax[1,0].set_xlabel('f / GHz')
ax[1,0].set_ylabel('S / dB')
ax[1,0].legend(['$Z_{out}$ = $Z_{in}$ = 50$\,$ohms', '$Z_{out}$ = 75$\,$ohms, $Z_{in}$ = 25$\,$ohms'])

ax[1,1].set_xlabel('f / GHz')
ax[1,1].set_ylabel('S / dB')
ax[1,1].legend(['$Z_{out}$ = $Z_{in}$ = 50$\,$ohms', '$Z_{out}$ = 75$\,$ohms, $Z_{in}$ = 25$\,$ohms'])

ax[0,0].grid(True, which='major', color='#C0C0C0', linestyle='-')
ax[0,0].grid(True, which='minor', color='#C5C9C7', linestyle='--')

ax[0,1].grid(True, which='major', color='#C0C0C0', linestyle='-')
ax[0,1].grid(True, which='minor', color='#C5C9C7', linestyle='--')

ax[1,0].grid(True, which='major', color='#C0C0C0', linestyle='-')
ax[1,0].grid(True, which='minor', color='#C5C9C7', linestyle='--')

ax[1,1].grid(True, which='major', color='#C0C0C0', linestyle='-')
ax[1,1].grid(True, which='minor', color='#C5C9C7', linestyle='--')

```

![impedance transformation](/img/impedance_transformation.png)