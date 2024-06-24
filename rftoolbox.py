import numpy as np
import regex as re
import matplotlib.pyplot as plt

class Signal:
    def __init__(self, time, amplitude):
        if len(time) < 2 and len(time) != len(amplitude):
            raise(ValueError)

        self.time = time
        self.amplitude = amplitude
        self.dt = time[1] - time[0]
        self.fs = 1 / self.dt
        self.length = len(time)

class Network:
    def __init__(self, filename):
        m = re.match(r".*s([0-9]+)p$", filename)
        if not m:
            raise Exception("invalid file type")

        self.ports = int(m.group(1))

        # Import data
        data = []
        nb_columns = None

        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()

                if not line or line.startswith("!"):
                    continue

                # Parse options
                if line.startswith("#"):
                    m = re.match(r"#\s+([A-Za-z]+)\s+S\s+([A-Z]+)\s+([A-Z]+)\s+([0-9.]+)$", line)
                    if not m:
                        raise(ValueError)

                    (unit, format, _, z0) = m.groups()
                    self.z0 = float(z0)

                    continue

                # Convert to float
                tmp = [float(i) for i in line.split()]

                # Use number of columns of first data line
                if nb_columns is None:
                    nb_columns = len(tmp)

                if len(tmp) < nb_columns:
                    data[-1].extend(tmp)
                else:
                    data.append(tmp)

        # No options found
        if not 'unit' in locals():
            raise(ValueError)

        # Scale frequency vector
        match unit:
            case 'Hz':
                scale = 1

            case 'kHz':
                scale = 1000

            case 'MHz':
                scale = 1000000

            case 'GHz':
                scale = 1000000000

            case _:
                raise(ValueError)

        freqs = np.zeros(len(data[:]))
        s = np.zeros((len(data[:]), self.ports, self.ports), dtype=complex)

        # Transform data
        for f in range(len(data[:])):
            freqs[f] = data[f][0] * scale
            for i in range(self.ports):
                for j in range(self.ports):
                    tmp = "s" + str(i+1) + str(j+1)

                    # Set column names and convert data
                    match format.upper():
                        case 'DB':
                            magnitude = data[f][2 * (i * self.ports + j) + 1]
                            magnitude = np.power(10, magnitude / 20)

                            angle = data[f][2 * (i * self.ports + j) + 2]

                            s[f][i][j] = magnitude * np.cos(2 * np.pi * angle / 360) + magnitude * np.sin(2 * np.pi * angle / 360) * 1j

                        case 'MA':
                            magnitude = data[f][2 * (i * self.ports + j) + 1]
                            angle = data[f][2 * (i * self.ports + j) + 2]

                            s[f][i][j] = magnitude * np.cos(2 * np.pi * angle / 360) + magnitude * np.sin(2 * np.pi * angle / 360) * 1j

                        case 'RI':
                            real = data[f][2 * (i * self.ports + j) + 1]
                            imag = data[f][2 * (i * self.ports + j) + 2]

                            s[f][i][j] = real + imag * 1j

                        case _:
                            raise(ValueError)

        # Set class values
        self.__s = s
        self.f = freqs
        self.df = freqs[1] - freqs[0]

        # Create abcd matrix, if port number is 2
        if self.ports == 2:
            self.__s_to_abcd()

    def s(self, port0, port1, frequencies=slice(None)):
        if port0 > self.ports or port0 < 1:
            raise(ValueError)

        if port1 > self.ports or port1 < 1:
            raise(ValueError)

        return self.__s[frequencies, port0 - 1, port1 - 1]

    def impulse_response(self, port0, port1):
        if port0 > self.ports or port0 < 1:
            raise(ValueError)

        if port1 > self.ports or port1 < 1:
            raise(ValueError)

        # Transform back to time domain
        h = np.fft.irfft(self.__s[:,port0 - 1, port1 - 1])

        # Get sampling frequency
        fs = 2 *self.f[-1]

        # Create time vector
        t = np.linspace(0, (len(h) - 1) / fs, len(h))

        return Signal(t, h)

    def __resample_impulse_response(self, port0, port1, fs):
        if port0 > self.ports or port0 < 1:
            raise(ValueError)

        if port1 > self.ports or port1 < 1:
            raise(ValueError)

        # Calculate maximum bandwidth
        fmax = fs / 2

        if fs % self.df:
            raise(ValueError)

        # Copy s-parameters
        s = self.__s[:,port0 - 1, port1 - 1]

        # Append zeros
        nb_zeros = int( (fmax - self.f[-1]) / self.df )
        s = np.append(s, np.zeros(nb_zeros) )

        # Calculate impulse response
        h = np.fft.irfft(s)

        # Create time vector
        t = np.linspace(0, (len(h) - 1) / fs, len(h))

        return Signal(t, h)

    def step_response(self, port0, port1):
        if port0 > self.ports or port0 < 1:
            raise(ValueError)

        if port1 > self.ports or port1 < 1:
            raise(ValueError)

        # Calculate impulse response
        impulse_response = self.impulse_response(port0, port1)

        # Calculate tdr
        # Convolve with step function, but zeros can be ignored since the result is zero, so the cumsum can be used
        a = np.cumsum(impulse_response.amplitude)

        # Create time vector, divide by two because of reflection (double length)
        t = np.linspace(0, impulse_response.dt / 2 * ( len(a) - 1 ), len(a))

        return Signal(t, a)

    def calculate_port_impedance(self, port0, port1, Z0):
        a = self.step_response(port0, port1)
        Z = Z0 * (1 + a.amplitude) / ( 1 - a.amplitude)

        return Signal(a.time, Z)

    def stimulate_port(self, port0, port1, vector : Signal, min_oversampling = 16, max_oversampling = 128):
        # Calculate target sampling frequency
        # oversampling = int( np.ceil(2 * self.f[-1] / vector.fs) )
        oversampling = min_oversampling

        smallest_dev = np.inf
        best_match = 0

        while True:
            remainder = oversampling * vector.fs % self.df
            if remainder:
                if smallest_dev > remainder:
                    smallest_dev = remainder
                    best_match = oversampling
                oversampling += 1
            else:
                break

            # Limit loops
            if oversampling == max_oversampling:
                oversampling = best_match
                break

        # Calculate (re-) sampling frequency
        fs = np.ceil( (oversampling * vector.fs) / self.df ) * self.df

        # Calculate impulse response
        h = self.__resample_impulse_response(port0, port1, fs)

        # Calculate simulated data_rate
        simulated_data_rate = fs / oversampling

        # Create time vector
        t = np.linspace(0, vector.length / simulated_data_rate, oversampling * vector.length)

        # Upsamling
        input = Signal(t, np.repeat(vector.amplitude, oversampling))

        # Convolution
        output = Signal(t, np.convolve(input.amplitude, h.amplitude))

        return simulated_data_rate, h, input, output

    # S to ABCD matrix transformation
    def __s_to_abcd(self):
        self.a = ( (1 + self.__s[:,0,0]) * (1 - self.__s[:,1,1]) + self.__s[:,0,1] * self.__s[:,1,0] ) / (2*self.__s[:,1,0])
        self.b = self.z0 * ( (1 + self.__s[:,0,0]) * (1 + self.__s[:,1,1]) - self.__s[:,0,1] * self.__s[:,1,0] ) / (2*self.__s[:,1,0])
        self.c = 1 / self.z0 * ( (1 - self.__s[:,0,0]) * (1 - self.__s[:,1,1]) - self.__s[:,0,1]*self.__s[:,1,0] ) / (2.*self.__s[:,1,0])
        self.d = ( (1 - self.__s[:,0,0]) * (1 + self.__s[:,1,1]) + self.__s[:,0,1] * self.__s[:,1,0] ) / (2*self.__s[:,1,0])


    # ABCD to S matrix transformation
    def __abcd_to_s(self, z01, z02):
        s = self.__s.copy()

        dAz = self.a*z02 + self.b + self.c*z01*z02 + self.d*z01

        s[:,0,0] = 1/dAz * (self.a *z02 + self.b - self.c *z01*z02 - self.d*z01)
        s[:,0,1] = 1/dAz * (2 * np.sqrt(z01) * np.sqrt(z02) * (self.a*self.d - self.b*self.c))
        s[:,1,0] = 1/dAz * (2 * np.sqrt(z01) * np.sqrt(z02))
        s[:,1,1] = 1/dAz * (-self.a*z02 + self.b - self.c*z01*z02 + self.d*z01)

        return s

    # Z01: output port impedance, Z02: input port impedance
    def impedance_transformaton(self, z01, z02):
        if self.ports == 2:
            s = self.__abcd_to_s(z01, z02)
            return s
        else:
            raise Exception("number of ports must be 2")

# End of class 'Network'

def create_digital_signal_vector(data_rate, nb_symbols, output_swing = 1.0, pulse_amplitudes = 2):
    random_integers = np.random.random_integers(0, pulse_amplitudes - 1, nb_symbols)
    bins = output_swing * (np.linspace(0, 1.0, pulse_amplitudes) - 0.5)
    vector = bins[random_integers]

    t = np.linspace(0, (nb_symbols - 1) / data_rate, nb_symbols)

    return Signal(t, vector)

def create_eye_diagram(h : Signal, input : Signal, output : Signal, resolution = 64, xticks = 4, yticks = 4):
    # Shift the input vector
    time_shift = np.argmax(h.amplitude)

    # input.time = input.time + input.dt * time_shift

    diff = np.diff(input.amplitude)
    peaks = np.where(np.abs(diff) > 0.0)

    peaks = np.reshape(peaks, np.size(peaks))
    peaks = peaks + time_shift + 1 # because the peak occurs one sample before the peak (diff operation)

    minimal_symbol_length = np.amin(np.diff(peaks))

    # Split vector and drop first symbol (settle time)
    symbols = np.split(output.amplitude, peaks)
    symbols = symbols[1:-1]

    # Create eye-diagram
    max_amplitude = np.around( np.amax( np.abs(output.amplitude) ), decimals = 1)
    # min_amplitude = np.amin(output.amplitude)

    eye_diagram = np.zeros( (minimal_symbol_length, resolution) )

    quantization_steps = np.flip(np.linspace(-max_amplitude, max_amplitude, resolution))

    for symbol in symbols:
        flipped_symbol = np.flip(symbol)
        for x in range(minimal_symbol_length):
            # Left aligned symbol
            y = np.digitize(symbol[x], quantization_steps) # index is ones-based
            eye_diagram[x, y - 1] += 1

            # Centered symbol
            shift = int(np.floor( (len(symbol) - minimal_symbol_length) / 2 ))
            y = np.digitize(symbol[x+shift], quantization_steps) # index is ones-based
            eye_diagram[x, y - 1] += 1

            # Right aligned symbol
            y = np.digitize(flipped_symbol[x], quantization_steps) # index is ones-based
            eye_diagram[minimal_symbol_length - x - 1, y - 1] += 1

    # Get eye-parameters (height, width)
    inv_hit_mask = np.zeros( (minimal_symbol_length, resolution) )
    index = np.where(eye_diagram == 0)
    inv_hit_mask[index] = 1

    # eye-height
    max_height = 0
    height = 0
    for i in range(resolution):
        if inv_hit_mask[int(minimal_symbol_length/2), i] == 1:
            height += 1
        else:
            if max_height < height:
                max_height = height
                height = 0

    if max_height < height:
        max_height = height

    # print('eye height: %1.4fV' % (max_height*2*max_amplitude/resolution))

    #eye-width
    max_width = 0
    width = 0
    for i in range(minimal_symbol_length):
        if inv_hit_mask[i, int(resolution/2)] == 1:
            width += 1
        else:
            if max_width < width:
                max_width = width
                width = 0

    if max_width < width:
        max_width = width

    # print('eye width: %1.4f%%' % (max_width/minimal_symbol_length*100))

    # Create scales
    x_scale = np.zeros( (xticks, 2) )
    y_scale = np.zeros( (yticks, 2) )

    x_scale[:, 0] = np.linspace(0, minimal_symbol_length - 1, xticks)
    x_scale[:, 1] = np.around( np.linspace(0, 1, xticks), decimals=2)

    y_scale[:, 0] = np.linspace(0, resolution - 1, yticks)
    y_scale[:, 1] = np.flip( np.around( np.linspace(-max_amplitude, max_amplitude, yticks) , decimals=2) )

    plt.imshow(eye_diagram.T, cmap='jet', interpolation='nearest')
    plt.xticks(x_scale[:,0],x_scale[:,1])
    plt.yticks(y_scale[:,0],y_scale[:,1])
    plt.xlabel('t / symbol')
    plt.ylabel('voltage')
    plt.title('eye diagram')
    plt.show()