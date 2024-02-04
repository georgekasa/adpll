import numpy as np
import matplotlib.pyplot as plt
import random
import time as time_python



def fft(frequencies, till, sample_rate, duration):

    t = np.linspace(0, duration, int(sample_rate * duration), endpoint=False)  # Time array

    frequencies = frequencies[0:till]

    # Generate a signal by summing sine waves of the specified frequencies
    signal = np.sum(np.sin(2 * np.pi * frequencies[:, np.newaxis] * t), axis=0)
    fft_result = np.fft.fft(signal)
    fft_frequencies = np.fft.fftfreq(len(fft_result), 1/sample_rate)
    positive_frequencies_mask = fft_frequencies >= 0
    # Plot the magnitude spectrum
    plt.figure(figsize=(10, 4))
    plt.plot(fft_frequencies[positive_frequencies_mask], np.abs(fft_result[positive_frequencies_mask]))
    plt.title('FFT Result (Positive Frequencies)')
    plt.title('FFT Result')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude')
    plt.show()


def subplot(y, y_name):

    # Create a figure and a grid of subplots
    fig, axs = plt.subplots(3, 1, figsize=(10, 8))

    # Subplot 1
    axs[0].plot(y[0], color='blue')
    axs[0].set_title(y_name[0])
    axs[0].grid(True)

    # Subplot 2
    axs[1].plot(y[1], color='orange')
    axs[1].set_title(y_name[1])
    axs[1].grid(True)
    # Subplot 3
    axs[2].plot(y[2], color='green')
    axs[2].set_title(y_name[2])
    axs[2].grid(True)
    # Adjust layout for better spacing
    plt.tight_layout()

    # Show the plot
    plt.show()




start_time = time_python.time()
# Set a seed for reproducibility (you can use any integer value)
random_seed = 42
random.seed(random_seed)

# Set the range for the random value
lower_bound = 2.35e9
upper_bound = 2.45e9

# Generate a random value between the specified range

Fref = 48e6
Fdco = 2.420e9
Fdco_init = random.uniform(lower_bound, upper_bound)#random.uniform(lower_bound, upper_bound)#2.41e9
Tdco_init = 1 / Fdco_init

tdc_res = 10e-12  # TDC time resolution

Tref = 1 / Fref
Tdco = 1 / Fdco

run_time = 100e-6
run_steps = round(run_time / Tref)
time = np.arange(1, run_steps + 1) * Tref

Fstep = 40e3  # Is this realistic?
Kdco = Fstep

FCW_int_bits = np.ceil(np.log2(Fdco / Fref))
FCW_frac_bits = np.ceil(np.log2(Fref / Fstep))
FCW_nob = FCW_int_bits + FCW_frac_bits

FCW = round((Fdco / Fref) * 2**FCW_frac_bits)
# acc_phase = np.arange(1, run_steps + 1) * FCW
vpa_scaling = 2**FCW_frac_bits

lambda_c = [2**(-3), 2**(-3), 2**(-3), 2**(-4)]
alpha = 2**(-7)
rho = 2**(-15)



# Initialize internal variables
vpa_out = np.zeros(run_steps)
Fdco_out = np.zeros(run_steps)
Tdco_out = np.zeros(run_steps + 2)  # +1 is to account for the initial value
acc_freq_ratio = np.zeros(run_steps)  # +1 is to account for the initial value
tdc_time_out = np.zeros(run_steps)
tdc_in = np.zeros(run_steps)
tdc_in_time = np.zeros(run_steps + 1)
tdc_quant_out = np.zeros(run_steps)
tdc_scaling = np.zeros(run_steps)
tdc_phase_out = np.zeros(run_steps)
acc_freq_ratio_fb = np.zeros(run_steps + 1)
freq_err = np.zeros(run_steps + 1)  # +4 is to account for the initial 0 lpf inputs
phase_err = np.zeros(run_steps + 1)
phase_err_acc = np.zeros(run_steps + 1)
iir1_out = np.zeros(run_steps + 1)
iir2_out = np.zeros(run_steps + 1)
iir3_out = np.zeros(run_steps + 1)
iir4_out = np.zeros(run_steps + 1)
dco_ctrl = np.zeros(run_steps + 1)

# Initial Fdco_out and Tdco_out
Fdco_out[0] = Fdco_init
Tdco_out[0] = Tdco_init
acc_freq_ratio[0] = Tref/Tdco_out[0]

run_steps = run_steps-10

for i in range(1, run_steps + 1):
    # Accumulate frequency ratio between Fref and DCO output
    acc_freq_ratio[i] = acc_freq_ratio[i - 1] + Tref / Tdco_out[i - 1]
    
    # Calculate integer part of accumulated frequency ratio (sampled VPA output)
    vpa_out[i - 1] = np.floor(acc_freq_ratio[i])
    
    # TDC input is the fractional part of the accumulated frequency ratio
    tdc_in[i - 1] = acc_freq_ratio[i] - vpa_out[i - 1]
    
    # TDC input frequency ratio to time conversion - TODO: Add TDC noise
    tdc_in_time[i - 1] = tdc_in[i - 1] * Tdco_out[i - 1]
    
    # TDC quantized output (depends on TDC time resolution)
    tdc_quant_out[i - 1] = np.floor(tdc_in_time[i - 1] / tdc_res)
    
    # TDC output scaling calculation - TODO: In reality, this calculation
    # is not absolutely accurate and must be estimated
    tdc_scaling[i - 1] = (tdc_res / Tdco_out[i - 1]) * 2**FCW_frac_bits
    
    # Calculate accumulated feedback path frequency ratio: VPA (integer) +
    # TDC (fractional) frequency ratio values after scaling. Final value is
    # rounded to the closest integer value due to digital fixed-point 
    # representation 
    acc_freq_ratio_fb[i] = round(vpa_out[i - 1] * vpa_scaling + tdc_quant_out[i - 1] * tdc_scaling[i - 1])
    
    # Calculate frequency error: Difference between FCW (desired frequency
    # value) and the difference of the current and previous accumulated 
    # feedback path frequency ratios
    freq_err[i - 1] = FCW - (acc_freq_ratio_fb[i] - acc_freq_ratio_fb[i - 1])
    
    # Accumulate frequency error to calculate phase error
    phase_err[i] = phase_err[i - 1] + freq_err[i - 1]
    
    # Digital loop filter implementation: 4 cascaded single-stage IIR
    # filters + PI control. The integrated term is used to suppress DCO 
    # flicker noise
    iir1_out[i] = lambda_c[0] * phase_err[i] + (1 - lambda_c[0]) * iir1_out[i - 1]
    iir2_out[i] = lambda_c[1] * iir1_out[i] + (1 - lambda_c[1]) * iir2_out[i - 1]
    iir3_out[i] = lambda_c[2] * iir2_out[i] + (1 - lambda_c[2]) * iir3_out[i - 1]
    iir4_out[i] = lambda_c[3] * iir3_out[i] + (1 - lambda_c[3]) * iir4_out[i - 1]
    
    phase_err_acc[i] = phase_err_acc[i - 1] + phase_err[i]
    
    dco_ctrl[i] = round(alpha * iir4_out[i] + rho * phase_err_acc[i])
    
    # DCO output frequency is calculated as base DCO frequency + DCO
    # control word multiplied by DCO frequency gain in Hz/LSB - TODO: Add
    # DCO noise
    Fdco_out[i] = Fdco_out[0] + Kdco * dco_ctrl[i]
    
    # DCO output period calculation
    Tdco_out[i] = 1 / Fdco_out[i]



print("FDCO start (GHz): " + str(np.round(Fdco_init/1e9, 5)))
end_time = time_python.time()

# Calculate the elapsed time
elapsed_time = end_time - start_time
print(f"Elapsed Time: {elapsed_time:.4f} seconds")




# Plot DCO output frequency - TODO: Add ADPLL phase noise plot once noise
# sources are added

# plt.plot(time[1:run_steps], Fdco_out[1:run_steps], '.-')
# plt.xlabel('Time (s)')
# plt.ylabel('DCO Output Frequency (Hz)')
# plt.title('DCO Output Frequency')
# plt.grid(True)
# plt.show()






subplot([Fdco_out[1:run_steps], dco_ctrl[1:run_steps], phase_err[1:run_steps]], ["FDCO", "DCO control", "phase error"])
#fft(Fdco_out, 500, 20e9, run_time*0.5)