import matplotlib.pyplot as plt
import numpy as np

# Define physical error rates p
p_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 5e-05]

# Running times and errors
surface_time_stim = [0.033, 0.042, 0.013, 0.051, 1.19, 6.96, 503, 515]
surface_error_stim = [0.007, 0.005, 0.004, 0.035, 0.98, 0.25, 39, 148]
surface_time_our = [2.85, 2.37, 2.16, 2.13, 3.01, 2.89, 2.87, 3.15]
surface_error_our = [0.08, 0.01, 0.11, 0.07, 0.31, 0.37, 0.36, 0.34]

repetition_time_stim = [0.0045, 0.008, 0.099, 0.56, 32.9, 269, None, None]
repetition_error_stim = [0.0004, 0.003, 0.014, 0.13, 4.4, 40, None, None]
repetition_time_our = [0.735, 0.711, 0.797, 0.781, 0.786, 0.822, 0.789, 0.787]
repetition_error_our = [0.09, 0.16, 0.25, 0.24, 0.19, 0.89, 0.25, 0.28]

square_time_stim = [0.06, 0.041, 0.015, 0.031, 0.68, 5.26, 290, 200]
square_error_stim = [0.01, 0.005, 0.002, 0.007, 0.21, 0.61, 19, 68]
square_time_our = [4.34, 3.45, 2.73, 2.76, 3.31, 4.89, 4.82, 4.40]
square_error_our = [0.26, 0.07, 0.04, 0.06, 0.10, 0.47, 0.28, 0.58]

hexagon_time_stim = [0.11, 0.071, 0.024, 0.03, 0.65, 4.84, 421, 494]
hexagon_error_stim = [0.019, 0.012, 0.005, 0.004, 0.23, 1.12, 17, 106]
hexagon_time_our = [7.00, 5.77, 4.95, 4.23, 4.68, 5.92, 8.99, 9.29]
hexagon_error_our = [0.05, 0.06, 1.23, 0.14, 0.12, 0.53, 2.37, 1.86]

# Create a figure with four subplots
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
axs = axs.ravel()

# Bar width and positions
bar_width = 0.35
x = np.arange(len(p_values))
x_stim = x - bar_width / 2
x_our = x + bar_width / 2

# Plot for Surface code
axs[0].bar(x_stim, surface_time_stim, width=bar_width, label='Surface (STIM)', color='#aec7e8', hatch='//', edgecolor='black')
axs[0].bar(x_our, surface_time_our, width=bar_width, label='Surface (Our)', color='#1f77b4', edgecolor='black')
axs[0].errorbar(x_stim, surface_time_stim, yerr=surface_error_stim, fmt='none', ecolor='black', capsize=3)
axs[0].errorbar(x_our, surface_time_our, yerr=surface_error_our, fmt='none', ecolor='black', capsize=3)
axs[0].set_xticks(x)
axs[0].set_xticklabels([f'{p:.5f}' for p in p_values], rotation=45)
axs[0].set_ylim(0, 600)  # Adjusted for max time (515 s)
axs[0].set_xlabel('Physical Error Rate $p$')
axs[0].set_ylabel('Running Time (s)')
axs[0].set_title('Surface Code')
axs[0].legend()
axs[0].grid(True, axis='y', ls="--")

# Plot for Repetition code
axs[1].bar(x_our, repetition_time_our, width=bar_width, label='Repetition (Our)', color='#2ca02c', edgecolor='black')
axs[1].bar([x_stim[i] for i in range(len(repetition_time_stim)) if repetition_time_stim[i] is not None], 
           [time for time in repetition_time_stim if time is not None], 
           width=bar_width, label='Repetition (STIM)', color='#98df8a', hatch='//', edgecolor='black')
axs[1].errorbar(x_our, repetition_time_our, yerr=repetition_error_our, fmt='none', ecolor='black', capsize=3)
axs[1].errorbar([x_stim[i] for i in range(len(repetition_time_stim)) if repetition_time_stim[i] is not None], 
                [time for time in repetition_time_stim if time is not None], 
                yerr=[err for err in repetition_error_stim if err is not None], fmt='none', ecolor='black', capsize=3)
axs[1].set_xticks(x)
axs[1].set_xticklabels([f'{p:.5f}' for p in p_values], rotation=45)
axs[1].set_ylim(0, 300)  # Adjusted for max time (269 s)
axs[1].set_xlabel('Physical Error Rate $p$')
axs[1].set_ylabel('Running Time (s)')
axs[1].set_title('Repetition Code')
axs[1].legend()
axs[1].grid(True, axis='y', ls="--")

# Plot for Square code
axs[2].bar(x_stim, square_time_stim, width=bar_width, label='Square (STIM)', color='#ff9896', hatch='//', edgecolor='black')
axs[2].bar(x_our, square_time_our, width=bar_width, label='Square (Our)', color='#d62728', edgecolor='black')
axs[2].errorbar(x_stim, square_time_stim, yerr=square_error_stim, fmt='none', ecolor='black', capsize=3)
axs[2].errorbar(x_our, square_time_our, yerr=square_error_our, fmt='none', ecolor='black', capsize=3)
axs[2].set_xticks(x)
axs[2].set_xticklabels([f'{p:.5f}' for p in p_values], rotation=45)
axs[2].set_ylim(0, 300)  # Adjusted for max time (290 s)
axs[2].set_xlabel('Physical Error Rate $p$')
axs[2].set_ylabel('Running Time (s)')
axs[2].set_title('Square Code')
axs[2].legend()
axs[2].grid(True, axis='y', ls="--")

# Plot for Hexagon code
axs[3].bar(x_stim, hexagon_time_stim, width=bar_width, label='Hexagon (STIM)', color='#c5b0d5', hatch='//', edgecolor='black')
axs[3].bar(x_our, hexagon_time_our, width=bar_width, label='Hexagon (Our)', color='#9467bd', edgecolor='black')
axs[3].errorbar(x_stim, hexagon_time_stim, yerr=hexagon_error_stim, fmt='none', ecolor='black', capsize=3)
axs[3].errorbar(x_our, hexagon_time_our, yerr=hexagon_error_our, fmt='none', ecolor='black', capsize=3)
axs[3].set_xticks(x)
axs[3].set_xticklabels([f'{p:.5f}' for p in p_values], rotation=45)
axs[3].set_ylim(0, 600)  # Adjusted for max time (494 s)
axs[3].set_xlabel('Physical Error Rate $p$')
axs[3].set_ylabel('Running Time (s)')
axs[3].set_title('Hexagon Code')
axs[3].legend()
axs[3].grid(True, axis='y', ls="--")

# Adjust layout to prevent overlap
plt.tight_layout()

# Display the plots
plt.show()