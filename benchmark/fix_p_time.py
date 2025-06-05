import matplotlib.pyplot as plt
import numpy as np

# Define physical error rates p
p_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 5e-05]

# Running times and errors
surface_time_stim = [0.033, 0.042, 0.013, 0.051, 1.19, 6.96, 503, 515]
surface_error_stim = [0.007, 0.005, 0.004, 0.035, 0.98, 0.25, 39, 148]
surface_time_ours = [2.85, 2.37, 2.16, 2.13, 3.01, 2.89, 2.87, 3.15]
surface_error_ours = [0.08, 0.01, 0.11, 0.07, 0.31, 0.37, 0.36, 0.34]

repetition_time_stim = [0.0045, 0.008, 0.099, 0.56, 32.9, 269, None, None]
repetition_error_stim = [0.0004, 0.003, 0.014, 0.13, 4.4, 40, None, None]
repetition_time_ours = [0.735, 0.711, 0.797, 0.781, 0.786, 0.822, 0.789, 0.787]
repetition_error_ours = [0.09, 0.16, 0.25, 0.24, 0.19, 0.89, 0.25, 0.28]

square_time_stim = [0.06, 0.041, 0.015, 0.031, 0.68, 5.26, 290, 200]
square_error_stim = [0.01, 0.005, 0.002, 0.007, 0.21, 0.61, 19, 68]
square_time_ours = [4.34, 3.45, 2.73, 2.76, 3.31, 4.89, 4.82, 4.40]
square_error_ours = [0.26, 0.07, 0.04, 0.06, 0.10, 0.47, 0.28, 0.58]

hexagon_time_stim = [0.11, 0.071, 0.024, 0.03, 0.65, 4.84, 421, 494]
hexagon_error_stim = [0.019, 0.012, 0.005, 0.004, 0.23, 1.12, 17, 106]
hexagon_time_ours = [7.00, 5.77, 4.95, 4.23, 4.68, 5.92, 8.99, 9.29]
hexagon_error_ours = [0.05, 0.06, 1.23, 0.14, 0.12, 0.53, 2.37, 1.86]

# Create a single plot
fig, ax = plt.subplots(figsize=(12, 8))

# Bar width and positions
num_bars = 8  # 4 codes * 2 methods (STIM and Ours)
offsets = [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]  # Offsets for 8 bars
x_positions = [[p * (1 + offset * 0.05) for p in p_values] for offset in offsets]

# Plot bars and error bars
ax.bar(x_positions[0], surface_time_stim, width=[p * 0.05 for p in p_values], label='Surface (STIM)', color='#aec7e8', hatch='//', edgecolor='black')
ax.bar(x_positions[1], surface_time_ours, width=[p * 0.05 for p in p_values], label='Surface (Ours)', color='#1f77b4', edgecolor='black')
# Filter None values for Repetition (STIM)
valid_indices = [i for i, time in enumerate(repetition_time_stim) if time is not None]
ax.bar([x_positions[2][i] for i in valid_indices], 
       [repetition_time_stim[i] for i in valid_indices], 
       width=[p * 0.05 for p in p_values][:len(valid_indices)], 
       label='Repetition (STIM)', color='#98df8a', hatch='//', edgecolor='black')
ax.bar(x_positions[3], repetition_time_ours, width=[p * 0.05 for p in p_values], label='Repetition (Ours)', color='#2ca02c', edgecolor='black')
ax.bar(x_positions[4], square_time_stim, width=[p * 0.05 for p in p_values], label='Square (STIM)', color='#ff9896', hatch='//', edgecolor='black')
ax.bar(x_positions[5], square_time_ours, width=[p * 0.05 for p in p_values], label='Square (Ours)', color='#d62728', edgecolor='black')
ax.bar(x_positions[6], hexagon_time_stim, width=[p * 0.05 for p in p_values], label='Hexagon (STIM)', color='#c5b0d5', hatch='//', edgecolor='black')
ax.bar(x_positions[7], hexagon_time_ours, width=[p * 0.05 for p in p_values], label='Hexagon (Ours)', color='#9467bd', edgecolor='black')

# Error bars
ax.errorbar(x_positions[0], surface_time_stim, yerr=surface_error_stim, fmt='none', ecolor='black', capsize=3)
ax.errorbar(x_positions[1], surface_time_ours, yerr=surface_error_ours, fmt='none', ecolor='black', capsize=3)
ax.errorbar([x_positions[2][i] for i in valid_indices], 
            [repetition_time_stim[i] for i in valid_indices], 
            yerr=[repetition_error_stim[i] for i in valid_indices], fmt='none', ecolor='black', capsize=3)
ax.errorbar(x_positions[3], repetition_time_ours, yerr=repetition_error_ours, fmt='none', ecolor='black', capsize=3)
ax.errorbar(x_positions[4], square_time_stim, yerr=square_error_stim, fmt='none', ecolor='black', capsize=3)
ax.errorbar(x_positions[5], square_time_ours, yerr=square_error_ours, fmt='none', ecolor='black', capsize=3)
ax.errorbar(x_positions[6], hexagon_time_stim, yerr=hexagon_error_stim, fmt='none', ecolor='black', capsize=3)
ax.errorbar(x_positions[7], hexagon_time_ours, yerr=hexagon_error_ours, fmt='none', ecolor='black', capsize=3)

# Set scales and labels
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-3, 1e3)
ax.set_xticks(p_values, [f'{p:.5f}' for p in p_values], rotation=45)
ax.set_xlabel('Physical Error Rate $p$')
ax.set_ylabel('Running Time (s)')
ax.set_title('Running Time Comparison for All Codes')
ax.legend()
ax.grid(True, which="both", ls="--", axis='y')

# Adjust layout
plt.tight_layout()

# Display the plot
plt.show()