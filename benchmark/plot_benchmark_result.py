import numpy as np
import matplotlib.pyplot as plt





def plot_1million_shots_different_d():

    # Distances for which time is measured
    distances = [3, 5, 7, 9, 11]

    # Time taken (seconds) for My code vs STIM at each distance, per code type
    # Replace these with your actual timing results
    my_time = {
        "Repetition": [0.12, 0.30, 0.70, 1.50, 3.40],
        "Surface": [0.10, 0.25, 0.60, 1.30, 2.90],
        "Square": [0.14, 0.32, 0.75, 1.60, 3.60],
        "Hexagon": [0.11, 0.28, 0.65, 1.40, 3.10],
    }

    stim_time = {
        "Repetition": [0.45, 2.70, 7.00, 13.60, 28.40],
        "Surface": [0.42, 2.50, 6.50, 13.00, 27.00],
        "Square": [0.48, 2.90, 7.20, 14.40, 30.00],
        "Hexagon": [0.43, 2.60, 6.80, 13.60, 28.00],
    }

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharey=True)

    for ax, ct in zip(axes.flat, my_time.keys()):
        x = [0, 1]  # 0 = My code, 1 = STIM
        for d, m, s in zip(distances, my_time[ct], stim_time[ct]):
            ax.plot(x, [m, s], marker='o', label=f"d={d}")
        ax.set_xticks(x)
        ax.set_xticklabels(["My code", "STIM"])
        ax.set_title(f"{ct} code")
        ax.set_ylabel("Time (s)")
        ax.grid(axis="y", linestyle="--", alpha=0.5)

    # Shared legend
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Distance", 
            bbox_to_anchor=(0.85, 0.5), loc="center left",fontsize=15)
    fig.suptitle("Slope Charts: Time for 1M Samples by Code Type")
    plt.tight_layout(rect=[0,0,0.85,0.95])
    plt.savefig("1millionsample.png", dpi=300, bbox_inches='tight')
    plt.show()










def plot_benchmark_result_different_d():
    # Extended distances
    distances = np.array([3, 5, 7, 9, 11, 13, 15])

    # Your measured sample counts (one list of length 7 per code type)
    my_samples = {
        "Repetition": [120, 300, 700, 1500, 3400, 6800, 13600],
        "Surface": [100, 250, 600, 1300, 2900, 5800, 11600],
        "Square": [140, 320, 750, 1600, 3600, 7200, 14400],
        "Hexagon": [110, 280, 650, 1400, 3100, 6200, 12400],
    }

    # STIM’s sample counts
    stim_samples = {
        "Repetition": [450, 2700, 7000, 13600, 28400, 56800, 113600],
        "Surface": [420, 2500, 6500, 13000, 27000, 54000, 108000],
        "Square": [480, 2900, 7200, 14400, 30000, 60000, 120000],
        "Hexagon": [430, 2600, 6800, 13600, 28000, 56000, 112000],
    }

    width = 0.35
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=True)

    for ax, code_type in zip(axes.flat, my_samples.keys()):
        # --- Bars ---
        ax.bar(distances - width/2, my_samples[code_type], width, label="My method")
        ax.bar(distances + width/2, stim_samples[code_type], width, label="STIM")
        
        # --- Exponential fit for My code: log(y) ~ m*x + c  ---
        y = np.array(my_samples[code_type])
        coeffs = np.polyfit(distances, np.log(y), 3)
        fit_my = np.exp(np.polyval(coeffs, distances))
        ax.plot(distances, fit_my,
                linestyle="--", label="Fit – My method")
        
        # --- Exponential fit for STIM ---
        y2 = np.array(stim_samples[code_type])
        coeffs2 = np.polyfit(distances, np.log(y2), 3)
        fit_stim = np.exp(np.polyval(coeffs2, distances))
        ax.plot(distances, fit_stim,
                linestyle="-.", label="Fit – STIM")
        
        # --- Formatting ---
        ax.set_title(code_type)
        ax.set_xlabel("Code Distance")
        ax.set_ylabel("Samples Required")
        ax.set_xticks(distances)
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        ax.legend()

    plt.suptitle("Sample Requirements & Exponential Fits: My Code vs. STIM", y=1.02)
    plt.tight_layout()
    plt.savefig("sample_requirements.png", dpi=300, bbox_inches='tight')
    plt.show()





def plot_benchmark_result_different_p():

    # p values as categories
    p_values = ["0.1", "0.05", "0.01", "0.005", "0.001", "0.0005", "0.0001"]
    positions = np.arange(len(p_values))

    # Placeholder sample counts at distance=7 for each code type
    # Replace these with your actual experimental results
    my_samples = {
        "Surface": [50, 100, 200, 400, 800, 1600, 3200],
        "Repetition": [40,  90, 180, 360, 720, 1440, 2880],
        "Hexagon": [60, 120, 240, 480, 960, 1920, 3840],
        "Square": [55, 110, 220, 440, 880, 1760, 3520],
    }

    stim_samples = {
        "Surface": [200, 400, 800, 1600, 3200, 6400, 12800],
        "Repetition": [180, 360, 720, 1440, 2880, 5760, 11520],
        "Hexagon": [240, 480, 960, 1920, 3840, 7680, 15360],
        "Square": [220, 440, 880, 1760, 3520, 7040, 14080],
    }

    width = 0.35
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=True)

    for ax, code_type in zip(axes.flat, my_samples.keys()):
        ax.bar(positions - width/2, my_samples[code_type], width, label="My code")
        ax.bar(positions + width/2, stim_samples[code_type], width, label="STIM")
        ax.set_xticks(positions)
        ax.set_xticklabels(p_values)
        ax.set_title(code_type)
        ax.set_xlabel("Physical error rate p")
        ax.set_ylabel("Samples Required")
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        ax.legend()

    plt.suptitle("Samples Needed vs. p at Code Distance 7: My Code vs. STIM", y=1.02)
    plt.tight_layout()
    plt.savefig("sample_requirement_different_p.png", dpi=300, bbox_inches='tight')
    plt.show()



def plot_benchmark_result_different_d_2():
    # Distances you tested
    distances = np.array([3, 5, 7, 9, 11, 13, 15])

    # Sample counts (replace with your real data)
    my_samples = {
        "Surface": [120, 300, 700, 1500, 3400, 6800, 13600],
        "Repetition": [100, 250, 600, 1300, 2900, 5800, 11600],
        "Square": [140, 320, 750, 1600, 3600, 7200, 14400],
        "Hexagon": [110, 280, 650, 1400, 3100, 6200, 12400],
    }

    stim_samples = {
        "Surface": [450, 2700, 7000, 13600, 28400, 56800, 113600],
        "Repetition": [420, 2500, 6500, 13000, 27000, 54000, 108000],
        "Square": [480, 2900, 7200, 14400, 30000, 60000, 120000],
        "Hexagon": [430, 2600, 6800, 13600, 28000, 56000, 112000],
    }

    # Base colors per code type
    base_colors = {
        "Surface": "#1f77b4",
        "Repetition": "#ff7f0e",
        "Square": "#2ca02c",
        "Hexagon": "#d62728",
    }

    n_dist = len(distances)
    code_types = list(my_samples.keys())
    n_series = len(code_types) * 2
    total_width = 0.8
    bar_width = total_width / n_series
    x_base = np.arange(n_dist)

    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot bars only
    for i, ct in enumerate(code_types):
        col = base_colors[ct]
        offset_my = (i*2 - (n_series-1)/2) * bar_width
        offset_stim = ((i*2+1) - (n_series-1)/2) * bar_width
        
        ax.bar(x_base + offset_my, my_samples[ct], width=bar_width, color=col, label=f"Our-{ct}")
        ax.bar(x_base + offset_stim, stim_samples[ct], width=bar_width,
            color=col, alpha=0.5, hatch='//', label=f"Stim-{ct}")

    # Formatting
    ax.set_xticks(x_base)
    ax.set_xticklabels(distances)
    ax.set_xlabel("Code Distance")
    ax.set_ylabel("Samples Required")
    ax.set_title("Samples vs. Code Distance: Our method vs. STIM")
    ax.grid(axis="y", linestyle="--", alpha=0.5)

    # Place legend inside top-left
    ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=18)

    plt.tight_layout()
    plt.savefig("sample_requirement_different_d.png", dpi=300, bbox_inches='tight')
    plt.show()







def plot_time_comparision():
    # Code distances you tested
    distances = np.array([3, 5, 7, 9, 11, 13, 15])
    labels = [str(d) for d in distances]

    # Time taken (in seconds) at each distance for each code type
    # Replace these with your actual timing results
    my_time = {
        "Surface": [0.12, 0.30, 0.70, 1.50, 3.40, 6.80, 13.60],
        "Repetition": [0.10, 0.25, 0.60, 1.30, 2.90, 5.80, 11.60],
        "Square": [0.14, 0.32, 0.75, 1.60, 3.60, 7.20, 14.40],
        "Hexagon": [0.11, 0.28, 0.65, 1.40, 3.10, 6.20, 12.40],
    }

    stim_time = {
        "Surface": [0.45, 2.70, 7.00, 13.60, 28.40, 56.80, 113.60],
        "Repetition": [0.42, 2.50, 6.50, 13.00, 27.00, 54.00, 108.00],
        "Square": [0.48, 2.90, 7.20, 14.40, 30.00, 60.00, 120.00],
        "Hexagon": [0.43, 2.60, 6.80, 13.60, 28.00, 56.00, 112.00],
    }

    # Base colors per code type
    base_colors = {
        "Surface": "#1f77b4",
        "Repetition": "#ff7f0e",
        "Square": "#2ca02c",
        "Hexagon": "#d62728",
    }

    n_dist = len(distances)
    code_types = list(my_time.keys())
    n_series = len(code_types) * 2
    total_height = 0.8
    bar_height = total_height / n_series
    y_base = np.arange(n_dist)

    fig, ax = plt.subplots(figsize=(8, 10))

    # Plot horizontal bars
    for i, ct in enumerate(code_types):
        col = base_colors[ct]
        # "My code" offset
        offset_my = (i*2 - (n_series-1)/2) * bar_height
        # "STIM" offset
        offset_stim = ((i*2+1) - (n_series-1)/2) * bar_height
        
        ax.barh(y_base + offset_my, my_time[ct], height=bar_height,
                color=col, label=f"My-{ct}")
        ax.barh(y_base + offset_stim, stim_time[ct], height=bar_height,
                color=col, alpha=0.5, hatch='//', label=f"S-{ct}")

    # Formatting
    ax.set_yticks(y_base)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Time Taken (seconds)")
    ax.set_ylabel("Code Distance")
    ax.set_title("Horizontal Benchmark: Time vs. Code Distance")
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    # Place legend inside
    ax.legend(loc='lower right', bbox_to_anchor=(0.98, 0.02), fontsize=18)

    plt.tight_layout()
    plt.savefig("time_requirement_different_d.png", dpi=300, bbox_inches='tight')
    plt.show()






if __name__ == "__main__":

    plot_time_comparision()
    #plot_1million_shots_different_d2()
    #plot_benchmark_result_different_p()