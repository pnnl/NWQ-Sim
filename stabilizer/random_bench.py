import os
import numpy as np
import matplotlib.pyplot as plt

def read_and_label_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    labels = ["sim_name", "sim_time", "distance", "rounds", "n_qubits"]
    data = {labels[i]: lines[i].strip() if i == 0 else float(lines[i]) for i in range(len(lines))}
    
    return data

def process_folder(folder_path):
    all_data = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"): 
            filepath = os.path.join(folder_path, filename)
            data = read_and_label_file(filepath)
            all_data.append(data)
    
    return sorted(all_data, key=lambda x: x["distance"])  # Sort by sim_time

def plot_data(all_data):
    sim_times_stab = []
    sim_times_stim = []
    distance_stab = []
    distance_stim = []
    sim_times_cpu = []
    distance_cpu = []
    sim_times_qiskit= []
    distance_qiskit = []

    for data in all_data:
        if data["sim_name"] == "nvgpu" and data["sim_time"] < 620:
            distance_stab.append(float(data["distance"]))
            sim_times_stab.append(float(data["sim_time"]))
        elif data["sim_name"] == "stim" and data["sim_time"] < 620:
            distance_stim.append(float(data["distance"]))
            sim_times_stim.append(float(data["sim_time"]))
        # elif data["sim_name"] == "cpu" and data["sim_time"] < 1300:
        #     distance_cpu.append(float(data["distance"]))
        #     sim_times_cpu.append(float(data["sim_time"]))
        elif data["sim_name"] == "qiskit_random" and data["sim_time"] < 620:
            distance_qiskit.append(float(data["distance"]))
            sim_times_qiskit.append(float(data["sim_time"]))

    plt.figure(figsize=(7, 5))
    plt.rcParams.update({
    'font.size': 18,          # Default font size
    # 'axes.titlesize': 20,     # Title font size
    # 'axes.labelsize': 18,     # Axis label font size
    # 'xtick.labelsize': 16,    # X-axis tick label size
    # 'ytick.labelsize': 16,    # Y-axis tick label size
    # 'legend.fontsize': 16,    # Legend font size
    })

    # plt.plot(distance_cpu, sim_times_cpu, "^-", label="STAB-CPU", linestyle="-", color="orange", markersize=5)
    plt.plot(distance_qiskit, sim_times_qiskit, "<-", label="Qiskit ('stabilizer')", linestyle="--", color="deepskyblue", markersize=6)
    plt.plot(distance_stim, sim_times_stim, ">-", label="Stim", linestyle="--", color="green", markersize=6)

    plt.plot(distance_stab, sim_times_stab, "v-", label="STABSim", linestyle="-", color="red", markersize=6)

    # Find where STAB-GPU overtakes Stim
    # crossover_distance = None
    # if len(distance_stab) > 5 and len(distance_stim) > 5:
    #     # Find the crossover point by comparing interpolated values
    #     min_dist = max(min(distance_stab), min(distance_stim))
    #     max_dist = min(max(distance_stab), max(distance_stim))
        
    #     # Check distances in the overlapping range
    #     for d in range(int(min_dist), int(max_dist) + 1):
    #         # Find closest points for interpolation
    #         stab_times_at_d = []
    #         stim_times_at_d = []
            
    #         for i, dist in enumerate(distance_stab):
    #             if abs(dist - d) < 2:  # Within 2 units
    #                 stab_times_at_d.append(sim_times_stab[i])
            
    #         for i, dist in enumerate(distance_stim):
    #             if abs(dist - d) < 2:  # Within 2 units
    #                 stim_times_at_d.append(sim_times_stim[i])
            
    #         if stab_times_at_d and stim_times_at_d:
    #             avg_stab = np.mean(stab_times_at_d)
    #             avg_stim = np.mean(stim_times_at_d)
                
    #             if avg_stab < avg_stim and crossover_distance is None:
    #                 crossover_distance = d
    #                 break
    
    # Add vertical line at crossover point
    # if crossover_distance:
    #     plt.axvline(x=crossover_distance, color='black', linestyle='--', alpha=0.7, linewidth=2)

    # Increase font size for axis labels and title
    plt.xlabel(r"# of Qubits = # of Gates")#, fontsize=22)
    plt.ylabel("Sim Time (s)")#, fontsize=22)
    # plt.title("Random Wide Circuit Simulation", fontsize=22)
    
    # Increase font size for ticks
    # plt.xticks(fontsize=22)
    # plt.yticks(fontsize=22)

    # Set the x-axis to log2 scale
    plt.ylim(0, 630)
    plt.xlim(5000, 31500)
    # plt.yscale("log", base=10)
    # plt.xscale("log", base=2)

    
    # Only show the x-axis ticks for 5, 15, 51, and 99 on the log2 scale
    # plt.xticks([5, 15, 27, 50, 100, 125], labels=[str(x) for x in [5, 15, 27, 50, 100, 125]])

    # Increase font size for the legend
    plt.legend(loc='upper left')

    plt.grid(True, linestyle=":", linewidth=0.7, which="both")

    # Adjust layout to minimize whitespace and prevent cutoff
    plt.tight_layout()

    # Save the plot as a PDF, ensuring all elements fit
    plt.savefig("/people/garn195/NWQ-Sim/stabilizer/graphics/random_improved.pdf", format="pdf", bbox_inches='tight', dpi=1000)

    plt.show()

folder_path = "/people/garn195/NWQ-Sim/stabilizer/random_bench_new"  
all_data = process_folder(folder_path)

print(f"Found {len(all_data)} files")

if not all_data:
    print("No data found.")

plot_data(all_data)
