from Bio.Seq import Seq
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog, messagebox
import csv
from datetime import datetime
import threading
import time

def get_user_info():
    name = simpledialog.askstring("Input", "Please enter your name:")
    if not name:
        messagebox.showerror("Error", "Name cannot be empty.")
        return None, None
    sequence = simpledialog.askstring("Input", "Please enter a DNA sequence:").upper()
    if sequence and all(nuc in "ATGC" for nuc in sequence):
        return name, Seq(sequence)
    else:
        messagebox.showerror("Error","Invalid DNA sequence. Please use only A, T, G, and C.")
        return None, None

def calculate_gc_content(sequence):
    return 100 * float(sequence.count("G") + sequence.count("C")) / len(sequence)

def calculate_melting_temperature(sequence):
    num_A = sequence.count("A")
    num_T = sequence.count("T")
    num_G = sequence.count("G")
    num_C = sequence.count("C")
    return 2 * (num_A + num_T) + 4 * (num_G + num_C)

def find_motif(sequence, motif):
    return [i for i in range(len(sequence)) if sequence.startswith(motif, i)]

def plot_combined_results(sequence, gc_content, tm, motif_positions):
    frequencies = {nuc: sequence.count(nuc) for nuc in "ATGC"}
    
    plt.figure(figsize=(20, 5))

    # Nucleotide Frequency Bar Chart
    plt.subplot(1, 4, 1)
    plt.bar(frequencies.keys(), frequencies.values(), color=['red', 'blue', 'green', 'yellow'])
    plt.xlabel("Nucleotides")
    plt.ylabel("Frequency")
    plt.title("Nucleotide Frequency in DNA Sequence")

    # GC Content Bar Chart
    plt.subplot(1, 4, 2)
    plt.bar(['GC Content'], [gc_content], color='green')
    plt.ylabel('Percentage (%)')
    plt.title('GC Content(gluccose content)')
    plt.ylim(0, 100)

    # Melting Temperature (Tm)
    plt.subplot(1, 4, 3)
    plt.axhline(y=tm, color='blue', linestyle='-')
    plt.title('Melting Temperature (Tm)')
    plt.ylabel('Temperature (°C)')
    plt.xticks([])
    plt.yticks(range(0, int(tm) + 20, 5))
    plt.grid(axis='y')

    # Nucleotide Frequency Pie Chart
    plt.subplot(1, 4, 4)
    plt.pie(frequencies.values(), labels=frequencies.keys(), autopct='%1.1f%%', colors=['red', 'blue', 'green', 'yellow'])
    plt.title('Nucleotide Composition')

    plt.tight_layout()
    plt.show()

def save_results(name, sequence, complementary_seq, reverse_complementary_seq, gc_content, tm, motif_positions, binary_representation, hashed_sequence, memory_addresses, compressed_sequence, parallel_chunks):
    with open("dna_analysis_results.csv", "a", newline='') as file:
        writer = csv.writer(file)
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        writer.writerow([f"Analysis Timestamp: {timestamp}"])
        writer.writerow(["Name", name])
        writer.writerow(["Original DNA sequence", str(sequence)])
        writer.writerow(["Complementary sequence", str(complementary_seq)])
        writer.writerow(["Reverse Complementary sequence", str(reverse_complementary_seq)])
        writer.writerow(["GC content", f"{gc_content:.2f}%"])
        writer.writerow(["Length of sequence", len(sequence)])
        writer.writerow(["Melting Temperature (Tm)", tm])
        writer.writerow(["Motif Positions", ', '.join(map(str, motif_positions)) if motif_positions else "Not found"])
        writer.writerow(["Binary Representation", binary_representation])
        writer.writerow(["Hashed Sequence", hashed_sequence])

        writer.writerow(["Memory Addresses", memory_addresses])
        writer.writerow(["Compressed Sequence", compressed_sequence])
        writer.writerow(["Parallel Processing Chunks", parallel_chunks])

        writer.writerow([])  # Blank row for separation

    # Print results to terminal
    print("\n--- DNA Sequence Analysis Results ---")
    print(f"Analysis Timestamp: {timestamp}")
    print(f"Name: {name}")
    print(f"Original DNA sequence: {sequence}")
    print(f"Complementary sequence: {complementary_seq}")
    print(f"Reverse Complementary sequence: {reverse_complementary_seq}")
    print(f"GC content(glucose content): {gc_content:.2f}%")
    print(f"Length of sequence: {len(sequence)}")
    print(f"Melting Temperature (Tm): {tm}")
    print(f"Motif Positions: {', '.join(map(str, motif_positions)) if motif_positions else 'Not found'}")
    print(f"Binary Representation: {binary_representation}")
    print(f"Hashed Sequence: {hashed_sequence}")
    print(f"Memory Addresses: {memory_addresses}")
    print(f"Compressed Sequence: {compressed_sequence}")
    print(f"Parallel Processing Chunks: {parallel_chunks}\n")

def dna_to_binary(sequence):
    binary_mapping = {'A': '00', 'T': '01', 'G': '10', 'C': '11'}
    return ''.join(binary_mapping[nuc] for nuc in sequence)

def hash_dna_sequence(sequence):
    return hash(str(sequence))

def assign_memory_addresses(sequence):
    return {f'Address {i}': nuc for i, nuc in enumerate(sequence)}

def compress_dna_sequence(sequence):
    compressed = []
    count = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            count += 1
        else:
            compressed.append((sequence[i - 1], count))
            count = 1
    compressed.append((sequence[-1], count))
    return compressed

def simulate_parallel_processing(sequence):
    chunk_size = len(sequence) // 4
    chunks = [sequence[i:i + chunk_size] for i in range(0, len(sequence), chunk_size)]
    return chunks

def simulate_pipeline(sequence):
    """Simulate a simple pipelining process by running stages sequentially."""
    print("\n\n")
    print("Pipeline stage 1: GC Content Calculation")
    time.sleep(1)
    gc_content = calculate_gc_content(sequence)
    
    print("Pipeline stage 2: Melting Temperature Calculation")
    time.sleep(1)
    tm = calculate_melting_temperature(sequence)
    
    print("Pipeline stage 3: Binary Conversion")
    time.sleep(1)
    binary_representation = dna_to_binary(sequence)
    
    return gc_content, tm, binary_representation

def run_multithreaded_analysis(sequence):
    """Simulate multithreading by running DNA calculations in separate threads and storing results."""
    results = {}

    def calculate_gc():
        results["gc_content"] = calculate_gc_content(sequence)

    def calculate_tm():
        results["tm"] = calculate_melting_temperature(sequence)
    
    # Start threads
    gc_thread = threading.Thread(target=calculate_gc)
    tm_thread = threading.Thread(target=calculate_tm)
    gc_thread.start()
    tm_thread.start()
    
    # Wait for threads to finish
    gc_thread.join()
    tm_thread.join()
    
    return results["gc_content"], results["tm"]

def main():
    root = tk.Tk()
    root.withdraw()
    
    name, dna_sequence = get_user_info()
    if not name or dna_sequence is None:
        return

    # Core calculations
    complementary_seq = dna_sequence.complement()
    reverse_complementary_seq = dna_sequence.reverse_complement()
    motif = simpledialog.askstring("Input", "Enter a motif to search for (or leave blank to skip):")
    motif_positions = find_motif(dna_sequence, motif) if motif else []

    # Simulate pipeline processing
    gc_content, tm, binary_representation = simulate_pipeline(dna_sequence)
    
    # Simulate multithreading for additional functions
    gc_content_threaded, tm_threaded = run_multithreaded_analysis(dna_sequence)
    
    # New attributes for analysis
    hashed_sequence = hash_dna_sequence(dna_sequence)
    memory_addresses = assign_memory_addresses(dna_sequence)
    compressed_sequence = compress_dna_sequence(dna_sequence)
    parallel_chunks = simulate_parallel_processing(dna_sequence)

    # Save results, print, and plot
    save_results(name, dna_sequence, complementary_seq, reverse_complementary_seq, gc_content, tm, motif_positions, binary_representation, hashed_sequence, memory_addresses, compressed_sequence, parallel_chunks)
    plot_combined_results(dna_sequence, gc_content, tm, motif_positions)

if __name__ == "__main__":
    main()