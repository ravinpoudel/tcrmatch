import os
import subprocess
import argparse
import pandas as pd
from Bio import SeqIO

def check_docker_image(image_name):
    try:
        # Check if the Docker image exists locally
        subprocess.run(['docker', 'image', 'inspect', image_name], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False

def build_docker_image(image_name):
    # Build the Docker image
    subprocess.run(['docker', 'build', '-t', image_name, '.'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def run_tcrmatch(tcrbetaseq_file):
    # Validate the input file
    if not os.path.isfile(tcrbetaseq_file):
        raise FileNotFoundError(f"File not found: {tcrbetaseq_file}")

    # Extract basename from the input_file_path
    input_file_basename = os.path.basename(tcrbetaseq_file)

    # Define the Docker image name
    docker_image_name = 'tcrmatch'

    # Check if the Docker image exists
    if not check_docker_image(docker_image_name):
        print(f"Docker image '{docker_image_name}' not found. Building the image...")
        build_docker_image(docker_image_name)

    # Define the Docker command
    docker_command = [
        'docker', 'run',
        '-v', f'{os.path.abspath(os.path.dirname(tcrbetaseq_file))}:/app/inputdata',
        docker_image_name,
        './tcrmatch',
        '-i', f'inputdata/{input_file_basename}',  # Use the extracted basename
        '-t', '4',
        '-d', 'data/IEDB_data.tsv',
        '-a'
    ]

    # Run the Docker command
    subprocess.run(docker_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Optionally, you can capture the output to a file
    output_file_path = f'output_{input_file_basename}.tsv'
    with open(output_file_path, 'w') as output_file:
        subprocess.run(docker_command, stdout=output_file, stderr=subprocess.PIPE)
    print(f"Output saved to {output_file_path}")
    return output_file_path


def get_unique_epitopes_from_tcrmatch(input_file):
    df = pd.read_csv(input_file, sep='\t')
    epitope_column = list(set(df['epitope'].tolist()))
    sublists = [element.split(',') if ',' in element else [element] for element in epitope_column]
    unique_epitopes = list(set([item for sublist in sublists for item in sublist]))
    return unique_epitopes

def generate_kmers(sequence, min_length=8, max_length=11):
    kmers = []
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            kmer = sequence[i:i + length]
            kmers.append(kmer)
    return kmers


def run_prodigal(input_genome, output_genes, output_proteins):
    # Define the prodigal command
    prodigal_command = [
        'pyrodigal',
        '-i', input_genome,
        '-o', output_genes,
        '-a', output_proteins
    ]
    # Run the prodigal command
    subprocess.run(prodigal_command)


def generate_kmers_from_genome_protein(proteins_file, min_kmer_length, max_kmer_length, output_kmers_file):
    protein_sequences = []
    # Use SeqIO to read the output_proteins_file
    for record in SeqIO.parse(proteins_file, "fasta"):
        protein_sequences.append(str(record.seq))
    # Generate kmers for each protein sequence
    all_kmers = []
    for protein_sequence in protein_sequences:
        kmers = generate_kmers(protein_sequence, min_kmer_length, max_kmer_length)
        all_kmers.extend(kmers)
    # Save kmers to the specified output file
    unique_all_kmers = list(set(all_kmers))
    with open(output_kmers_file, 'w') as file:
        for kmer in unique_all_kmers:
            file.write(f"{kmer}\n")
    print(f"Total number of elements in all_kmers: {len(all_kmers)}")
    print(f"Total number of elements in unique_all_kmers: {len(unique_all_kmers)}")
    print(f"All unique kmers saved to {output_kmers_file}")
    return unique_all_kmers

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run tcrmatch with Docker")
    parser.add_argument("-tcrbetaseq_file", help="Path to the tcrbeta sequence file")
    parser.add_argument("-input_genome", help="Input genome file")
    parser.add_argument("-mismatch", help="Number of mismatch allow in kmer match")

    args = parser.parse_args()

    output_file_path = run_tcrmatch(args.tcrbetaseq_file)
    
    # Extract basename without extension from the input_file_path
    input_file_basename, _ = os.path.splitext(os.path.basename(args.tcrbetaseq_file))

    print("Input TCR sequnce file: ", input_file_basename)

    unique_epitopes_from_TCRmatch = get_unique_epitopes_from_tcrmatch(output_file_path)

    # generate 8 - 11 bp kmers from tcrmatch results
    kmers_8_11bp = [generate_kmers(x) for x in unique_epitopes_from_TCRmatch]
    all_kmers = list(set([item for sublist in kmers_8_11bp for item in sublist]))

    # Save the modified_list to a file
    output_file_path_all_kmers_TCR = f'{input_file_basename}_all_kmers_from_tcrmatch.txt'
    with open(output_file_path_all_kmers_TCR, 'w') as output_file:
        for item in all_kmers:
            output_file.write(f"{item}\n")

    # Run prodigal
    input_genome_basename, _ = os.path.splitext(os.path.basename(args.input_genome))
    print("Input genome: ",input_genome_basename)
    output_genes_file = f'{input_genome_basename}.genes'
    output_proteins_file = f'{input_genome_basename}.faa'
    run_prodigal(args.input_genome, output_genes_file, output_proteins_file)

    # Generate kmers from genome protein
    min_length_param = 8
    max_length_param = 11
    output_kmers_file_for_genome = f'{input_genome_basename}_all_kmers.txt'
    kmers_list = generate_kmers_from_genome_protein(output_proteins_file, min_length_param, max_length_param, output_kmers_file_for_genome)

    # running tcr match
    kmers_list_TCRmatch = output_file_path_all_kmers_TCR
    kmers_list_genome = output_kmers_file_for_genome
    cross_match_output = f'{input_genome_basename}_cross_match.txt'
    cross_match_summary_output = f'{input_genome_basename}_cross_match_summary.txt'
    mismatch = args.mismatch
    cross_match_only_crags_output = f'{input_genome_basename}_cross_match_only_crags.txt'

    # Run the R script using subprocess
    subprocess.run(['Rscript', 'kmer_match_detailed_cl_many_kmers.R', kmers_list_TCRmatch, kmers_list_genome,
                    cross_match_output, cross_match_summary_output, mismatch,
                    cross_match_only_crags_output])
        