import streamlit as st
from Bio import SeqIO
import requests
import time
from io import StringIO
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt


# Function to calculate molecular weight of a protein sequence
def calculate_molecular_weight(protein_sequence):
    # Placeholder implementation - replace with your actual calculation logic
    # Here, we're assuming a fixed weight for each amino acid
    amino_acid_weights = {'A': 71.07, 'R': 156.18, 'N': 114.08, 'D': 115.08, 'C': 103.10,
    'Q': 128.13, 'E': 129.11, 'G': 57.05, 'H': 137.14, 'I': 113.15,
    'L': 113.15, 'K': 128.17, 'M': 131.19, 'F': 147.17, 'P': 97.11,
    'S': 87.08, 'T': 101.10, 'W': 186.20, 'Y': 163.17, 'V': 99.13}
    weight = sum(amino_acid_weights.get(aa, 0) for aa in protein_sequence)
    return weight

# Function to retrieve protein data (sequence and gene) with retry
def retrieve_protein_data_with_retry(identifier, max_retries=3):
    url = f"https://www.uniprot.org/uniprot/{identifier}.fasta"
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            response.raise_for_status()
            # Create a StringIO object from the response text
            fasta_text = StringIO(response.text)
            # Use SeqIO.read with the StringIO object
            protein_record = SeqIO.read(fasta_text, "fasta")
            # Retrieve gene information from UniProt
            gene_url = f"https://www.uniprot.org/uniprot/{identifier}.txt"
            gene_response = requests.get(gene_url)
            gene_response.raise_for_status()
            gene_data = gene_response.text
            gene_info = ""
            for line in gene_data.splitlines():
                if line.startswith("GN"):
                    gene_info += line + "\n"
            return protein_record, gene_info
        except requests.RequestException as e:
            if attempt < max_retries - 1:
                delay = 2 ** attempt
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                raise e

# Set page title and icon
st.set_page_config(page_title="Protein Data Explorer", page_icon=":microscope:")

# Main title and description
st.title("Protein Data Explorer")
st.write("This web app allows you to explore protein data and interaction networks.")

# Input section
input_option = st.radio("Select Input Option:", ("Uniprot ID", "Protein Sequence"))

uniprot_id = None  # Define uniprot_id variable outside the conditionals

if input_option == "Uniprot ID":
    uniprot_id = st.text_input("Enter Uniprot ID:")
else:
    protein_sequence = st.text_area("Enter Protein Sequence (FASTA format):", height=200)

# Define protein_record and gene_info variables outside the input conditionals
protein_record = None
gene_info = ""

if input_option == "Uniprot ID":
    if uniprot_id:
        protein_record, gene_info = retrieve_protein_data_with_retry(uniprot_id)
        if protein_record:
            st.write("Protein Sequence:")
            st.code(str(protein_record.seq))
            st.write("Gene Information:")
            st.code(gene_info)
        else:
            st.error("Failed to retrieve protein data. Please check your Uniprot ID.")
else:
    if protein_sequence:
         # Create a StringIO object from the text area input
        fasta_text = StringIO(protein_sequence)
        protein_record = SeqIO.read(fasta_text, "fasta")
        st.write("Protein Sequence:")
        st.code(protein_sequence)

# Display protein characteristics and fetch protein-protein interaction network only if protein_record is not None
if protein_record:
    st.write("Protein Characteristics:")
    st.write(f"- Length: {len(protein_record.seq)} amino acids")
    st.write(f"- Molecular Weight: {calculate_molecular_weight(protein_record.seq):,.0f} Daltons")
    st.write(f"- Gene Information:")
    st.code(gene_info)

    # Add more characteristics as needed

    # Fetch protein-protein interaction network
    def fetch_interaction_network(uniprot_id):
        url = f"https://string-db.org/api/tsv/network?identifiers={uniprot_id}"
        response = requests.get(url)
        if response.ok:
            network_data = response.text
            network_df = pd.read_csv(StringIO(network_data), sep='\t')
            return network_df
        else:
            return None

    st.write("Protein-Protein Interaction Network:")
    interaction_network = fetch_interaction_network(uniprot_id)
    if interaction_network is not None:
        st.write(interaction_network)
        ppi_graph = nx.from_pandas_edgelist(interaction_network, "preferredName_A", "preferredName_B")
        pos = nx.spring_layout(ppi_graph)  # Define node positions
        fig, ax = plt.subplots()  # Create a new figure
        nx.draw(ppi_graph, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=8, ax=ax)
        st.pyplot(fig)  # Display the graph in Streamlit with the specified figure
    else:
        st.error("Failed to fetch interaction network. Please try again later.")


# Footer
st.sidebar.markdown("Bioinformatics II - Assignment")
st.sidebar.markdown("---")
st.sidebar.markdown("Jeliza Justine A/P Sebastin - A21EC0034")
st.sidebar.markdown("Harchana Arulappan - A21EC0028")
st.sidebar.markdown("Siti Nurkamilah Binti Saiful Bahari - A21EC0131")
