# HaploCluster

## Under Development

A Python tool for **haplotype clustering** and **noise reassignment** from DNA sequencing reads in FASTQ files. HaploCluster processes sequencing reads using k-mer profiles, applies **DBSCAN clustering**, and ensures that all reads—including noise—are assigned to meaningful clusters. The output consists of clustered reads saved into individual FASTQ files.

## Key Features
- **DBSCAN Clustering**: Groups similar DNA sequences using Jaccard distance and user-defined parameters (`k` and `eps`).
- **Noise Reassignment**: Reads identified as noise by DBSCAN are reassigned to the nearest cluster based on distance to cluster centroids.
- **Flexible K-mer Analysis**: Extracts k-mer profiles with a default k-mer size of `k=8` to represent sequence features effectively.
- **Output FASTQ Files**: Saves each cluster's reads into separate FASTQ files with descriptive filenames, including cluster ID and the number of reads in each cluster.

## How It Works

1. **Input**:
   - HaploCluster scans the working directory for FASTQ files and processes them one by one.

2. **K-mer Profile Creation**:
   - Each sequence is represented as a k-mer profile to capture sequence-level features.

3. **Jaccard Distance Matrix**:
   - Pairwise distances between sequences are calculated using the Jaccard similarity of k-mer profiles.

4. **DBSCAN Clustering**:
   - Sequences are clustered using DBSCAN with customizable parameters (`eps=0.15` by default).

5. **Noise Reassignment**:
   - Reads identified as noise are reassigned to the nearest cluster based on the distance to cluster centroids.

6. **Output**:
   - Each cluster's reads are saved to a separate FASTQ file named in the format:
     ```
     Cluster_XX-YYseqs-OriginalFileName.fastq
     ```
     - `XX`: Cluster ID
     - `YY`: Number of reads in the cluster
     - `OriginalFileName`: Name of the input FASTQ file.

## Usage

### Requirements
- Python 3.8+
- Install dependencies using `pip`:
  ```bash
  pip install biopython numpy scikit-learn

## Running the Script

1. Place your FASTQ files in the working directory.
2. Run the script:
   ```bash
   python cluster.py
3. The script will process each FASTQ file in the directory and generate clustered FASTQ files for each cluster.

## Output Example

For an input file `reads.fastq`:

### Output FASTQ Files:
- `Cluster_1-120seqs-reads.fastq`
- `Cluster_2-80seqs-reads.fastq`
- `...`

### Console Output:
```plaintext
Processing file: reads.fastq
Cluster 1: 120 reads saved to Cluster_1-120seqs-reads.fastq
Cluster 2: 80 reads saved to Cluster_2-80seqs-reads.fastq
 ```
## Customization

### Adjust Parameters

Modify the `k` value (k-mer size) and `eps` value (DBSCAN clustering tightness) directly in the script:

```python
cluster_and_save_reads(input_fastq, k=8, eps=0.15)
 ```
## Applications

- **Haplotype Clustering**: Organize sequencing reads into distinct haplotype clusters based on sequence similarity.
- **Noise Reduction**: Eliminate sequencing noise by reassigning noisy reads to appropriate clusters.
- **Environmental Sequencing**: Cluster reads from complex microbial or fungal communities.

## Future Enhancements

- Parameter optimization for automatic `k` and `eps` selection.
- Support for additional distance metrics (e.g., edit distance).
- Visualization of clustering results.

