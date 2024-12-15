import os
import itertools
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import squareform

def kmer_profile(seq, k=8):
    """Generate a k-mer profile for a given sequence."""
    profile = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        profile[kmer] = profile.get(kmer, 0) + 1
    return profile

def jaccard_distance(profile_a, profile_b):
    """Compute Jaccard distance between two k-mer profiles."""
    set_a = set(profile_a.keys())
    set_b = set(profile_b.keys())
    inter = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))
    if union == 0:
        return 1.0
    return 1 - float(inter) / union

def compute_centroid(profiles):
    """Compute the centroid of a cluster (average k-mer profile)."""
    centroid = {}
    num_profiles = len(profiles)
    for profile in profiles:
        for kmer, count in profile.items():
            centroid[kmer] = centroid.get(kmer, 0) + count / num_profiles
    return centroid

def assign_noise_to_clusters(noise_reads, noise_profiles, cluster_profiles):
    """Assign noise reads to the nearest cluster based on distance to centroids."""
    # Compute centroids for each cluster
    centroids = [compute_centroid(profiles) for profiles in cluster_profiles]

    reassigned_clusters = [[] for _ in centroids]
    for i, noise_profile in enumerate(noise_profiles):
        # Find the closest centroid for each noise read
        distances = [jaccard_distance(noise_profile, centroid) for centroid in centroids]
        closest_cluster = np.argmin(distances)
        reassigned_clusters[closest_cluster].append(noise_reads[i])
    return reassigned_clusters

def cluster_and_save_reads(input_fastq, k=8, eps=0.15):
    """Cluster sequences using DBSCAN, reassign noise reads, and save reads by cluster."""
    # Read sequences from the FASTQ file
    reads = list(SeqIO.parse(input_fastq, "fastq"))
    read_seqs = [str(r.seq) for r in reads]

    # Generate k-mer profiles for each sequence
    profiles = [kmer_profile(seq, k=k) for seq in read_seqs]
    
    # Compute the Jaccard distance matrix
    n = len(read_seqs)
    dist_matrix = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        d = jaccard_distance(profiles[i], profiles[j])
        dist_matrix[i, j] = d
        dist_matrix[j, i] = d

    # Apply DBSCAN clustering
    db = DBSCAN(eps=eps, min_samples=2, metric="precomputed")
    labels = db.fit_predict(dist_matrix)

    # Group reads by cluster and identify noise
    cluster_dict = {}
    noise_reads = []
    noise_profiles = []
    for idx, lbl in enumerate(labels):
        if lbl == -1:  # Noise
            noise_reads.append(reads[idx])
            noise_profiles.append(profiles[idx])
        else:
            cluster_dict.setdefault(lbl, []).append((reads[idx], profiles[idx]))

    # Extract reads and profiles for each cluster
    cluster_reads = [[item[0] for item in cluster_dict[c]] for c in sorted(cluster_dict.keys())]
    cluster_profiles = [[item[1] for item in cluster_dict[c]] for c in sorted(cluster_dict.keys())]

    # Reassign noise reads to the nearest clusters
    reassigned_clusters = assign_noise_to_clusters(noise_reads, noise_profiles, cluster_profiles)

    # Append reassigned noise reads to the clusters
    for i, reassigned in enumerate(reassigned_clusters):
        cluster_reads[i].extend(reassigned)

    # Save clusters to FASTQ files
    base_name = os.path.splitext(os.path.basename(input_fastq))[0]
    for cluster_id, cluster in enumerate(cluster_reads, start=1):
        output_file = f"Cluster_{cluster_id}-{len(cluster)}seqs-{base_name}.fastq"
        SeqIO.write(cluster, output_file, "fastq")
        print(f"Cluster {cluster_id}: {len(cluster)} reads saved to {output_file}")

def main():
    # Look for a FASTQ file in the current working directory
    fastq_files = [f for f in os.listdir() if f.endswith(".fastq")]
    if len(fastq_files) == 0:
        print("No FASTQ files found in the working directory.")
        return

    # Process each FASTQ file in the directory
    for input_fastq in fastq_files:
        print(f"Processing file: {input_fastq}")
        cluster_and_save_reads(input_fastq, k=8, eps=0.15)

if __name__ == "__main__":
    main()





