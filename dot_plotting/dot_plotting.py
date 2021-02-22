from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt

def dot_plot(seq_record,comparison_sequence,complement=False,window=3):
    subject_strand = seq_record
    seq_two = comparison_sequence
    data = np.array([[int((subject_strand[i:i + window] != seq_two[j:j + window]))
                       for i in range(len(subject_strand) - window)]
                      for j in range(len(seq_two) - window)])
    return data

def make_quick_plot(window, sequence_1, sequence_2):
    print( 'window size: ' + str(window))
    plt.rcParams['figure.figsize'] = 10,10
    plt.imshow(dot_plot(sequence_1,sequence_2,complement=False,window=window),cmap="Greys_r",interpolation='none')
    plt.show()

def count_matches(seq_A,seq_B):
    count = 0
    for i in range(0,len(seq_A)):
        if seq_A[i] == seq_B[i]:
            count += 1
    return count

def dot_plot_tolerant(seq_record,comparison_sequence,window=3):
    subject_strand = seq_record
    seq_two = comparison_sequence
    data = np.array([[count_matches(subject_strand[i:i + window], seq_two[j:j + window])
                       for i in range(len(subject_strand) - window)]
                      for j in range(len(seq_two) - window)])
    print( data)
    return data
def make_quick_plot_tolerant(window, sequence_1,sequence_2):
    print( 'window size: ' + str(window))
    plt.rcParams['figure.figsize'] = 10,10
    plt.imshow(dot_plot_tolerant(sequence_1,sequence_2,window=window),cmap="viridis",interpolation='none')
    plt.show()

if __name__ == "__main__":
    first_sequence = Seq('GCTAGCTAGTAGCTTAGGATGATCGTACGTAGCTAGCTGATTATAGAGAGAGAAGGAGAA')
    second_sequence = Seq('TTCGCTTGCTCTCTCTATAATCAGTTAGCTTCGTACGATCATCCTAAGGTACTAGCTAGC').reverse_complement()
    make_quick_plot(window=6, sequence_1=first_sequence, sequence_2=second_sequence)
    plt.show()
    make_quick_plot_tolerant(window=10, sequence_1=first_sequence, sequence_2=second_sequence)
