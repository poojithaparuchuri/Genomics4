"""
ORF HMM Script
"""
#Import all the necessary libraries)
import numpy as np
import matplotlib.pyplot as plt
import random


# Define a Hidden Markov Model (HMM) class for generating Open Reading Frame (ORF) sequences
class HMM:
  # Initialize HMM with states, start probabilities, transition probabilities, and emission probabilities
  def __init__(self):
    self.states = ['Start_Codon', 'Coding_Sequence', 'Stop_Codon']
    self.start_probability = np.array([1.0, 0.0, 0.0])
    self.transition_probability = np.array([
        [0.05, 0.94, 0.01],
        [0.0, 0.95, 0.05],
        [0.0, 0.0, 1.0]
    ])

    self.emission_probs = {
        'Start_Codon': {'ATG': 1.0},
        'Coding_Sequence': {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25},
        'Stop_Codon': {'TAA': 0.34, 'TAG': 0.33, 'TGA': 0.33}
    }

  # Generate a single ORF sequence using the HMM with a minimum length condition
  def sample(self):
    while True:
      orf = ''
      state = np.random.choice(self.states, p=self.start_probability)

      # Continue generating the sequence until a stop codon is encountered
      while True:
        if state == 'Start_Codon':
          orf += 'ATG'
        elif state == 'Coding_Sequence':
          # Emit a nucleotide based on the emission probabilities
          nucleotide = random.choices(list(self.emission_probs[state]), weights=list(self.emission_probs[state].values()))[0]
          orf += nucleotide
        else:
          # Emit a stop codon and break out of the loop
          stop_codon = random.choices(list(self.emission_probs[state]), weights=list(self.emission_probs[state].values()))[0]
          orf += stop_codon
          break

        # Update the state based on transition probabilities
        transition_probability = self.transition_probability[self.states.index(state)]
        state = np.random.choice(self.states, p=transition_probability)

      # Check if the generated ORF length is greater than 10
      if len(orf) > 10:
        return orf


# Main execution
if __name__ == '__main__':
  # Create an instance of the HMM class
  hmm = HMM()

  # Generate a specified number of ORF sequences
  num_samples = 250
  sequences = []
  while len(sequences) < num_samples:
    orf = hmm.sample()
    sequences.append(orf)

  lengths = [len(seq) for seq in sequences]

  # Plot a histogram of ORF lengths
  plt.hist(lengths, bins=range(min(lengths), max(lengths) + 4, 3), edgecolor='black')
  plt.title('Length Distribution of ORFs (Length > 10)')
  plt.xlabel('Length of the ORFs')
  plt.ylabel('Frequency')
  plt.show()

  # Write ORF lengths and sequences to a file (ensure proper indentation)
  with open('SreePoojitha_ORFoutput.txt', 'w') as f:  # Open in write mode
    for l, s in zip(lengths, sequences):
      f.write(f'{l}, {s}\n')
