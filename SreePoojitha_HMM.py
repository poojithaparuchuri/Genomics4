"""
Sree Poojitha Paruchuri
Problem Set-5
HMM model to train and compare DNA sequences.
"""

import math
import matplotlib.pyplot as plt

# Define a dictionary to map nucleotides to their corresponding indices
baseIDx = {"A": 0, "C": 1, "G": 2, "T": 3}


def main():
    # Define file paths for input sequences
    spaciiFA = "MSpacii.fa"
    pathogenFA = "pathogen.fa"
    spaciiFA_T = "MSpacii_training.fa"
    pathogenFA_T = "pathogen_training.fa"
    # Retrieve sequences from input files
    spaciiID2seq = getSeq(spaciiFA)
    pathogenID2seq = getSeq(pathogenFA)

    # Initialize and train Hidden Markov Models for Spacii and Pathogen sequences
    spaciiTrainModel = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    pathTrainModel = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    spaciiTrainModel = trainModel(spaciiTrainModel, spaciiFA_T)
    pathTrainModel = trainModel(pathTrainModel, pathogenFA_T)

    # Calculate and store log-likelihood scores for Spacii and Pathogen sequences
    markovScoresSpacii = []
    markovScoresPath = []

    for ID in spaciiID2seq.keys():
        markovScoresSpacii.append(getLogLike(spaciiTrainModel, pathTrainModel, spaciiID2seq[ID]))
    for ID in pathogenID2seq.keys():
        markovScoresPath.append(getLogLike(spaciiTrainModel, pathTrainModel, pathogenID2seq[ID]))

    # Plot a histogram of log-likelihood scores for Spacii and Pathogen sequences
    plt.hist([markovScoresPath, markovScoresSpacii], bins=20, label=['pathogen', 'spacii'], rwidth=1, density=True)
    plt.show()

    # Write log-likelihood scores to a file
    scoresOutputText(markovScoresSpacii, markovScoresPath)


def scoresOutputText(markovScoresSpacii, markovScoresPath):
    # Write log-likelihood scores to a tab-delimited file
    f = open("nucleotide_composition_HMM_results.tab", "w")
    f.write("SpaciiScores\tpathogenScores\n")
    for i in range(len(markovScoresSpacii)):
        f.write(str(markovScoresSpacii[i])+"\t"+str(markovScoresPath[i])+"\n")
    f.close()


def getLogLike(model1, model2, seq):
    # Calculate the log-likelihood score for a given sequence based on two models
    Pmod1 = 1
    Pmod2 = 1

    for i in range(len(seq) - 1):
        x = baseIDx[seq[i]]
        y = baseIDx[seq[i + 1]]
        Pmod1 *= model1[x][y]
        Pmod2 *= model2[x][y]

    score = math.log(Pmod1) - math.log(Pmod2)
    return score


def trainModel(model, data):
    # Train a Hidden Markov Model by counting transitions in the training data
    id2seq = getSeq(data)

    for seq in id2seq.values():
        for i in range(len(seq) - 1):
            x = baseIDx[seq[i]]
            y = baseIDx[seq[i + 1]]
            model[x][y] += 1

    # Normalize to turn counts into probabilities
    for row in model:
        s = sum(row)
        for i in range(len(row)):
            row[i] = row[i] / s
    print(model)
    return model


def getSeq(filename):
    # Read sequences from a FASTA file and store them in a dictionary
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.find(">") == 0:
            currkey = line.rstrip()[1:]
            id2seq[currkey] = ""
        else:
            id2seq[currkey] = id2seq[currkey] + line.rstrip()
    return id2seq


# Execute the main function
main()
