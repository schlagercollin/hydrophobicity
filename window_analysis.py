CTD = "KHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAA\
AQAALQSSWGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGWGSASNAGSGS\
GFNGGFGSSMDSKSSGWGM"

fensterSize = 5 #Define the window size for analsysis

hydroDict = {
    "F": 1.79,
    "I": 1.8,
    "L": 1.7,
    "M": 1.23,
    "V": 1.22,
    "W": 2.25,
    "Y": 0.96,
    "A": 0.31,
    "P": 0.72,
    "N": -0.6,
    "Q": -0.22,
    "T": 0.26,
    "C": 1.54,
    "G": 0.00,
    "S": -0.04,
    "D": -0.77,
    "E": -0.64,
    "H": 0.13,
    "K": -0.99,
    "R": -1.01
}
# Values obtained from original publication Fauchere and Pliska, \ 1983

def analyze(sequence):
    Ntotal = len(sequence)
    Npos = sequence.count("R") + sequence.count("K")
    Nneg = sequence.count("D") + sequence.count("E")
    
    fPos = Npos / Ntotal
    fNeg = Nneg / Ntotal

    FCR = fPos + fNeg
    NCPR = fPos - fNeg
    # ChargeAsymmetry = (fPos + fNeg)**2/(fPos - fNeg) # div by zero error?

    hydrophobicity = 0
    for aminoAcid in sequence:
        hydrophobicity += hydroDict[aminoAcid]
    
    relativeHydrophobicity = hydrophobicity/Ntotal

    return {
            "sequence" : sequence, 
            "relHydro" : relativeHydrophobicity, 
            "NCPR" : NCPR,
            "fPos" : fPos,
            "fNeg" : fNeg
            }

def window_analysis(sequence, windowSize):
    outputData = []
    for index, value in enumerate(sequence):
        window = sequence[index:index+windowSize]
        if len(window) != windowSize:
            break
            # Terminate if at end. Does this not mess up analysis of end??
        
        outputData.append(analyze(window))
    return outputData

def main(protein_seq, window_size=5):
    window_size = int(window_size)
    output = window_analysis(protein_seq, window_size)
    x_data = []
    y_data = []
    for index, datapoint in enumerate(output):
        x_data.append(index)
        y_data.append(datapoint["relHydro"])
    return x_data, y_data, protein_seq



