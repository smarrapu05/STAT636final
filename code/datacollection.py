import numpy as np
from Bio.PDB import PDBList, PDBParser
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
from Bio.Align import substitution_matrices
from scipy.spatial.distance import pdist, squareform
import os

PDBIDS = {'Mb': '1MBN', 'Hb': '4HHB'}
Chains = {'Mb': 'A', 'Hb': 'A'}
DistThreshold = 7

def getCalphaCoords(pdb_id,chain_id):
    #get pdb
    pdbl = PDBList()
    pdb_filename = pdbl.retrieve_pdb_file(pdb_id,pdir='.',file_format='pdb')
    #parse pdb
    parser = PDBParser()
    structure = parser.get_structure(pdb_id,pdb_filename)
    chain = structure[0][chain_id]
    coords = []
    sequence = ''
    res_indices = []
    #create sequence to find CA
    for residue in chain:
        if residue.id[0] != ' ':
            continue

        if 'CA' in residue:
            coords.append(residue['CA'].get_coord())
            sequence += seq1(residue.get_resname())
            res_indices.append(residue.id[1])

    return np.array(coords), sequence, np.array(res_indices)

print("getting calpha")
coordsMb, seqMb, idxMb = getCalphaCoords(PDBIDS['Mb'], Chains['Mb'])
coordsHb, seqHb, idxHb = getCalphaCoords(PDBIDS['Hb'], Chains['Hb'])

print("sequence align")
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

alignment = aligner.align(seqMb,seqHb)[0]

alignedSeqMb = alignment[0]
alignedSeqHb = alignment[1]

alignedIndMb = []
alignedIndHb = []
idxCountMb = 0
idxCountHb = 0

for i in range(len(alignedSeqMb)):
    charMb = alignedSeqMb[i]
    charHb = alignedSeqHb[i]

    notGapMb = charMb != '-'
    notGapHb = charHb != '-'

    if notGapMb and notGapHb:
        alignedIndMb.append(idxCountMb)
        alignedIndHb.append(idxCountHb)
    if notGapMb:
        idxCountMb += 1
    if notGapHb:
        idxCountHb += 1

print(f"aligned core size = {len(alignedIndMb)}")

finalCoordsMb = coordsMb[alignedIndMb]
finalCoordsHb = coordsHb[alignedIndHb]

def getAdjacency(coords, threshold):
    dists = pdist(coords,metric='euclidean')
    distMat = squareform(dists)
    A = (distMat < threshold).astype(int)
    np.fill_diagonal(A,0)
    return A

Amb = getAdjacency(finalCoordsMb,DistThreshold)
Ahb = getAdjacency(finalCoordsHb,DistThreshold)

print("adjacency mats")
print(f"Amb Shape = {Amb.shape}")
print(f"Ahb Shape = {Ahb.shape}")

np.savetxt("adjMb.csv",Amb,delimiter=",",fmt="%d")
np.savetxt("adjHb.csv",Ahb,delimiter=",",fmt="%d")

import csv
print("node mapping")

mappingData = []
for i in range(len(alignedIndMb)):
    mbIdx = alignedIndMb[i]
    hbIdx = alignedIndHb[i]
    mbRes = seqMb[mbIdx]
    hbRes = seqHb[hbIdx]
    mbPDBnum = idxMb[mbIdx]
    hbPDBnum = idxHb[hbIdx]

    mappingData.append([i+1,mbRes,mbPDBnum,hbRes,hbPDBnum])

with open('nodeMapping.csv','w',newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Node_Index', 'Mb_Residue', 'Mb_PDB_ID', 'Hb_Residue', 'Hb_PDB_ID'])
    writer.writerows(mappingData)

print(f"saved mapping size {len(mappingData)} !")
