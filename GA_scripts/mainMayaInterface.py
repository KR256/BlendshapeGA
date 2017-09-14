
maya_dir = '/Users/kylereed/Documents/maya/scripts/'
owd = "/Users/kylereed/Documents/MATLAB/BlendshapeGeneticAlgorithm/"
TEST_CASES = [6,8,11,16]
TEST_ID = 4
NUM_GENERATIONS = 100
TESTER_NAME = 'Kyle'
coarse2Fine = [5,10,20,25,25,25]
automated = False

### 1.
### -- Run this to evolve Identies by sampling from Population --- ###
import sys
sys.path.append(maya_dir)
import samplePopIdentities
reload(samplePopIdentities)
samplePopIdentities.setupBlendshapes(TEST_ID,owd,TEST_CASES,TESTER_NAME)
samplePopIdentities.create18Faces(TEST_ID)
samplePopIdentities.createGeneUI()
if(automated):
    samplePopIdentities.automatedEvolution(NUM_GENERATIONS)


### 2.
### -- Run this to evolve Identity by sampling from PPCA on identity  --- ###
import sys
sys.path.append(maya_dir)
import samplePPCA_identity
reload(samplePPCA_identity)
samplePPCA_identity.setupBlendshapes(TEST_ID,owd,TEST_CASES,TESTER_NAME)
samplePPCA_identity.create18Faces(TEST_ID,coarse2Fine)
samplePPCA_identity.createGeneUI()
if(automated):
    samplePPCA_identity.automatedEvolution(NUM_GENERATIONS)


### 3.
### -- Run this to evolve Smiles by sampling from Population --- ###
import sys
sys.path.append(maya_dir)
import samplePopSmiles
reload(samplePopSmiles)
samplePopSmiles.setupBlendshapes(TEST_ID,owd,TEST_CASES,TESTER_NAME)
samplePopSmiles.create18Faces(TEST_ID)
samplePopSmiles.createGeneUI()
if(automated):
    samplePopSmiles.automatedEvolution(NUM_GENERATIONS)

### 4.
### -- Run this to evolve Smiles by sampling from PPCA on smiles  --- ###
import sys
sys.path.append(maya_dir)
import samplePPCA_smiles
reload(samplePPCA_smiles)
samplePPCA_smiles.setupBlendshapes(TEST_ID,owd,TEST_CASES,TESTER_NAME)
samplePPCA_smiles.create18Faces(TEST_ID,coarse2Fine)
samplePPCA_smiles.createGeneUI()
if(automated):
    samplePPCA_smiles.automatedEvolution(NUM_GENERATIONS)



### - BLENDSHAPE WEIGHT SOLVING - ###


### 5.
### -- Run this to sample a blendshape rig for an identity  --- ###
import sys
sys.path.append(maya_dir)
import sampleFaceBSRig
reload(sampleFaceBSRig)
sampleFaceBSRig.setupBlendshapes(TEST_ID,owd,TEST_CASES,TESTER_NAME)
sampleFaceBSRig.create18Faces(TEST_ID,[10,20,30,40,52,52,52])
sampleFaceBSRig.createGeneUI()
if(automated):
    sampleFaceBSRig.automatedEvolution(NUM_GENERATIONS)


### 6.
### -- Run this to sample solved blendshape smiles in a dataset  --- ###
import sys
sys.path.append(maya_dir)
import sampleSolvedBsSmiles
reload(sampleSolvedBsSmiles)
sampleSolvedBsSmiles.setupBlendshapes(TEST_ID,owd,TEST_CASES,TESTER_NAME)
sampleSolvedBsSmiles.create18Faces(TEST_ID)
sampleSolvedBsSmiles.createGeneUI()
if(automated):
    sampleSolvedBsSmiles.automatedEvolution(NUM_GENERATIONS)

### 7.
### -- Run this to sample from PPCA of solved blendshape smiles in a dataset  --- ###
import sys
sys.path.append(maya_dir)
import samplePPCABSsolves
reload(samplePPCABSsolves)
samplePPCABSsolves.setUpInterface(TEST_ID,owd,TEST_CASES,coarse2Fine,TESTER_NAME)
samplePPCABSsolves.setupBlendshapes(TEST_ID)
samplePPCABSsolves.create18Faces(TEST_ID)
if(automated):
    samplePPCABSsolves.automatedEvolution(NUM_GENERATIONS)