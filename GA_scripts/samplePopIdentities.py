import maya
import maya.cmds as cmds
import os
import sys
import numpy as np
import glob
import maya.mel as mel
import scipy.io as spio
from functools import partial
import itertools
import random

global BlendshapeWeightMatrix
global startingID
global randomTargetList




def setupBlendshapes(FACE_ID,envName,tIDs,oN,*args):

	global numVerts
	global identityMesh
	global smileScores
	global avgSmileVec
	global matData
	global eigenvalues
	global components
	global W
	global v
	global SAMPLE_RANGE
	global SmileLog
	global startingID
	global testIDs
	global owd
	global numPCs
	global outName
	global SmileChosenLog
	global generation
	global num_Bs

	testIDs = tIDs
	owd = envName
	os.chdir(owd)
	startingID = testIDs[FACE_ID-1]
	SAMPLE_RANGE = 0.05
	pcNum = 25
	outName = oN
	SmileLog = []
	SmileChosenLog = []
	numVerts = 4906
	generation = 0

	smileScores = [None] * 18
	gen = 0
	outName = oN
	# neutral_Filename = "neutrals/FACE"+str(FACE_ID)+"_Neutral.obj"
	# neutral_Filename = "neutrals/bathFemale000.obj"
	filePath = "Resources/averageIdentity.obj"
	print filePath
	fileO = cmds.file(filePath,i=True, rnn=True, ignoreVersion = True, options = "mo=0",  loadReferenceDepth  = "all")
	cmds.rename(fileO[4], 'Mesh')
	cmds.select('Mesh1')
	mel.eval('polySetToFaceNormal')

	os.chdir("Dataset/OBJs_neutrals")
	obj_Filenames_full = glob.glob("*.obj")
	obj_Filenames = obj_Filenames_full[:]
	for tr in range(len(testIDs)):
		obj_Filenames.remove(obj_Filenames_full[testIDs[tr]-1])
	num_Bs = len(obj_Filenames)

	for i in range(num_Bs):
		bsN = 'BLENDSHAPE'+str(i+1)
		cmds.scriptEditorInfo(suppressWarnings=True)
		filePath = obj_Filenames[i]
		filePathWtExt = filePath.rsplit('.', 1)[0]
		print 'Loading Mesh file: ' + filePath
		faceObj = cmds.file(filePath,i=True, rnn=True, ignoreVersion = True, options = "mo=0",  loadReferenceDepth  = "all")
		print faceObj[0]
		cmds.rename(faceObj[0], bsN)
		cmds.select(bsN)
		mel.eval('polySetToFaceNormal')
		cmds.move(500 + 200*((i)%10),-300.0*((i)/10),0,r=True)

	cmds.group('BLENDSHAPE*',n='blendshapes',w=True)
	# cmds.select('faces',all=True)

	## Pressing Blendshape button
	mel.eval('select -r "BLENDSHAPE*"')
	mel.eval('select -tgl "Mesh1"')
	mel.eval('blendShape;')



def useIDsBlendShapeWeights(inputID,bsID,startingID,*args):

	bSLog = []
	for bs in range(1,Num_BlendShapes+1):
		bsName = 'blendShape' +str(inputID)+ '.' + 'FACE' + str(startingID) + '_' + BlendshapeNames[bs-1] + '_Mesh'
		bsWeight = BlendshapeWeightMatrix[bsID-1,bs-1]
		# print bsName
		# print bsWeight
		#mel.eval('setAttr "blendShape1.FACE1_DUAL_smile_CTL_Mesh" 0.405634')
		scriptName = "setAttr" + " \"" + bsName + "\" " + str(bsWeight) + ";"
		# scriptName2 = "setAttr" + " \"" + "Mesh." + bsName + "\" " + str(bsWeight) + ";"
		# print scriptName
		# mel.eval(scriptName)
		# mel.eval('setAttr "blendShape1.FACE1_DUAL_smile_CTL_Mesh" 0.405634')
		cmds.setAttr(bsName,bsWeight)

def create18Faces(startID,*args):
	global startingID
	global numPCs
	global generation
	global randomTargetList
	global smileScores

	randomTargetList = setUpRandomTarget()
	print randomTargetList

	for face in range(1,19):
		cmds.select('Mesh1')
		faceN = 'face'+str(face)
		faceN2 = 'face'+str(face+1)
		cmds.duplicate('Mesh1',n=faceN)
		print 'Face should be' + str(faceN)
		mel.eval('select -r "BLENDSHAPE*"')
		mel.eval('select -tgl '+faceN)
		mel.eval('blendShape;')

		smileOut = addRandomSmile(face)
		smileScores[face-1] = smileOut
		print 'Face ' + str(face) + ' spawned at random'

		
		cmds.select('face'+str(face))
		cmds.move(200*((face-1)%6),1200 + (-300.0*((face-1)/6)),0,r=True)
	cmds.group('face*',n='faces',w=True)
	cmds.select('faces')
	cmds.viewFit()

	# createGeneUI()

def addRandomSmile(faceID,*args):
	global randomTargetList
	global identityMesh
	global smileScores
	global startingID

	try:
		randomSmile = randomTargetList[0]
	except IndexError:
		print 'WARNING'
		randomTargetList = setUpRandomTarget()
		randomSmile = randomTargetList[0]

	print 'Random smile chosen as:' + str(randomTargetList[0]+1)
	smileOut = np.zeros(num_Bs)
	smileOut[randomSmile] = 1.0
	# smileScores[faceID-1] = smileOut
	setBlendShapePCs(faceID,smileOut,startingID)
	randomTargetList = np.delete(randomTargetList,0)

	return smileOut


def createGeneUI(*args):
	global startingID
	print startingID
	myWindow = "Smile Selector"
	if maya.cmds.window(myWindow,ex=True):
		maya.cmds.deleteUI(myWindow)
	cmds.window(myWindow, width = 300, height = 200) 
	cmds.columnLayout("columnLayoutName01", adjustableColumn=True )
	cmds.text("3DFITtitle", label = "Check boxes of most similar faces:", width = 20, height = 20, backgroundColor = [0.2, 0.2, 0.2])
	cmds.gridLayout( numberOfColumns=6, nr = 3, cellWidthHeight=(50, 50),parent = "columnLayoutName01" )
	cBox = []
	for cB in range(18):
		cBox.append(cmds.checkBox( label = str(cB+1) ,align='right'))
		#cBvals[cB-1] = cmds.checkBox(cBox, query = True, value = True)
	cmds.separator(parent = "columnLayoutName01")
	cmds. rowLayout("nameRowLayout01", numberOfColumns = 2, parent = "columnLayoutName01")
	cmds.text( label = "Select ID of MOST similar face:", width = 200, height = 20, backgroundColor = [0.2, 0.2, 0.2],parent = "nameRowLayout01")
	bestFace = cmds.intField( "bestFace", minValue = 1, maxValue = 18, value = 18,editable=True, parent = "nameRowLayout01")
	cmds.separator(parent = "columnLayoutName01")
	cmds.button(label = "Next Generation", command = partial(generateNextGen,bestFace,cBox,True),parent = "columnLayoutName01")
	cmds.button(label = "Save Face", command = partial(saveFace,bestFace),parent = "columnLayoutName01")
	cmds.showWindow()

def generateNextGen(bF,cBvals,interFaceCond,*args):

	global identityMesh
	global smileScores
	global SmileLog
	global SmileChosenLog
	global num_Bs
	global generation


	generation += 1

	oldSmileScores = smileScores[:]

	chosenFaces = []

	if(interFaceCond==True):
		for cB in range(18):
			if cmds.checkBox( cBvals[cB],q=True, v=True):
				chosenFaces.append(cB+1)

		bestFace = `cmds.intField( bF, query = True, value = True)`
	else:
		bestFace = str(bF)
		chosenFaces = cBvals

	chosenFaces2 = chosenFaces[:]

	if int(bestFace) in chosenFaces:
		chosenFaces2.remove(int(bestFace))
	print chosenFaces2

	chosenCombos = itertools.combinations(chosenFaces2, 2)
	chosenCombos = list(chosenCombos)

	chosenSmiles = []
	for j in range(len(chosenFaces2)):
		chosenSmiles.append(oldSmileScores[int(chosenFaces2[j])-1])
	eliteFaceScore = oldSmileScores[int(bestFace)-1]
	chosenSmiles.insert(0,eliteFaceScore)
	SmileChosenLog.append(chosenSmiles[:])

	temp = [None] * 18
	for i in range(1,19):
		if i==1:
			eliteFaceScore = list(oldSmileScores[int(bestFace)-1])
			temp[0] = eliteFaceScore
			setBlendShapePCs(i,eliteFaceScore,startingID)
			print 'Face 1 chosen as ELITE face: '+ bestFace
			print eliteFaceScore
		elif i==2:
			eliteFaceScore = np.array(oldSmileScores[int(bestFace)-1])
			chosenSum = eliteFaceScore
			for t in range(len(chosenFaces2)):
				chosenSum += np.array(oldSmileScores[chosenFaces2[t]-1])

			avgFace = chosenSum / (len(chosenFaces2)+1)
			temp[1] = avgFace
			setBlendShapePCs(i,avgFace,startingID)
			print 'Face 2 chosen as average of: '+ bestFace + ' and ' + str(chosenFaces2)
			print avgFace

		elif i<(len(chosenFaces2)+len(chosenCombos)+3):
			if i < len(chosenFaces2)+3:
				f1 = bestFace
				f2 = chosenFaces2[i-3]
				print 'Face '+ str(i) +' chosen as average of ELITE and Face ' + str(f2)
			else:
				j = i - (len(chosenFaces2)+3)
				chosenTemp = chosenCombos[j]
				f1 = chosenTemp[0]
				f2 = chosenTemp[1]
				print 'Face '+ str(i) +' chosen as average of Face ' + str(f1) +' and Face ' + str(f2)
			f1PCs = np.array(oldSmileScores[int(f1)-1])
			f2PCs = np.array(oldSmileScores[int(f2)-1])
			meanPCs = (f1PCs + f2PCs) / 2

			print meanPCs
			temp[i-1] = meanPCs
			setBlendShapePCs(i,meanPCs,startingID)
		elif i <12:
			print 'i equals' + str(i)
			f1 = random.choice(chosenFaces)

			numChanges = 3

			chosenPCs = random.sample(range(num_Bs),numChanges)
			oS = list(oldSmileScores)
			cFace = list(oS[f1-1])
			print 'cFace:'
			print cFace
			ranFace = mutateBlendshapes(cFace,0.5,chosenPCs)
			# print ranFace
			temp[i-1] = ranFace
			setBlendShapePCs(i,ranFace,startingID)
			print 'Face '+ str(i) + ' chosen as mutation of Face ' + str(f1) +' for PC ' + str(chosenPCs) 
			print ranFace
		else:
			smileOut = addRandomSmile(i)
			temp[i-1] = smileOut
			print 'Face ' + str(i) + ' spawned at random'

		# print temp

	smileScores = temp[:]
	SmileLog.append(smileScores)

	# global groundTruth
	# global neutralBSdeltas
	# faceVec = blendshapesToVec(list(oldSmileScores[int(bestFace)-1]),neutralBSdeltas)
	# print faceVec
	# res = np.array(groundTruth) - np.array(faceVec)
	# print res
	# sse = np.sqrt(np.sum(res**2))
	# print sse
	# sseList.append(sse)
	# print sseList
	return temp[:]


def getBlendShapePCs(faceID,startingID,*args):
	cFaceList = []
	for pc in range(num_Bs):
		bsName = 'blendShape' +str(faceID+1)+ '.' + 'BLENDSHAPE' + str(pc+1)
		cFaceList.append(cmds.getAttr( bsName))
	return cFaceList
    
def setBlendShapePCs(faceID,bshapeList,startingID,*args):
	for pc in range(len(bshapeList)):
		bsName = 'blendShape' +str(faceID+1)+ '.' + 'BLENDSHAPE' + str(pc+1)
		cmds.setAttr( bsName,bshapeList[pc])

def createFaceBSLog(*args):
	global startingID
	FaceBSLog = []
	for i in range(1,19):
		cFaceList = getBlendShapePCs(i,startingID)
		FaceBSLog.append(cFaceList)
	return FaceBSLog

def mutateBlendshapes(mutateFace,mutateRate,chosenPCs,*args):

	for i in range(len(chosenPCs)):
		mutateFace[chosenPCs[i]] = random.uniform(0,mutateRate)
	return mutateFace


def randBlendshapes(mutateFace,chosenPCs,*args):
	global Num_BlendShapes
	randWeights = []

	randBlends = np.zeros(Num_BlendShapes)
	
	#Smile blendshape

	for i in range(len(chosenPCs)):
		randBlends[chosenPCs[i]] = random.uniform(0, mutateFace)

	return randBlends

def saveFace(bF,*args):
	global SmileLog
	global startingID
	global owd
	global smileScores
	global SmileChosenLog
	global outName

	os.chdir(owd + '/resultingLogs')
	bestFace = `cmds.intField( bF, query = True, value = True)`
	oldSmileScores = smileScores[:]
	eliteFaceScore = oldSmileScores[int(bestFace)-1]
	SmileChosenLog.append(eliteFaceScore)
	# eliteFace = getFaceVec(int(bestFace),'glob')
	# fileName = 'Composite_FACE' + str(startingID) + '_ppcaSmiles_' + bestFace + '.mat'
	# spio.savemat(fileName, {"eliteFace": eliteFace})
	# print 'Saved End Face: ' + fileName
	fileName2 = 'SmileLog_FACE' + str(startingID) + '_samplePopSmileBSv_' + str(outName) + '.mat'
	fileName3 = 'SmileChosenLog_FACE' + str(startingID) + '_samplePopSmileBSv_' + str(outName) + '.mat'
	spio.savemat(fileName2, {"SmileLog": SmileLog})
	spio.savemat(fileName3, {"SmileChosenLog": SmileChosenLog})
	print 'Saved Smile Log:' + fileName2
	print 'Saved Chosen Log:' + fileName3

	faceID = int(bestFace)
	cmds.select('face' + bestFace)
	cmds.move(-200*((faceID-1)%6),300.0*((faceID-1)/6),0,r=True)

	eliteFace = getFaceVec(int(bestFace),'glob')
	fileName = 'FACE' + str(startingID) + '_samplePopSmileBSv_' + outName + '.mat'
	spio.savemat(fileName, {"smileResult":eliteFace})


def setUpRandomTarget(*args):

	global num_Bs
	global startingID
	global testIDs

	rTList = np.random.permutation(num_Bs)
	SELECTED_TEST_FACES = [x - 1 for x in testIDs]
	rTList = np.setdiff1d(rTList,SELECTED_TEST_FACES)
	rTList = np.random.permutation(rTList)

	return rTList


def automatedEvolution(NUM_GENERATIONS,*args):

	global FaceBSLog
	global nextBest
	global nextElite
	global smileScores
	global groundTruth
	global neutralBSdeltas
	global sseList
	global owd
	global testIDs
	global startingID

	sseList = []

	SSE_gens = []

	os.chdir(owd)

	filePath = "Blendshapes/neutral/neutralBSdeltas.mat"
	matData = spio.loadmat(filePath)
	neutralBSdeltas = matData['neutralBSDeltas']
	groundTruth = neutralBSdeltas[startingID-1]
	newIDs = [x - 1 for x in testIDs]
	neutralBSdeltas = np.delete(neutralBSdeltas,newIDs,0)



	for gen in range(NUM_GENERATIONS):

		genResults = compareWithGroundTruth(groundTruth,neutralBSdeltas,smileScores)
		SSE_gens.append(genResults['smallestSSE'])

		nextElite = genResults['eliteIndex']
		nextBest = genResults['bestFaces']

		smileScores = generateNextGen(nextElite,nextBest,False)

		print 'SMILES UPDATED'


	print SSE_gens

	print 'Final SSE: ' + str(SSE_gens[-1])



def compareWithGroundTruth(gTruth,neutralBSdeltas,SmileData,*args):

	global startingID
	# global smileScores


	# filePath = "Resources/GA_neutrals.mat"
	# matData = spio.loadmat(filePath)
	# neutrals = matData['neutrals']
	# groundTruth = neutrals[startingID-1]

	# print SmileData

	sseList = []
	for face in range(1,19):
		faceBlendshapes = SmileData[face-1]
		print 'Face Blendshape:' + str(face)
		print faceBlendshapes
		faceVec = blendshapesToVec(faceBlendshapes,neutralBSdeltas)
		res = np.array(gTruth) - np.array(faceVec)
		sse = np.sqrt(np.sum(res**2))
		sseList.append(sse)

	print sseList
	sortedSSE = np.argsort(sseList)
	smallestSSE = sseList[sortedSSE[0]]
	sortedSSE = [x+1 for x in sortedSSE]
	eliteIndex = sortedSSE[0]
	bestFaces = sortedSSE[:4]

	print 'Elite Face:' + str(eliteIndex)
	print 'Best Faces:' + str(bestFaces[:4])
	print 'Smallest SSE' + str(smallestSSE)

	return {'eliteIndex':eliteIndex, 'bestFaces':bestFaces ,'smallestSSE':smallestSSE }

def blendshapesToVec(blendshapes,neutralBSdeltas,*args):
	faceVec = np.matmul(blendshapes, neutralBSdeltas)

	return faceVec












