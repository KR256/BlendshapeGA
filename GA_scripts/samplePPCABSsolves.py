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

def setUpInterface(sID,envName,tIDs,nPCs,oN):
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
	global owd
	global numPCs
	global outName
	global SmileChosenLog
	global generation
	global BlendshapeNames
	global BlendshapeNamesFull
	global Num_BlendShapes
	global Num_BlendShapesFull
	global owd
	global testIDs

	testIDs = tIDs
	startingID = testIDs[sID-1]
	owd = envName
	numPCs = nPCs
	SAMPLE_RANGE = 2.5
	pcNum = 25
	outName = oN
	SmileLog = []
	SmileChosenLog = []
	numVerts = 4906
	generation = 0
	os.chdir(owd)
	filePath = "Blendshapes/smile/weightMatrixPPCA.mat"
	matData = spio.loadmat(filePath)

	avgSmileVec = matData['mu']
	avgSmileVec = avgSmileVec[0]
	eigenvalues = matData['pcvar']
	components = matData['coeff']
	W = matData['W']
	v = matData['v']
	smileScores = [None] * 18

	with open('TedTargets/blendshape_Fnames.txt') as f:
		BlendshapeNames = f.read().splitlines()

	with open('TedTargets/blendshape_names.txt') as f:
		BlendshapeNamesFull = f.read().splitlines()

	Num_BlendShapes = len(BlendshapeNames)
	Num_BlendShapesFull = len(BlendshapeNamesFull)


def setupBlendshapes(FACE_ID,*args):
	global owd
	global startingID
	os.chdir(owd)
	neutral_Filename = "Blendshapes/neutral/bsNeutrals/FACE"+str(startingID)+"_Neutral.obj"
	# neutral_Filename = "neutrals/bathFemale000.obj"
	# neutral_Filenames = glob.glob("processedNeutrals/*.obj")
	# neutral_Filename = neutral_Filenames[FACE_ID-1]
	print neutral_Filename
	cmds.file(neutral_Filename,i=True,groupName="blendshapes")
	cmds.select('Mesh')
	mel.eval('polySetToFaceNormal')


	os.chdir("Blendshapes/smile/deformedTestTargets/FACE_"+str(startingID)+"/")
	obj_Filenames = glob.glob("*.obj")
	num_Bs = len(obj_Filenames)

	for i in range(num_Bs):
		filePath = obj_Filenames[i]
		filePathWtExt = filePath.rsplit('.', 1)[0]
		print 'Loading Mesh file: ' + filePath
		cmds.file(filePath,i=True,groupName="blendshapes")
		cmds.select(filePathWtExt + '_Mesh')
		mel.eval('polySetToFaceNormal')
		cmds.move(500 + 200*((i)%10),-300.0*((i)/10),0,r=True)

	cmds.group('FACE*',n='blendshapes',w=True)
	# cmds.select('faces',all=True)

	## Pressing Blendshape button
	mel.eval('select -r "FACE*"')
	mel.eval('select -tgl "Mesh"')
	mel.eval('blendShape;')

def create18Faces(startID,*args):
	global startingID
	global randomTargetList
	global SmileLog
	global smileScores

	randomTargetList = np.random.permutation(64)
	newIDs = [x - 1 for x in testIDs]
	randomTargetList = np.setdiff1d(randomTargetList,newIDs)
	randomTargetList = np.random.permutation(randomTargetList)

	for face in range(1,19):
		cmds.select('Mesh')
		faceN = 'face'+str(face)
		cmds.duplicate('Mesh',n=faceN)
		mel.eval('select -r "FACE*"')
		mel.eval('select -tgl '+faceN)
		mel.eval('blendShape;')
		# print 'Face ' + str(face) + ' using Target smile ' + str(randomTarget)
		
		addRandomSmile(face,generation)

		cmds.select('face'+str(face))
		cmds.move(200*((face-1)%6),1200 + (-300.0*((face-1)/6)),0,r=True)
	cmds.group('face*',n='faces',w=True)
	cmds.select('faces')
	cmds.viewFit()

	smileScoresCopy = smileScores[:]
	SmileLog.append(smileScoresCopy)
	# print SmileLog

	createGeneUI()
	# createGeneUI()

def addRandomSmile(faceID,generation,*args):
	global smileScores
	global startingID

	personalisationsScore = getPersonalisationScore(generation)

	print "Random Smile:" + str(faceID)
	print personalisationsScore
	smileScores[faceID-1] = personalisationsScore
	smileFaceOut = scoreToMesh(faceID,personalisationsScore,generation)
	setBlendShapePCs(faceID,smileFaceOut,startingID)

def getPersonalisationScore(pcLim,*args):

	global SAMPLE_RANGE
	global numPCs

	componentsTemp = components.transpose()
	componentsTemp = componentsTemp[0:numPCs[pcLim]]
	componentsTemp = componentsTemp.transpose()

	randomSample = np.random.uniform(-SAMPLE_RANGE, SAMPLE_RANGE, 25)
	WTW = np.dot(W.transpose(),W)
	WlinWTW = np.linalg.solve(WTW,W.transpose())
	brig = np.multiply(np.identity(25),v)
	WTWbrig = WTW + brig
	WTWrand = np.dot(WTWbrig,randomSample.transpose())
	WTWrand = np.array(WTWrand)
	Y_hat = np.dot(WlinWTW.transpose(),WTWrand)
	score = np.dot(Y_hat.transpose(),componentsTemp)

	return score

def scoreToMesh(faceID,score,pcLim,*args):

	global avgSmileVec
	global components
	global identityMesh

	componentsTemp = components.transpose()
	componentsTemp = componentsTemp[0:numPCs[pcLim]]
	componentsTemp = componentsTemp.transpose()

	addVec = np.dot(componentsTemp,score.transpose())
	meshVec = addVec + avgSmileVec 

	return meshVec

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
	global startingID
	global randomTargetList
	global FaceBSLog
	global generation
	global numPCs

	generation += 1
	print 'Generation: ' + str(generation)
	if generation > len(numPCs)-1:
		numPCs.append(25)

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
	print chosenCombos
	print len(chosenCombos)

	chosenSmiles = []
	for j in range(len(chosenFaces2)):
		chosenSmiles.append(oldSmileScores[int(chosenFaces2[j])-1])
	eliteFaceScore = oldSmileScores[int(bestFace)-1]
	chosenSmiles.insert(0,eliteFaceScore)
	SmileChosenLog.append(chosenSmiles[:])

	for i in range(1,19):
		if i==1:
			eliteFaceScore = list(oldSmileScores[int(bestFace)-1])
			eliteFaceScore = padWithZeros(eliteFaceScore,generation)
			smileScores[i-1] = eliteFaceScore
			eliteMesh = scoreToMesh(i,eliteFaceScore,generation)
			setBlendShapePCs(1,eliteMesh,startingID)
			print 'Face 1 chosen as ELITE face: '+ bestFace
			print eliteFaceScore
		elif i==2:
			eliteFaceScore = np.array(oldSmileScores[int(bestFace)-1])
			chosenSum = eliteFaceScore
			for t in range(len(chosenFaces2)):
				chosenSum += np.array(oldSmileScores[chosenFaces2[t]-1])
			avgFace = chosenSum / (len(chosenFaces2)+1)
			avgFaceScore = padWithZeros(avgFace,generation)
			smileScores[i-1] = avgFaceScore
			avgFaceMesh = scoreToMesh(i,avgFaceScore,generation)
			setBlendShapePCs(i,avgFaceMesh,startingID)
			print 'Face 2 chosen as average of: '+ bestFace + ' and ' + str(chosenFaces2)
			print avgFaceScore
		elif i==3:
			eliteFaceScore = np.array(oldSmileScores[int(bestFace)-1])
			for t in range(len(eliteFaceScore)):
				if( generation == 1):
					eliteFaceScore[t] = eliteFaceScore[t] * 1.5
				elif(t >= numPCs[generation-2] and generation > 1):
					eliteFaceScore[t] = eliteFaceScore[t] * 1.5
			exagFaceScore = padWithZeros(eliteFaceScore,generation)
			smileScores[i-1] = exagFaceScore
			exagFaceMesh = scoreToMesh(i,exagFaceScore,generation)
			setBlendShapePCs(i,exagFaceMesh,startingID)
			print 'Face 3 chosen as capped sum of: '+ bestFace + ' and ' + str(chosenFaces2)
			print exagFaceScore
		elif i<(len(chosenFaces2)+4):
			if i < len(chosenFaces2)+4:
				f1 = bestFace
				f2 = chosenFaces2[i-4]
				print 'Face '+ str(i) +' chosen as average of ELITE and Face ' + str(f2)
			# else:
			# 	j = i - (len(chosenFaces2)+4)
			# 	chosenTemp = chosenCombos[j]
			# 	f1 = chosenTemp[0]
			# 	f2 = chosenTemp[1]
			# 	print 'Face '+ str(i) +' chosen as average of Face ' + str(f1) +' and Face ' + str(f2)
			f1PCs = np.array(oldSmileScores[int(f1)-1])
			f2PCs = np.array(oldSmileScores[int(f2)-1])
			meanPCs = (f1PCs + f2PCs) / 2
			meanPCs = padWithZeros(meanPCs,generation)
			print meanPCs
			smileScores[i-1] = meanPCs
			aMesh = scoreToMesh(i,meanPCs,generation)
			setBlendShapePCs(i,aMesh,startingID)
		elif i<12:
			f1 = random.choice(chosenFaces)
			numChanges = 3
			if(numPCs[generation-1] < 25):
				chosenPCs = random.sample(range(numPCs[generation-1],numPCs[generation]),numChanges)
			else:
				chosenPCs = random.sample(range(25),numChanges)
			cFace = list(oldSmileScores[f1-1])
			cFace = padWithZeros(cFace,generation)
			# print 'cFace:'
			# print cFace
			ranFace = mutateBlendshapes(f1,0.7,chosenPCs,f1,cFace)
			smileScores[i-1] = ranFace
			muMesh = scoreToMesh(i,ranFace,generation)
			setBlendShapePCs(i,muMesh,startingID)
			print 'Face '+ str(i) + ' chosen as mutation of Face ' + str(f1) +' for PC ' + str(chosenPCs) 
			print ranFace
		else:
			# randomFace = randomTargetList[0]
			# randomTargetList = randomTargetList[1:]
			# ranFace = mutateBlendshapes(randomFace)
			# setBlendShapePCs(i,ranFace,startingID)
			addRandomSmile(i,generation)
			print 'Face ' + str(i) + ' spawned at random'

	smileScoresCopy = smileScores[:]
	SmileLog.append(smileScoresCopy)

	# global groundTruth
	# global bsVertsMat

	# smileBSscores = createFaceBSLog()

	# faceVec = blendshapesToVec(smileBSscores[0],bsVertsMat)
	# print faceVec
	# res = np.array(groundTruth) - np.array(faceVec)
	# print res
	# sse = np.sqrt(np.sum(res**2))
	# print sse
	# sseList.append(sse)
	# print sseList

def getBlendShapePCs(faceID,startingID,*args):
	cFaceList = []
	for pc in range(1,Num_BlendShapes+1):
		bsName = 'blendShape' +str(faceID+1)+ '.' + 'FACE' + str(startingID) + '_' + BlendshapeNames[pc-1] + '_Mesh'
		cFaceList.append(cmds.getAttr( bsName))
	return cFaceList
    
def setBlendShapePCs(faceID,bshapeList,startingID,*args):
	bshapeList = np.clip(bshapeList,0,1)
	for pc in range(1,Num_BlendShapes+1):
		bsName = 'blendShape' +str(faceID+1)+ '.' + 'FACE' + str(startingID) + '_' + BlendshapeNames[pc-1] + '_Mesh'
		# print bsName
		# print bshapeList[pc-1]
		cmds.setAttr( bsName,bshapeList[pc-1])

def createFaceBSLog(*args):
	global startingID
	FaceBSLog = []
	for i in range(1,19):
		cFaceList = getBlendShapePCs(i,startingID)
		FaceBSLog.append(cFaceList)
	return FaceBSLog

def mutateBlendshapes(face,mutateRate,chosenPCs,f1,oldScores,*args):

	for i in range(len(chosenPCs)):
		ranPCVal = oldScores[chosenPCs[i]]
		ranNewVal = random.uniform(-4*mutateRate*np.sqrt(eigenvalues[chosenPCs[i]-1]), 4*mutateRate*np.sqrt(eigenvalues[chosenPCs[i]]))
		oldScores[chosenPCs[i]] = ranNewVal
	return oldScores

def padWithZeros(scores,gen,*args):

	global numPCs

	if(numPCs[gen-1] < 25):
		zerosToAdd = [0] * (numPCs[gen] - numPCs[gen-1])
	else:
		zerosToAdd = []
	returnMesh = np.concatenate((scores, zerosToAdd), axis=0)

	return returnMesh

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
	fileName2 = 'SmileLog_FACE' + str(startingID) + '_ppcaBlendshapeSmiles_' + str(outName) + '.mat'
	fileName3 = 'SmileChosenLog_FACE' + str(startingID) + '_ppcaBlendshapeSmiles_' + str(outName) + '.mat'
	spio.savemat(fileName2, {"SmileLog": SmileLog})
	spio.savemat(fileName3, {"SmileChosenLog": SmileChosenLog})
	print 'Saved Smile Log:' + fileName2
	print 'Saved Chosen Log:' + fileName3
	print SmileLog
	print SmileChosenLog

	eliteFace = getBlendShapePCs(int(bestFace),startingID)
	fileName = 'FACE' + str(startingID) + '_' + 'samplePPCABlendshapes_' + outName + '.mat'
	spio.savemat(fileName, {"BlendshapeWeights":eliteFace})


def automatedEvolution(NUM_GENERATIONS,*args):

	global FaceBSLog
	global nextBest
	global nextElite
	global smileScores
	global groundTruth
	global smileBSdeltas
	global sseList
	global bsVertsMat
	global startingID
	sseList = []


	SSE_gens = []

	os.chdir(owd)

	filePath = "Blendshapes/smile/deformedTestTargets/FACE_"+str(startingID)+"/bsVerts52_"+str(startingID)+".mat"
	matData = spio.loadmat(filePath)
	bsVertsMat = matData['bsVertsMat']


	filePath = "Blendshapes/neutral/bsNeutrals.mat"
	matData = spio.loadmat(filePath)
	testBSneutrals = matData['bsNeutrals']
	bsNeutral = testBSneutrals[startingID-1]

	# filePath = "Resources/testBlendshapeSolveDeltas.mat"
	# matData = spio.loadmat(filePath)
	# testBlendshapeSolveDeltas = matData['testBlendshapeSolveDeltas']
	# groundTruth = testBlendshapeSolveDeltas[FACE_ID-1]

	filePath = "Resources/smileDeltasTraining.mat"
	matData = spio.loadmat(filePath)
	smileDeltas = matData['exprDeltas']
	groundTruth = smileDeltas[startingID-1]

	bsVertsMat = np.array(bsVertsMat) - np.array(bsNeutral)

	# print 'Ground Truth'
	# print groundTruth

	for gen in range(NUM_GENERATIONS):

		smileBSscores = createFaceBSLog()

		genResults = compareWithGroundTruth(groundTruth,bsVertsMat,smileBSscores)
		SSE_gens.append(genResults['smallestSSE'])

		nextElite = genResults['eliteIndex']
		nextBest = genResults['bestFaces']

		generateNextGen(nextElite,nextBest,False)


	print SSE_gens

	print 'Final SSE: ' + str(SSE_gens[-1])



def compareWithGroundTruth(groundTruth,smileBSdeltas,FaceMatrix,*args):

	global startingID
	

	sseList = []
	for face in range(1,19):
		faceBlendshapes = FaceMatrix[face-1]
		print 'Face Blendshape:' + str(face)
		print faceBlendshapes
		faceVec = blendshapesToVec(faceBlendshapes,smileBSdeltas)
		print 'Face vec:' + str(face)
		print faceVec
		res = np.array(groundTruth) - np.array(faceVec)
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














