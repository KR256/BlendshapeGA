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



def setupBlendshapes(FACE_ID,envName,tIDs,oN,*args):
	global outName
	global gen
	global BlendshapeWeightMatrix
	global BlendshapeNames
	global BlendshapeNamesFull
	global Num_BlendShapes
	global Num_BlendShapesFull
	global startingID
	global randomTargetList
	global owd
	global testIDs

	testIDs = tIDs
	startingID = testIDs[FACE_ID-1]
	owd = envName
	os.chdir(owd)

	with open('TedTargets/blendshape_Fnames.txt') as f:
	    BlendshapeNames = f.read().splitlines()

	with open('TedTargets/blendshape_names.txt') as f:
	    BlendshapeNamesFull = f.read().splitlines()


	Num_BlendShapes = len(BlendshapeNames)
	Num_BlendShapesFull = len(BlendshapeNamesFull)
	gen = 0
	outName = oN
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

def create18Faces(startID,pcN,*args):
	global startingID
	global pcNum
	pcNum = pcN

	for face in range(1,19):
		cmds.select('Mesh')
		faceN = 'face'+str(face)
		cmds.duplicate('Mesh',n=faceN)
		mel.eval('select -r "FACE*"')
		mel.eval('select -tgl '+faceN)
		mel.eval('blendShape;')

		numChanges = 3
		chosenPCs = random.sample(range(pcNum[0]),numChanges)
		ranFace = randBlendshapes(0.8,chosenPCs)
		setBlendShapePCs(face,ranFace,startingID)
		print 'Face ' + str(face) + ' spawned at random' + '... BSs chosen:' + str(chosenPCs)
		for j in range(numChanges):
			print 'Changed: ' + BlendshapeNames[chosenPCs[j]]
		
		cmds.select('face'+str(face))
		cmds.move(200*((face-1)%6),1200 + (-300.0*((face-1)/6)),0,r=True)
	cmds.group('face*',n='faces',w=True)
	cmds.select('faces')
	cmds.viewFit()

	# createGeneUI()

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
	global FaceBSLog
	global pcNum
	global gen

	gen += 1
	print 'Generation: ' + str(gen)
	if gen > len(pcNum)-1:
		pcNum.append(52)
	FaceBSLog = createFaceBSLog()

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

	for i in range(1,19):
		if i==1:
			eliteFace = list(FaceBSLog[int(bestFace)-1])
			setBlendShapePCs(1,eliteFace,startingID)
			print 'Face 1 chosen as ELITE face: '+ bestFace
		elif i==2:
			eliteFace = np.array(FaceBSLog[int(bestFace)-1])
			chosenSum = eliteFace
			for t in range(len(chosenFaces2)):
				chosenSum += np.array(FaceBSLog[chosenFaces2[t]-1])
			avgFace = chosenSum / (len(chosenFaces2)+1)
			setBlendShapePCs(i,avgFace,startingID)
			print 'Face 2 chosen as average of: '+ bestFace + ' and ' + str(chosenFaces2)
			print avgFace
		elif i==3:
			eliteFace = np.array(FaceBSLog[int(bestFace)-1])
			chosenSum = eliteFace
			for t in range(len(chosenFaces2)):
				for r in range(len(FaceBSLog[0])):
					tempBS = np.array(FaceBSLog[chosenFaces2[t]-1])
					if(chosenSum[r] == 0):
						chosenSum[r] += tempBS[r]
			mergedFace = np.clip(chosenSum,0,1)
			setBlendShapePCs(i,mergedFace,startingID)
			print 'Face 3 chosen as capped sum of: '+ bestFace + ' and ' + str(chosenFaces2)
			print mergedFace
		elif i<(len(chosenFaces2)):
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
			f1PCs = np.array(FaceBSLog[int(f1)-1])
			f2PCs = np.array(FaceBSLog[int(f2)-1])
			meanPCs = (f1PCs + f2PCs) / 2
			setBlendShapePCs(i,meanPCs,startingID)
			print meanPCs
		elif i<12:
			f1 = random.choice(chosenFaces)
			numChanges = 3
			if(pcNum[gen-1] == 52):
				chosenPCs = random.sample(range(52),numChanges)
			else:
				chosenPCs = random.sample(range(pcNum[gen-1],pcNum[gen]),numChanges)
			ranFace = mutateBlendshapes(f1,0.5,chosenPCs)
			setBlendShapePCs(i,ranFace,startingID)
			print 'Face '+ str(i) + ' chosen as mutation of Face ' + str(f1)  + '... BSs chosen:' + str(chosenPCs)
		else:
			# randomFace = randomTargetList[0]
			# randomTargetList = randomTargetList[1:]
			# ranFace = mutateBlendshapes(randomFace)
			# setBlendShapePCs(i,ranFace,startingID)
			numChanges = 5
			chosenPCs = random.sample(range(pcNum[gen]),numChanges)
			ranFace = randBlendshapes(0.7,chosenPCs)
			setBlendShapePCs(i,ranFace,startingID)
			print 'Face ' + str(i) + ' spawned at random' + '... BSs chosen:' + str(chosenPCs)

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
	for pc in range(1,Num_BlendShapes+1):
		bsName = 'blendShape' +str(faceID+1)+ '.' + 'FACE' + str(startingID) + '_' + BlendshapeNames[pc-1] + '_Mesh'
		cmds.setAttr( bsName,bshapeList[pc-1])

def createFaceBSLog(*args):
	global startingID
	FaceBSLog = []
	for i in range(1,19):
		cFaceList = getBlendShapePCs(i,startingID)
		FaceBSLog.append(cFaceList)
	return FaceBSLog

def mutateBlendshapes(face,mutateRate,chosenPCs,*args):
	global Num_BlendShapes
	global FaceBSLog
	global startingID
	mutateFace = getBlendShapePCs(face,startingID)

	for i in range(len(chosenPCs)):
		mutateFace[chosenPCs[i]] = random.uniform(0,mutateRate)
	print mutateFace
	return mutateFace

def randBlendshapes(mutateFace,chosenPCs,*args):
	global Num_BlendShapes
	randWeights = []

	randBlends = np.zeros(Num_BlendShapes)
	
	#Smile blendshape

	randBlends[2] = random.uniform(0.5, 1)
	for i in range(len(chosenPCs)):
		randBlends[chosenPCs[i]] = random.uniform(0, mutateFace)

	return randBlends

def saveFace(bF,*args):
	global startingID
	global outName
	global owd

	os.chdir(owd + '/resultingLogs')
	FaceBSLog = createFaceBSLog()
	print FaceBSLog
	bestFace = `cmds.intField( bF, query = True, value = True)`
	eliteFace = FaceBSLog[int(bestFace)-1]
	fileName = 'FACE' + str(startingID) + '_' + 'sampleRig_' + outName + '.mat'
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
	global owd
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

	# BS TARGET instead
	# filePath = "Resources/testBlendshapeSolveDeltas.mat"
	# matData = spio.loadmat(filePath)
	# testBlendshapeSolveDeltas = matData['testBlendshapeSolveDeltas']
	# groundTruth = testBlendshapeSolveDeltas[FACE_ID-1]

	filePath = "Resources/smileDeltasFull.mat"
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












