#pragma rtGlobals=1		// Use modern global access method.

// Blockade_detection_and_measurements_v7
// 06/2019 by Fabien
// detects and measures properties of blockades from a current vs time recording
// what has changed in this version compared to Blockade_detection_and_measurements_v6:
// store values of CrossTh1 for "real" blockades for future display operation

Function BlockDetMeas_v7(SourceWave, Th1, Th2)
	Wave SourceWave
	Variable Th1, Th2
	
	NewDataFolder/O root:BlockDetMeas_v7
	SetDataFolder root:BlockDetMeas_v7
	
	PauseUpdate ; Silent 1
	
	WaveStats/Q SourceWave

	// creates EvCrossTh1 2D wave such as:
	// row number = blockade n°
	// column 0 of row R = n° of last point of SourceWave before crossing Th1 by decrease for blockade n° R
	// column 1 of row R = n° of first point of SourceWave after crossing Th1 by increase for blockade n° R
	Make/O/I/N=(1, 2) EvCrossTh1 // needs the /I flag to create a 32-bit signed integer wave (to avoid problems above 2^24=16 777 216 as single-precision floating point waves have 24 bits of precision)
	EvCrossTh1[0][0]=0
	EvCrossTh1[0][1]=0
	
	Variable i=0, j=0
	
	do
		if((SourceWave[i] >= Th1) && (SourceWave[i+1] < Th1))
			EvCrossTh1[j][0] = i
		elseif ((SourceWave[i] < Th1) && (SourceWave[i+1] >= Th1))
			EvCrossTh1[j][1] = i+1
			j +=1
			Redimension/N=(j+1, -1) EvCrossTh1
		endif
		i +=1
	while (i<(V_npnts-1)) // at the end of this loop, j = total nb of events (possible blockades) from event n° 0 to event n° (j-1)
	Redimension/N=(j, -1) EvCrossTh1
	Variable EvNb = j
	
	// an event begins when the current increases again for the first time after having crossed Th1 by decrease
	// an event ends when the current decreases for the last time before crossing Th1 by increase
	Make/O/I/N=(EvNb, 2) EvEdges // needs the /I flag to create a 32-bit signed integer wave (to avoid problems above 2^24=16 777 216 as single-precision floating point waves have 24 bits of precision)
	// initializes all terms of EvEdges wave to 0
	j=0
	do
		EvEdges[j][0]=0
		EvEdges[j][1]=0
		j +=1
	while (j<EvNb)
	
	Variable test1
	j=0
	do
		test1 = 0
		i = EvCrossTh1[j][0]
		do
			if (SourceWave[i] <= SourceWave[i+1]) // finds the beginning of the event n°j
				EvEdges[j][0] = i
				test1=1
			endif
			i += 1
		while (test1 == 0)
		j += 1
	while (j<EvNb)
	
	j=0
	do
		test1 = 0
		i = EvCrossTh1[j][1]
		do
			if (SourceWave[i-1] >= SourceWave[i]) // finds the end of the event n°j
				EvEdges[j][1] = i
				test1=1
			endif
			i -= 1
		while (test1 == 0)
		j += 1
	while (j<EvNb)
	
	Make/O/D/N=(EvNb) EvIb // needs the /D flag to create a double-precision floating point wave to avoid problems when using values of 32-bits integer waves (EvEdges, BlockEdges) that are greater than 16 777 216
	// initializes all terms of EvIb wave to 0
	j=0
	do
		EvIb[j]=0
		j +=1
	while (j<EvNb)
	
	// calculates Ib for each event (possible blockade)
	j=0
	do
		i = EvEdges[j][0]
		do
			EvIb[j] += SourceWave[i]
			i += 1
		while (i <= EvEdges[j][1]) 
		EvIb[j] = EvIb[j]/(EvEdges[j][1]-EvEdges[j][0]+1) // mean current value during blockade (between EvEdges[j][0] and EvEdges[j][1] both included)
		j +=1
	while (j<EvNb)
	
	// compares Ib of each event (possible blockade) to Th2 to keep only real blockades
	Make/O/I/N=(1, 2) BlockCrossTh1 // needs the /I flag to create a 32-bit signed integer wave (to avoid problems above 2^24=16 777 216 as single-precision floating point waves have 24 bits of precision)
	BlockCrossTh1[0][0]=0
	BlockCrossTh1[0][1]=0
	Make/O/I/N=(1, 2) BlockEdges // needs the /I flag to create a 32-bit signed integer wave (to avoid problems above 2^24=16 777 216 as single-precision floating point waves have 24 bits of precision)
	BlockEdges[0][0]=0
	BlockEdges[0][1]=0
	Make/O/D/N=1 BlockIb // needs the /D flag to create a double-precision floating point wave to avoid problems when using values of 32-bits integer waves (EvEdges, BlockEdges) that are greater than 16 777 216
	BlockIb[0]=0
	
	j=0
	Variable k=0
	do
		if (EvIb[j]<=Th2)
			BlockCrossTh1[k][0]=EvCrossTh1[j][0]
			BlockCrossTh1[k][1]=EvCrossTh1[j][1]
			BlockEdges[k][0]=EvEdges[j][0]
			BlockEdges[k][1]=EvEdges[j][1]
			BlockIb[k]=EvIb[j]
			k +=1
			Redimension/N=(k+1, -1) BlockCrossTh1
			Redimension/N=(k+1, -1) BlockEdges
			Redimension/N=(k+1) BlockIb
		endif
		j +=1
	while (j<EvNb) // at the end of this loop, k = total nb of blockades from blockade n° 0 to blockade n° (k-1)
	Redimension/N=(k, -1) BlockCrossTh1
	Redimension/N=(k, -1) BlockEdges
	Redimension/N=(k) BlockIb
	Variable BlockNb=k
	
	// // needs the /D flag to create double-precision floating point waves to avoid problems when using values of 32-bits integer waves (EvEdges, BlockEdges) that are greater than 16 777 216
	Make/O/D/N=(BlockNb) BlockStdDev // square root of the mean of the square of the distance to the mean
	Make/O/D/N=(BlockNb) BlockBeginTime
	Make/O/D/N=(BlockNb) BlockEndTime
	Make/O/D/N=(BlockNb) BlockDura
	// initializes all terms of BlockStdDev, BlockBeginTime and BlockEndTime waves to 0
	j=0
	do
		BlockStdDev[j]=0
		BlockBeginTime[j]=0
		BlockEndTime[j]=0
		BlockDura[j]=0
		j +=1
	while (j<BlockNb)
	
	// calculates BlockStdDev, BlockBeginTime, BlockEndTime and BlockDura for each blockade
	j=0
	do
		i = BlockEdges[j][0]
		do
			BlockStdDev[j] += (SourceWave[i]-BlockIb[j])^2
			i += 1
		while (i <= BlockEdges[j][1])
		BlockStdDev[j] = BlockStdDev[j]/(BlockEdges[j][1]-BlockEdges[j][0]+1)
		BlockStdDev[j] = sqrt(BlockStdDev[j])
		BlockBeginTime[j] = leftx(SourceWave) + BlockEdges[j][0]*deltax(SourceWave)
		BlockEndTime[j] = leftx(SourceWave) + BlockEdges[j][1]*deltax(SourceWave)
		BlockDura[j] = BlockEndTime[j] - BlockBeginTime[j]		
		j += 1
	while (j<BlockNb)
	
	SetDataFolder root:

End