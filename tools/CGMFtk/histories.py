##
#
# @mainpage CGMFtk class to read & analyze CGMF history files
#
# -----------------------------------------------------------------------------
#  CGMF-1.1
#  Copyright TRIAD/LANL/DOE - see file LICENSE
#  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
# -----------------------------------------------------------------------------
#
# @section Description
# The class Histories provides routines to read complete CGMF Monte Carlo history
# files, return arrays or lists of important quantities for further analyses, 
# compute average quantities, distributions and correlations, and provide a 
# summary table of the main average quantities.
#

# Imports
import numpy as np
import sys

# Functions
class Histories:
	
	def __init__ (self, filename, nevents=None):
		"""! Initializes the Histories class

		filename -- file path/name
		nevents -- number of fission events to read in
		"""

		# check that the files exists, otherwise, exit
		try:
			f = open(filename)
			self.filename = filename
			f.close()
		except IOError:
			sys.exit('that file does not exist')

		# check that nevents is > 0, if given
		try:
			val = int(nevents)
			if (val <= 0):
				print ('nevents must be greater than zero')
				sys.exit()
		except ValueError:
			sys.exit('nevents must be a number greater than zero')
		except TypeError:
			if (nevents is not None):
				sys.exit('nevents must be a number')
                                

			
		# read the history file	
		self.histories = self._readHistoryFileFromCGMF (filename,nevents)

		self.numberFragments = len(self.histories)
		self.numberEvents = int(self.numberFragments/2)
		
		# print a warning if nevents is greater than the number of events read
		if (nevents is not None and nevents>self.numberEvents):
			print ('WARNING')
			print ('You asked for ',int(nevents),' events and there are only ',self.numberEvents,' in this history file')
			
		self.A = self.histories[:,0].astype(int)
		self.Z = self.histories[:,1].astype(int)
		self.N = self.A-self.Z
		self.U = self.histories[:,2].astype(np.float64)
		self.J = self.histories[:,3].astype(np.float64)
		self.P = self.histories[:,4].astype(int)
		self.KEpre = self.histories[:,5].astype(np.float64)
		self.KEpost = self.histories[:,24].astype(np.float64)

		numDataOut = len(self.histories[0,:])

		# neutron multiplicities
		self.nu    = self.histories[:,6].astype(int)
		self.nuLF  = self.nu[::2]
		self.nuHF  = self.nu[1::2]
		self.nutot = self.nuLF+self.nuHF
		
		# gamma multiplicities
		self.nug    = self.histories[:,7].astype(int)
		self.nugLF  = self.nug[::2]
		self.nugHF  = self.nug[1::2]
		self.nugtot = self.nugLF+self.nugHF

		# neutron and gamma energies in the center-of-mass of fragments
		self.nEcm   = self.histories[:,8]
		self.nEcmLF = self.nEcm[::2]
		self.nEcmHF = self.nEcm[1::2]

		# neutron and gamma energies in the lab
		self.nElab   = self.histories[:,9]
		self.nElabLF = self.nElab[::2]
		self.nElabHF = self.nElab[1::2]
		
		self.gElab   = self.histories[:,11]
		self.gElabLF = self.gElab[::2]
		self.gElabHF = self.gElab[1::2]

		self.nEcmFragments  = self.nEcm  #-- just fragments (no pre-fission neutrons)
		self.nElabFragments = self.nElab #-- just fragments (no pre-fission neutrons)

		# timing information for gammas - if recorded
		self.gTimes = self.histories[:,12]

		# neutron directional cosines in the lab frame
		self.cmNeutronDircos = self.histories[:,19]
		self.labNeutronDircos = self.histories[:,20]

		#-- pre-fission neutron energies and directional cosines
		self.preFissionNu            = self.histories[1::2,21]
		self.preFissionNeutronElab   = self.histories[1::2,22]
		self.preFissionNeutronDircos = self.histories[1::2,23]

		# Fission fragment momentum vectors (pre-neutron emission)
		pfx = self.histories[:,13]
		pfy = self.histories[:,14]
		pfz = self.histories[:,15]
		self.preFF = np.dstack ((pfx,pfy,pfz))
		self.preFF = self.preFF.reshape((self.numberFragments,3))

		# Fission fragment momentum vectors (post-neutron emission)
		pfx = self.histories[:,16]
		pfy = self.histories[:,17]
		pfz = self.histories[:,18]
		self.postFF = np.dstack ((pfx,pfy,pfz))
		self.postFF = self.postFF.reshape((self.numberFragments,3))

		# Total Kinetic Energy (MeV)
		self.TKEpre = self.KEpre[::2]+self.KEpre[1::2]
		# Total Kinetic Energy (MeV) after neutron emission
		self.TKEpost = self.KEpost[::2]+self.KEpost[1::2]
	
		# Total eXcitation Energy (MeV)
		self.TXE = self.U[::2]+self.U[1::2]

		# Light Fragments
		self.Al = self.A[::2]
		self.Zl = self.Z[::2]
		self.Ul = self.U[::2]
		self.Jl = self.J[::2]
		self.Pl = self.P[::2]
		self.KEl = self.KEpre[::2]

		# Heavy Fragments
		self.Ah = self.A[1::2]
		self.Zh = self.Z[1::2]
		self.Uh = self.U[1::2]
		self.Jh = self.J[1::2]
		self.Ph = self.P[1::2]
		self.KEh = self.KEpre[1::2]


	def _readHistoryFileFromCGMF (self,filename,nevents):
		""" Reads CGMF history file (filename) and returns a list of simulations """

		f = open (filename)

		if (nevents is not None):
			maxNumberEvents = 2*int(nevents)
			hasMaxNumber = True
		else:
			hasMaxNumber = False

		A  = [] #-- fragment masses
		Z  = [] #-- fragment charges
		U  = [] #-- fragment initial excitation energies
		J  = [] #-- fragment initial angular momentum (hbar)
		P  = [] #-- fragment initial parity

		KEpre  = [] #-- pre-neutron emission fragment kinetic energy (MeV)
		KEpost = [] #-- post-neutron emission fragment kinetic energy (MeV)

		nmult    = [] #-- neutron multiplicity (n/f)
		gmult    = [] #-- gamma multiplicity (g/f)
		prenmult = [] #-- pre-fission neutron multiplicity (n/f)

		cmNeutronDircos    = [] #-- neutron directional cosines in center-of-mass frame of the fragment
		cmNeutronEnergies  = [] #-- neutron energies in CM frame (MeV)
		labNeutronDircos   = [] #-- neutron directional cosines in LAB frame
		labNeutronEnergies = [] #-- neutron energies in LAB frame (MeV)

		labPreFissionNeutronEnergies = [] #-- pre-fission neutron energies in LAB frame (MeV)
		labPreFissionNeutronDircos   = [] #-- pre-fission neutron directional cosines in LAB frame

		cmGammaDircos    = [] #-- gamma directional cosines in CM frame
		cmGammaEnergies  = [] #-- gamma energies in CM frame (MeV)
		labGammaDircos   = [] #-- gamma directional cosines in LAB frame
		labGammaEnergies = [] #-- gamma energies in LAB frame (MeV)

		photonAges = [] #-- timing information for gamma rays

		preFragments  = [] #-- pre-neutron emission fragment momentum vectors
		postFragments = [] #-- post-neutron emission fragment momentum vectors

		preFragmentsX = [] #-- pre-neutron emission x-momentum vector
		preFragmentsY = [] #-- pre-neutron emission y-momentum vector
		preFragmentsZ = [] #-- pre-neutron emission z-momentum vector

		postFragmentsX = [] #-- post-neutron emission x-momentum vector
		postFragmentsY = [] #-- post-neutron emission y-momentum vector
		postFragmentsZ = [] #-- post-neutron emission z-momentum vector

		##############################################################################

		line = f.readline()
		if (line.find("#")!=0):
			sys.exit("ERROR: FIRST LINE OF OUTPUT FILE SHOULD CONTAIN '#'")
		else:
			data = line[1:].strip().split()
			ZAIDc   = int(data[0])
			Zc = int(ZAIDc/1000)
			Ac = int(ZAIDc)-1000*Zc + 1 #-- add incident neutron to get CN 
			Asym = int(float(Ac)/2.0)
			Einc    = float(data[1])
			if (Einc==0.0):
				Ac = Ac - 1 #-- no incident neutron for spontaneous fission
			if (len(data)>2):
				self.time = float(data[-1]) # timing cut-off for gammas, -1 records all times of emission
			else: 
				self.time = 0

		isLightFragment=True

		c=0
		index=0
		nfragments=0
		while True:
			c+=1
			# if we only want to read a certain number of events, check here
			if (hasMaxNumber and c>maxNumberEvents):
				break
			line = f.readline()
			if (len(line)==0):
				break
	
			data = line.split()
			nfragments+=1
			
			A.append(int(data[0]))
			Z.append(int(data[1]))
			U.append(float(data[2]))
			J.append(float(data[3]))
			P.append(int(data[4]))
			KEpre.append(float(data[5]))
			KEpost.append(float(data[6]))

			if (isLightFragment):
				Al=int(data[0])
			else:
				Ah=int(data[0])
		
			nn=int(data[7])
			ng=int(data[8])
			npn=int(data[9])

			nmult.append(nn)
			gmult.append(ng)
			prenmult.append(npn) #-- (should always be 0 for LF)

			# read pre- and post-neutron emission fission fragment momentum vectors
			data = f.readline().split()
			preFF=[]
			postFF=[]
			preFF.append  (np.array([float(data[0]),float(data[1]),float(data[2])]))
			postFF.append (np.array([float(data[3]),float(data[4]),float(data[5])]))
			
			pfx = float(data[0])
			pfy = float(data[1])
			pfz = float(data[2])

			pfx2 = float(data[3])
			pfy2 = float(data[4])
			pfz2 = float(data[5])

			#-- read neutron momentum vectors
			# in center of mass frame
			cmDn=[]
			cmEn=[]
			labDn=[]
			labEn=[]
			labDpren=[]
			labEpren=[]

			if (nn>0):

				data=f.readline().split()
				for i in range(nn):
					cmDn.append(np.array([float(data[0+i*4]),float(data[1+i*4]),float(data[2+i*4])]))
					cmEn.append(float(data[3+i*4]))
	
				# in lab frame
				data=f.readline().split()
				for i in range(nn):
					labDn.append(np.array([float(data[0+i*4]),float(data[1+i*4]),float(data[2+i*4])]))
					labEn.append(float(data[3+i*4]))

			#-- read photon momentum vectors (lab=center-of-mass frame in this case, for now)
			Dg=[]
			Eg=[]
			Tg=[]

			if (ng>0 and self.time<0.):
				data=f.readline().split()
				for i in range(ng):
					Dg.append(np.array([float(data[0+i*5]),float(data[1+i*5]),float(data[2+i*5])]))
					Eg.append(float(data[3+i*5]))
					Tg.append(float(data[4+i*5]))
			elif (ng>0):
				data=f.readline().split()
				for i in range(ng):
					Dg.append(np.array([float(data[0+i*4]),float(data[1+i*4]),float(data[2+i*4])]))
					Eg.append(float(data[3+i*4]))
					Tg.append(0.0)

			#-- read pre-fission neutron data (if any)
			if (isLightFragment):
				if (npn!=0):
					print ("npn should be zero for light fragment! -- ABORT!\n")
					exit(-1)
				isLightFragment=False
			else:
				if (Ac-Al-Ah!=npn):
					print ("Incorrect number of pre-fission neutrons! -- ABORT!\n")
				isLightFragment=True

				if (npn>0):
					data = f.readline().split()
					for i in range(npn):
						labDpren.append(np.array([float(data[0+i*4]),float(data[1+i*4]),float(data[2+i*4])]))
						labEpren.append(float(data[3+i*4]))
				else:
					labDpren.append(np.array([0.0,0.0,0.0]))
					labEpren.append(0.0)

			#-- store in final list for output
			
			cmNeutronEnergies.append(cmEn)
			cmNeutronDircos.append(cmDn)
			labNeutronEnergies.append(labEn)
			labNeutronDircos.append(labDn)
			labPreFissionNeutronEnergies.append(labEpren)
			labPreFissionNeutronDircos.append(labDpren)
											 
			cmGammaEnergies.append(Eg)
			cmGammaDircos.append(Dg)
			labGammaEnergies.append(Eg)
			labGammaDircos.append(Dg)
			photonAges.append(Tg)

			preFragmentsX.append (pfx)
			preFragmentsY.append (pfy)
			preFragmentsZ.append (pfz)

			postFragmentsX.append (pfx2)
			postFragmentsY.append (pfy2)
			postFragmentsZ.append (pfz2)

		f.close()
		
		nevents = int(nfragments/2)

		data = np.dstack((A,Z,U,J,P,KEpre,nmult,gmult,cmNeutronEnergies,labNeutronEnergies,cmGammaEnergies, labGammaEnergies,photonAges,preFragmentsX,preFragmentsY,preFragmentsZ,postFragmentsX,postFragmentsY,postFragmentsZ,cmNeutronDircos,labNeutronDircos,prenmult,labPreFissionNeutronEnergies,labPreFissionNeutronDircos,KEpost))


		data = data[0,:,:]

		return (data)

	# Functions to return all quantities recorded
	
	def getFissionHistories(self):
		"""Returns a list with the full simulation history"""
		return (self.histories)

	def getNumberFragments(self):
		"""Returns the total number of fission fragments recorded in a CGMF run
			"""
		return (self.numberFragments)

	def getTimeCoincidenceWindow(self):
		"""Returns the time coincidence window for the gamma rays"""
		return (self.time)
	
	def getNumberEvents(self):
		"""Returns the number of simulated events"""
		return (self.numberEvents)
	
	def getA(self):
		"""Returns a list of the mass (A) for each fragment"""
		return (self.A)

	def getZ(self):
		"""Returns a list of the charge (Z) for each fragment"""
		return (self.Z)

	def getN(self):
		"""Returns a list of the neutron number (N) for each fragment"""
		return (self.N)

	def getU(self):
		"""Returns a list of excitation energy (U) for each fragment"""
		return (self.U)
	
	def getJ(self):
		"""Returns a list of spin (J) for each fragment"""
		return (self.J)

	def getP(self):
		"""Returns a list of parity, +1 or -1, (P) for each fragment"""
		return (self.P)

	def getKEpre(self):
		"""Returns a list of kinetic energy (KE) for each fragment before neutron emission"""
		return (self.KEpre)

	def getKEpost(self):
		"""Returns a list of kinetic energies (KE) for each fragment after neutron emission"""
		return (self.KEpost)

	def getNu(self):
		"""Returns a list of the number of neutrons for each fission fragment (no pre-fission)"""
		return (self.nu)
	
	def getNuEvent(self):
		"""Returns a list of the total number of neutrons for each fission event (no pre-fission)"""
		return (self.nuLF+self.nuHF)

	def getNuEnergyCut(self, Eth):
		"""Returns an array of neutron multiplicity per fragment, with an energy threshold Eth (MeV)"""
		nuCut=[]
		for i in range(self.numberFragments):
			En=np.asarray(self.nElab[i])
			En=En[En>Eth]
			nuCut.append(np.size(En))
		nuCut=np.asarray(nuCut)
		return (nuCut)

	def getNuLF(self):
		"""Returns a list of the number of neutrons for each light fission fragment"""
		return (self.nuLF)

	def getNuHF(self):
		"""Returns the number of neutrons for each heavy fission fragment"""
		return (self.nuHF)

	def getNutot(self):
		"""Returns a list of the total number of neutrons for each fission event including pre-fission neutrons"""
		return (self.nuLF+self.nuHF+self.preFissionNu)
	
	def getPreFissionNu(self):
		"""Returns a list of the number of pre-fission neutrons for each fission event"""
		return (self.preFissionNu)

	def getNug(self):
		"""Returns a list of the number of gammas emitted for each fission fragment"""
		return (self.nug)

	def getNugEnergyCut(self, Eth):
		"""Returns an array of gamma multiplicity per fragment, with an energy threshold Eth (MeV)"""
		nugCut=[]
		for i in range(self.numberFragments):
			Eg=np.asarray(self.gElab[i])
			Eg=Eg[Eg>Eth]
			nugCut.append(np.size(Eg))
		nugCut=np.asarray(nugCut)
		return (nugCut)

	def getNugLF(self):
		"""Returns the number of gammas emitted for each light fission fragment"""
		return (self.nugLF)

	def getNugHF(self):
		"""Returns a list of the number of gammas emitted for each heavy fission fragment"""
		return (self.nugHF)
	
	def getNugtot(self):
		"""Returns a list of the total number of gammas per fission event"""
		return (self.nugLF+self.nugHF)
	
	def getTKEpre(self):
		"""Returns a list of the total kinetic energy per fission event (pre neutron emission)"""
		return (self.TKEpre)

	def getTKEpost(self):
		"""Returns a list of the total kinetic energy per fission event (post neutron emission)"""
		return (self.TKEpost)

	def getTXE(self):
		"""Returns a list of the total excitation energy per fission event"""
		return (self.TXE)

	def getNeutronElab (self):
		"""Returns a list of lists of the neutron energies in the lab frame for each fission fragment"""
		return (self.nElab)

	def getNeutronEcm (self):
		"""Returns a list of lists of the neutron energies in the cm frame for each fission fragment"""
		return (self.nEcm)

	def getGammaElab (self):
		"""Returns a list of lists of the gamma energies for each fission fragment"""
		return (self.gElab)

	def getGammaAges (self):
		"""Returns a list of the gamma times"""

		return (self.gTimes)

	def getPreFissionNeutronElab (self):
		"""Returns a list of lists of the neutrons eneriges before fission for each fission event"""
		return (self.preFissionNeutronElab)

	def getPreFissionNeutronDircos (self):
		"""Returns a list of lists of the neutron directional cosines before fission for each fission event"""
		return (self.preFissionNeutronDircos)
	
	def getFragmentMomentumPre(self):
		"""Returns a list of momenta vectors for each fission fragment before neutron and gamma emission"""
		return (self.preFF)
	
	def getFragmentMomentumPost(self):
		"""Returns a list of momenta vectors for each fission fragment after neutron and gamma emission"""
		return (self.postFF)
	
	def getLabNeutronDircos (self):
		"""Returns a list of neutron directional cosines for each fission fragment in the lab frame"""
		return (self.labNeutronDircos)
	
	def getcmNeutronDircos (self):
		"""Returns a list of neutron directional cosines for each fission fragment in the center of mass frame"""
		return (self.cmNeutronDircos)

	# quantities for the light fragments

	def getALF(self):
		"""Returns a list of the masses of the light fission fragments"""
		return (self.Al)

	def getZLF(self):
		"""Returns a list of the charges of the light fission fragments"""
		return (self.Zl)

	def getNLF(self):
		"""Returns a list of the number of neutrons of the light fission fragments"""
		return (self.Al-self.Zl)

	def getULF(self):
		"""Returns a list of the excitation energy of the light fission fragments"""
		return (self.Ul)

	def getJLF(self):
		"""Returns a list of the spin of the light fission fragments"""
		return (self.Jl)

	def getPLF(self):
		"""Returns a list of the parity of the light fission fragments"""
		return (self.Pl)

	def getKELF(self):
		"""Returns a list of the kinetic energies of the light fission fragments"""
		return (self.KEl)

	# quantities for the heavy fragments

	def getAHF(self):
		"""Returns a list of the masses of the heavy fission fragments"""
		return (self.Ah)

	def getZHF(self):
		"""Returns a list of the charges of the heavy fission fragments"""
		return (self.Zh)

	def getNHF(self):
		"""Returns a list of the number of neutrons of the heavy fission fragments"""
		return (self.Ah - self.Zh)

	def getUHF(self):
		"""Returns a list of the excitation energy of the heavy fission fragments"""
		return (self.Uh)

	def getJHF(self):
		"""Returns a list of the spin of the heavy fission fragments"""
		return (self.Jh)

	def getPHF(self):
		"""Returns a list of the parity of the heavy fission fragments"""
		return (self.Ph)

	def getKEHF(self):
		"""Returns a list of the kinetic energy of the heavy fission fragments"""
		return (self.KEh)

	#################################################################
	#-- average quantities						#
	#################################################################

	def nubar(self):
		"""Returns average neutron multiplicity per fragment, without pre-fission neutrons"""
		return (np.mean(self.nu))

	def nubartot(self):
		"""Returns average neutron multiplicity per fission event, taking pre-fission neutrons into account"""
		return (np.mean(self.nuLF+self.nuHF+self.preFissionNu))

	def nubarg(self,timeWindow=None,Eth=None):
		"""Returns the average gamma multiplicity, per fission fragment

		timeWindow - logical, if included, returns the multiplicity as a function of time
		
		OR
		
		timeWindow - numpy array or list of times at which to calculate the average gamma multiplicity, s

		Eth - lower threshold for gamma energy, MeV
		"""

		if (timeWindow is not None):
			if (timeWindow == True):
				times = np.logspace(-9,0,15)
			else:
				times = timeWindow
			ages = self.getGammaAges()
			gE = self.getGammaElab()

			if (Eth is not None):
				Eth = Eth
			else:
				Eth = 0.
				
			photonAges = []
			photonEnergies = []
			for i in range(len(ages)):
				photonAges += ages[i]
				photonEnergies += gE[i]
			photonAges = np.array(photonAges)
			photonEnergies = np.array(photonEnergies)

			nug = []
			for i in times:	
				mask = np.logical_and(photonAges<=i,photonEnergies>=Eth)
				nug.append(float(len(photonAges[mask])))
			nug = np.array(nug)
			nug = nug/float(self.numberEvents)
			return (times,nug)

		else:
			return (np.mean(self.nug))
	
	def nubargtot(self,timeWindow=None,Eth=None):
		"""Returns the average gamma multiplicity, per fission event

		timeWindow - logical, if included, returns the multiplicity as a function of time
		
		OR
		
		timeWindow - numpy array or list of times at which to calculate the average gamma multiplicity, s

		Eth - lower threshold for gamma energies, MeV
		"""

		if (timeWindow is not None):
			if (timeWindow == True):
				times = np.logspace(-9,0,15)
			else:
				times = timeWindow
			ages = self.getGammaAges()
			gE = self.getGammaElab()
			ages = ages[::2]+ages[1::2]
			gE = gE[::2]+gE[1::2]

			if (Eth is not None):
				Eth = Eth
			else:
				Eth = 0.
				
			photonAges = []
			photonEnergies = []
			for i in range(len(ages)):
				photonAges += ages[i]
				photonEnergies += gE[i]
			photonAges = np.array(photonAges)
			photonEnergies = np.array(photonEnergies)

			nug = []
			for i in times:	
				mask = np.logical_and(photonAges<=i,photonEnergies>=Eth)
				nug.append(float(len(photonAges[mask])))
			nug = np.array(nug)
			nug = nug/float(self.numberEvents)
			return (times,nug)

		else:
			return (np.mean(self.nugLF+self.nugHF))

	def preFissionNubar (self):
		"""Returns the average neutron multiplicity of the pre-fission neutrons"""
		return (np.mean(self.preFissionNu))

	#################################################################
	#-- compute the mean value of a list of lists			#
	#################################################################
	
	def meanList (self, listOfValues):
		"""Returns the mean of a list of lists (such as neutron energy)

		listofValues -- list of lists, e.g. neutron energies, gamma energies
		"""
		l=[]
		for x in listOfValues:
			l+=x
		return (np.mean(np.asarray(l)))

	#-- mean energies of neutrons in the lab frame			
	def meanNeutronElab(self):
		"""Returns the mean neutron energy in the lab frame, with pre-fission neutrons"""
		nE = self.nElab
		nEPE = self.preFissionNeutronElab
		nuPE = self.preFissionNu
		l = []
		for x in nE:
			l += x	
		for i in range(len(nEPE)):
			if (nuPE[i]>0):
				l += nEPE[i]
		
		return (np.mean(l))

	def meanNeutronElabLF (self):
		"""Returns the mean energy of neutrons emitted from the light fragment in the lab frame"""
		return (self.meanList(self.nElabLF))
	
	def meanNeutronElabHF (self):
		"""Returns the mean energy of neutrons emitted from the heavy fragment in the lab frame"""
		return (self.meanList(self.nElabHF))

	def meanNeutronElabFragments (self):
		"Returns the mean neutron energy in the lab frame, without pre-fission neutrons"""
		return (self.meanList(self.nElabFragments))

	#-- mean energies of neutrons in the cm frame

	def meanNeutronEcmFragments (self):
		"""Returns the mean neutron energy in the center of mass frame, without pre-fission neutrons"""
		return (self.meanList(self.nEcmFragments))

	def meanNeutronEcmLF (self):
		"""Returns the mean energy of neutrons emitted from the light fragment in the center of mass frame"""
		return (self.meanList(self.nEcmLF))
	
	def meanNeutronEcmHF (self):
		"""Returns the mean energy of neutrons emitted from the heavy fragment in the center of mass frame"""
		return (self.meanList(self.nEcmHF))

	#-- mean energies of pre-fission neutrons in the lab frame
	def meanPreFissionNeutronElab(self):
		"""Returns the mean energy of pre-fission neutrons in the lab frame"""
		nuPE = self.preFissionNu
		nEPE = self.preFissionNeutronElab
		l = []
		for i in range(len(nEPE)):
			if (nuPE[i]>0):
				l += nEPE[i]
		if (len(l)>0.):
			return (np.mean(l))
		else:
			return(0.0)
	
	#-- mean energies of gammas in the lab frame
	def meanGammaElab(self):
		"""Returns the mean energies of the gammas in the lab frame"""
		return (self.meanList(self.gElab))

	def meanGammaElabLF (self):
		"""Returns the mean energies of the gammas emitted from the light fragments in the lab frame"""
		return (self.meanList(self.gElabLF))
	
	def meanGammaElabHF (self):
		"""Returns the mean energies of the gammas emitted from the heavy fragments in the lab frame"""
		return (self.meanList(self.gElabHF))

	#################################################################
	#-- P(nu)							#
	#################################################################

	def Pnu (self, Eth=None):
		"""Returns a list of probability as a function of neutron multiplicity

		Eth -- optional energy threshold, MeV
		"""

		if (Eth is not None):
			Eth = float(Eth)
			nutot = []
			#nE = self.nElabLF + self.nElabHF + self.preFissionNeutronElab
			nE = self.nElab + self.preFissionNeutronElab
			for x in nE:
				x = np.array(x)
				nutot.append(len(x[x>=Eth]))
			nutot = np.array(nutot)
		else:
			#nutot = self.nuLF + self.nuHF + self.preFissionNu
			nutot = self.nutot +self.preFissionNu

		numax = np.max(nutot)
		nu = np.arange(0,numax+1,1)
		proba = np.zeros(numax+1)
		s=float(self.numberEvents)
		for i in nu:
			proba[i] = len(nutot[nutot==i])/s
		return (nu,proba)

	#################################################################
	#-- P(nug)							#
	#################################################################

	def Pnug (self, Eth=None):
		"""Returns a list of probability as a function of neutron multiplicity

		Eth -- optional energy threshold, MeV
		"""

		if (Eth is not None):
			Eth = float(Eth)
			nugtot = []
			nE = self.gElabLF + self.gElabHF
			for x in nE:
				x = np.array(x)
				nugtot.append(len(x[x>=Eth]))
			nugtot = np.array(nugtot)
		else:
			nugtot = self.nugtot

		nugmax = np.max(nugtot)
		nug = np.arange(0,nugmax+1,1)
		proba = np.zeros(nugmax+1)
		s=float(self.numberEvents)
		for i in nug:
			proba[i] = len(nugtot[nugtot==i])/s
		return (nug,proba)

	#################################################################
	#-- quantities as a function of the fragment mass		#
	#################################################################

	def _qA (self, dummy, quantity):
		"""Returns two arrays, for the fragment mass and the quantity of interest
		
		quantity -- observable of interest as a function of mass, currently supported:
		nuA -- neutron multiplicity 
		nugA -- gamma multiplicity
		UA -- excitation energy
		spinA -- spin
		TKEA -- total kinetic energy
		"""
		xA = np.arange(np.min(self.A),np.max(self.A))
		nA = np.size(xA)
		yA = np.zeros(nA)
		if (quantity=='nuA'):
			q=self.nu
		elif (quantity=='nugA'):
			q=self.nug
		elif (quantity=='UA'):
			q=self.U
		elif (quantity=='TKEA'):
			KE = self.KEpre
			TKE = np.zeros(len(KE))
			TKE[::2] = KE[::2]+KE[1::2]
			TKE[1::2] = TKE[::2]
			q=TKE
		elif (quantity=='spinA'):
			q=self.J
		for i in range(nA):
			r=q[self.A==xA[i]]
			if (r.size):
				yA[i]=np.mean(r)
			else:
				yA[i]=0.0
		return (xA, yA)

	def nubarA (self):
		"""Returns a two-dimensional array for neutron multiplicity as a function of A of format [A,nu]"""
		return (self._qA(self,'nuA'))
	
	def nubargA (self):
		"""Returns a two-dimensional array for gamma multiplicity as a function of A of format [A,nug]"""
		return (self._qA(self, 'nugA'))
	
	def UA (self):
		"""Returns a two-dimensional array for excitation energy as a function of A of format [A,U]"""
		return (self._qA(self, 'UA'))
	
	def TKEA (self):
		"""Returns a two-dimensional array for total kinetic energy as a function of A of format [A,TKE]"""
		return (self._qA(self, 'TKEA'))
	
	def spinA (self):
		"""Returns a two-dimensional array for spin as a function of A of format [A,J]"""
		return (self._qA(self, 'spinA'))

	#################################################################
	#-- quantities as a function of TKE				#
	#################################################################

	def _qTKE (self, dummy, quantity):
		"""Returns two arrays, for the fission event total kinetic energy and the quantity of interest
		
		quantity -- observable of interest as a function of total kinetic energy, currently supported:
		nuTKE -- neutron multiplicity per fission event
		nugTKE -- gamma multiplicity per fission event
		UTKE -- excitation energy
		spinTKE -- spin (J)
		"""
		xTKE = np.arange(np.min(self.TKEpre),np.max(self.TKEpre))
		nTKE = np.size(xTKE)
		yTKE = np.zeros(nTKE)
		if (quantity=='nuTKE'):
			q=self.nutot
			TKE = self.TKEpre
		elif (quantity=='nugTKE'):
			q=self.nugtot
			TKE = self.TKEpre
		elif (quantity=='UTKE'):
			q=self.U
			TKE1 = cgmf.TKEpre
			TKE = []
			for x in TKE1:
				TKE.append(x)
				TKE.append(x)
			TKE = np.array(TKE)
		elif (quantity=='spinTKE'):
			q=self.J
			TKE1 = cgmf.TKEpre
			TKE = []
			for x in TKE1:
				TKE.append(x)
				TKE.append(x)
			TKE = np.array(TKE)
		for i in range(nTKE-1):
			r=q[np.logical_and(TKE>=xTKE[i],TKE<xTKE[i+1])]
			if (r.size):
				yTKE[i]=np.mean(r)
			else:
				yTKE[i]=0.0
		return (xTKE, yTKE)

	def nubarTKE (self):
		"""Returns a two-dimensional array for total neutron multiplicity as a function of TKE of format [TKE,nu]"""
		return (self._qTKE(self, 'nuTKE'))

	def nubargTKE (self):
		"""Returns a two-dimensional array for total gamma multiplicity as a function of TKE of format [TKE,nug]"""
		return (self._qTKE(self, 'nugTKE'))

	def UTKE (self):
		"""Returns a two-dimensional array for excitation energy as a function of TKE of format [TKE,U]"""
		return (self._qTKE(self, 'UTKE'))

	def spinTKE (self):
		"""Returns a two-dimensional array for excitation energy as a function of TKE of foramt [TKE,J]"""
		return (self._qTKE(self, 'spinTKE'))

	#################################################################
	#-- Average Prompt Fission Neutron Spectrum			#
	#################################################################

	def pfns (self,Eth=None):
		"""Returns two arrays, one for the energy grid, and one for the corresponding prompt fission neutron spectrum
		
		Eth - neutron threshold energy, MeV 
		"""
		if (Eth is not None):
			Eth = float(Eth)
		else:
			Eth = 0.0
		# Outgoing energy grid
		egrid = np.logspace(-3,2,100)
		degrid = egrid[1:]-egrid[:-1]
		nElab = self.nElabHF + self.nElabLF + self.preFissionNeutronElab
		l=[]
		for x in nElab:
			l+=x
		l = np.array(l)
		neutronEnergies=l[l>=Eth]
		hlab, binEdges = np.histogram(neutronEnergies,bins=egrid)
		nn=np.sum(hlab)
		binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
		return (binCenters,hlab/degrid/nn)

	#################################################################
	#-- Average Prompt Fission Gamma Spectrum			#
	#################################################################

	def pfgs (self,Eth=None,minTime=None,maxTime=None):
		"""Returns two arrays, one for the energy grid, and one for the corresponding prompt fission gamma spectrum
		
		Eth - optional gamma threshold energy, MeV

		minTime/maxTime - optional lower/upper bound for the time window for the photon emission, s
		"""
		if (Eth is not None):
			Eth = float(Eth)
		else:
			Eth = 0.0
		# Outoing energy grid
		egrid = np.logspace(-1,1,200)
		degrid = egrid[1:]-egrid[:-1]
		l=[]
		gElab = self.gElab
		for x in gElab:
			l+=x
		l = np.array(l)

		if (minTime is not None and maxTime is not None):
			ages = self.getGammaAges()
			a = []
			for x in ages:
				a+=x
			a = np.array(a)
			mask1 = np.logical_and(a>=minTime,a<=maxTime)
			mask2 = np.logical_and(mask1,l>=Eth)
			gammaEnergies = l[mask2]
		else:
			gammaEnergies=l[l>=Eth]
		hlab, binEdges = np.histogram(gammaEnergies,bins=egrid)
		ng=np.sum(hlab)
		binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
		return (binCenters,hlab/degrid/ng*np.mean(self.nugtot))

	#################################################################
	#-- Fission Fragment Angles (relative to beam axis)		#
	#################################################################

	def FFangles (self,afterEmission=None):
		"""Returns cos(theta) of the fission fragments with respect to the z(beam)-axis

		afterEmission - angles after neutron emission (True, default), before neutron emission (False)
		"""

		if (afterEmission is not None):
			afterNeutronEmission = afterEmission
		else:
			afterNeutronEmission = True
		
		if (afterNeutronEmission):
			mom = self.getFragmentMomentumPost()
		else:
			mom = self.getFragmentMomentumPre()

		cosThetaFF = mom[:,2]/np.sum(mom**2,axis=1)**(0.5)
		cosThetaFF = np.array(cosThetaFF.tolist())

		return (cosThetaFF)

	#################################################################
	#-- Neutron Angles (relative to beam axis)			#
	#################################################################

	def nangles (self,Eth=None,lab=None,includePrefission=None):
		"""Returns cos(theta) of the neutrons with respect to the z(beam)-axis
		first for all neutrons, then from light fragments, then from heavy fragments

		Eth - neutron threshold energy, MeV

		lab - True: neutron energy in the lab frame (default), False: in cm

		includePrefission - include (True, default) or don't include (False) pre-fission neutrons
			only for the calculations in the lab frame
		"""
		
		if (Eth is not None):
			Eth = float(Eth)
		else:
			Eth = 0.0
	
		if (lab is not None):
			inLabFrame = lab
		else:
			inLabFrame = True

		if (inLabFrame):
			nE = self.getNeutronElab()
			ncos = self.getLabNeutronDircos()
			nuPre = self.getPreFissionNu()
			nPreE = self.getPreFissionNeutronElab()
			nPrecos = self.getPreFissionNeutronDircos()
			if (includePrefission is not None):
				includePreFissionNeutrons = includePrefission
			else:
				includePreFissionNeutrons = True
		else:
			nE = self.getNeutronEcm()
			ncos = self.getcmNeutronDircos()
			includePreFissionNeutrons = False
			

		ncosAll = []
		ncosLight = []
		ncosHeavy = []
		nEAll = []
		nELight = []
		nEHeavy = []

		Etemp = nE[::2]
		costemp = ncos[::2]
		for i in range(len(costemp)):
			nELight += Etemp[i]
			ncosLight += costemp[i]
			nEAll += Etemp[i]
			ncosAll += costemp[i]

		Etemp = nE[1::2]
		costemp = ncos[1::2]
		for i in range(len(costemp)):
			nEHeavy += Etemp[i]
			ncosHeavy += costemp[i]
			nEAll += Etemp[i]
			ncosAll += costemp[i]

		# have to also include pre-fission neutrons which are seperate now
		if (includePreFissionNeutrons):
			for i in range(len(nuPre)):
				if (nuPre[i]>0):
					nEAll += nPreE[i]
					ncosAll += nPrecos[i]

		nELight = np.array(nELight)
		nEHeavy = np.array(nEHeavy)
		nEAll = np.array(nEAll)
		ncosLight = np.array(ncosLight)
		ncosHeavy = np.array(ncosHeavy)
		ncosAll = np.array(ncosAll)

		ncosAll = ncosAll[nEAll>=Eth][:,2]
		ncosLight = ncosLight[nELight>=Eth][:,2]
		ncosHeavy = ncosHeavy[nEHeavy>=Eth][:,2]

		return (ncosAll,ncosLight,ncosHeavy)

	#################################################################
	#-- Neutron-Fragment Angles 					#
	#################################################################

	def nFangles (self,Eth=None,afterEmission=None,includePrefission=None):
		"""Returns cos(theta) between the neutrons and the fragments

		Eth - neutron threshold energy, MeV
		
		afterEmission - angles after neutron emission (True, default), before neutron emission (False)

		includePrefission - include (True, default) or don't include (False) pre-fission neutrons
			only for the calculations in the lab frame
		"""
		
		if (Eth is not None):
			Eth = float(Eth)
		else:
			Eth = 0.0

		nE = self.getNeutronElab()
		ncos = self.getLabNeutronDircos()
		nuPre = self.getPreFissionNu()
		nPreE = self.getPreFissionNeutronElab()
		nPrecos = self.getPreFissionNeutronDircos()
		if (includePrefission is not None):
			includePreFissionNeutrons = includePrefission
		else:
			includePreFissionNeutrons = True

		if (afterEmission is not None):
			afterNeutronEmission = afterEmission
		else:
			afterNeutronEmission = True
		
		if (afterNeutronEmission):
			mom = self.getFragmentMomentumPost()
		else:
			mom = self.getFragmentMomentumPre()

		mom = mom[::2] # want angles with respect to the light fragments

		cosThetaFF = ((mom.T/np.sum(mom**2,axis=1)**(0.5)).T).tolist()
		cosThetaFF = np.array(cosThetaFF)
			

		ncosAll = []
		ncosLight = []
		ncosHeavy = []
		nEAll = []
		nELight = []
		nEHeavy = []
		FcosLight = []
		FcosHeavy = []
		FcosAll = []

		Etemp = nE[::2]
		costemp = ncos[::2]
		for i in range(len(costemp)):
			nELight += Etemp[i]
			ncosLight += costemp[i]
			nEAll += Etemp[i]
			ncosAll += costemp[i]
			FcosLight += [cosThetaFF[i]]*len(costemp[i])
			FcosAll += [cosThetaFF[i]]*len(costemp[i])

		Etemp = nE[1::2]
		costemp = ncos[1::2]
		for i in range(len(costemp)):
			nEHeavy += Etemp[i]
			ncosHeavy += costemp[i]
			nEAll += Etemp[i]
			ncosAll += costemp[i]
			FcosHeavy += [cosThetaFF[i]]*len(costemp[i])
			FcosAll += [cosThetaFF[i]]*len(costemp[i])

		# have to also include pre-fission neutrons which are seperate now
		if (includePreFissionNeutrons):
			for i in range(len(nuPre)):
				if (nuPre[i]>0):
					nEAll += nPreE[i]
					ncosAll += nPrecos[i]
					FcosAll += [cosThetaFF[i]]*len(nPrecos[i])

		nELight = np.array(nELight)
		nEHeavy = np.array(nEHeavy)
		nEAll = np.array(nEAll)
		ncosLight = np.array(ncosLight)
		ncosHeavy = np.array(ncosHeavy)
		ncosAll = np.array(ncosAll)
		FcosLight = np.array(FcosLight)
		FcosHeavy = np.array(FcosHeavy)
		FcosAll = np.array(FcosAll)

		nFcosAll = ncosAll[nEAll>=Eth]*FcosAll[nEAll>=Eth]
		nFcosAll = np.sum(nFcosAll,axis=1)
		nFcosLight = ncosLight[nELight>=Eth]*FcosLight[nELight>=Eth]
		nFcosLight = np.sum(nFcosLight,axis=1)
		nFcosHeavy = ncosHeavy[nEHeavy>=Eth]*FcosHeavy[nEHeavy>=Eth]
		nFcosHeavy = np.sum(nFcosHeavy,axis=1)

		return (nFcosAll,nFcosLight,nFcosHeavy)	

	#################################################################
	#-- Angular distribution of of n-n opening angle		#
	#################################################################

	def nnangles (self,**keyword_parameters):
		"""Returns an array of the cosine of the angle between each neutron pair of neutrons emitted from the same fragment, first for the light, then heavy, then all fragments

		lab -- neutrons in the lab frame
		
		cm -- neutrons in the center of mass frame
		"""
		nLF = self.nuLF
		nHF = self.nuHF
		ncos = self.labNeutronDircos
		if ('lab' in keyword_parameters):
			ncos = self.labNeutronDircos
		if ('cm' in keyword_parameters):
			ncos = self.cmNeutronDircos

		nncosLF = []
		nncosHF = []
		nncosAll = []

		tempcos = 0.0

		c = 0
		for i in range(len(nLF)):
			for j in range(nLF[i]):
				for k in range(j,nLF[i]):
					tempcos = ncos[c][j][0]*ncos[c][k][0]+ncos[c][j][1]*ncos[c][k][1]+ncos[c][j][2]*ncos[c][k][2]
					nncosLF.append(tempcos)
					nncosAll.append(tempcos)
			c=c+1
			for j in range(nHF[i]):
				for k in range(j,nHF[i]):
					tempcos = ncos[c][j][0]*ncos[c][k][0]+ncos[c][j][1]*ncos[c][k][1]+ncos[c][j][2]*ncos[c][k][2]
					nncosHF.append(tempcos)
					nncosAll.append(tempcos)
			c=c+1

		return (nncosLF,nncosHF,nncosAll)

	#################################################################
	#-- Calculates fragment distribution for a gamma-ray energy	#
	#################################################################
	def gammaSpec(self,Egamma,dEgamma,post=True):

		"""Calculates the distribution of fission fragments (or products), given a gamma-ray energy and energy resolution

		Egamma - array of gamma-ray energies (in lab frame) [MeV]
		
		dEgamma - array of gamma-ray energy resolutions [MeV]
		
		post - True (False) indicates that post- (pre-) neutron emission fission products (fragments) will be returned
		"""

		# Create the list of Z,A pairs
		Zlist = []
		Alist = []

		# Check if post is True or False
		if post:
			multiplier = 1
		else:
			multiplier = 0

		# First validate the gamma-ray energies and their resolutions
		eglist = np.array(Egamma)
		deglist = np.array(dEgamma)

		if np.any(eglist<0.) or np.any(dEgamma<0.):
			print("Gamma-ray energies and/or resolutions should be positive...")
			return

		if not np.shape(eglist) == np.shape(deglist):
			print("Gamma-ray energies and resolutions must be of same shape...")
			return

		if np.shape(eglist)==():
			eglist = np.array([eglist])
			deglist = np.array([deglist])

		# Next, go through the gamma-ray emissions for each event
		for n,gammas in enumerate(self.gElab):

			# Next, check if all of the gamma-rays specified in eglist (within deglist window) are found in gamma rays emitted
			log = np.zeros(len(eglist)).astype(bool)
			for m,eg in enumerate(eglist):

				# Change the log to True for a given gamma-ray window
				if np.any([g <= eg + deglist[m] and g >= eg - deglist[m] for g in gammas]):
					log[m] = True

			# If all gamma-ray gates are satisfied then add this event
			if np.all(log):
				Zlist.append(self.Z[n])
				Alist.append(self.A[n]-multiplier*self.nu[n])

		Zlist = np.array(Zlist)
		Alist = np.array(Alist)
		return np.column_stack((Zlist,Alist))

	#################################################################
	#-- Calculate gamma multiplicity as a function of time 		#
	#################################################################

	def gammaMultiplicity(self,Eth=None,Afragment=None,Zfragment=None,minTime=None,maxTime=None):
		"""Returns the gamma-ray multiplicity as a function of time since the fission event for a specific isotope

		Afragment - mass of a given isotope
		
		Zfragment - charge of a given isotope
		
		If both Afragment and Zfragment are given, the multiplicity is calculated for that isotope (default all isotopes)

		Eth - threshold gamma-ray energy for the calculation, MeV

		minTime/maxTime - define the window over which the number of gamma rays is counted, s
		"""
	
		if (maxTime is not None and minTime is not None):
			tMin = minTime
			tMax = maxTime
			nSteps = 50*np.abs(np.log10(tMax)-np.log10(tMin))+1
		else:
			tMin = 1e-8
			tMax = 100
			nSteps = 501
		timeBins = np.linspace(tMin,tMax,int(nSteps))
		binCenters = 0.5*(timeBins[:-1]+timeBins[1:])

		if (Eth is not None):
			Eth = float(Eth)
		else:
			Eth = 0.0

		Aall = self.getA()
		Zall = self.getZ()
		nu = self.getNu()
		Aall = Aall - nu

		photonAges = []
		photonEnergies = []
		ages = self.getGammaAges()
		gE = self.getGammaElab()

		if (Afragment is not None and Zfragment is not None):
			A = Afragment
			Z = Zfragment
			mask = np.logical_and(Aall==A,Zall==Z)
			for i in range(len(mask)):
				if (mask[i]):
					photonAges += ages[i]
					photonEnergies += gE[i]
				
		else:
			# we calculate the multiplicity from all of the fission fragments
			for i in range(len(ages)):
				photonAges += ages[i]
				photonEnergies += gE[i]

		photonAges = np.array(photonAges)
		photonEnergies = np.array(photonEnergies)
		photonAges = photonAges[photonEnergies>=Eth]

		gMultiplicity,binEdges = np.histogram(photonAges,bins=timeBins)

		return (binCenters,gMultiplicity)

	#################################################################
	#-- Calculate the isomeric ratio of a given state 		#
	#################################################################

	def isomericRatio(self,thresholdTime,A,Z,Jm,Jgs):
		"""Calculates the isomeric ratio for the given state, relative to the ground state
		
		thresholdTime - time to separate isomeric state from ground state, s
		
		A - mass of fragment of interest
		
		Z - charge of fragment of interest
		
		Jm - spin of the isomeric state
		
		Jgs - spin of the ground state
		"""
		
		Aall = self.getA()
		Zall = self.getZ()
		nu = self.getNu()
		Aall = Aall - nu
		ages = self.getGammaAges()
		mask = np.logical_and(Aall==A,Zall==Z)
		
		photonAges = []
		nIsomer = 0
		nGroundState = 0

		for i in range(len(ages)):
			if mask[i]:
				a = np.array(ages[i])
				if (len(a[a>=thresholdTime])>0):
					nIsomer += 1
				else:
					nGroundState += 1

		r = float(nIsomer)/(nIsomer+nGroundState)

		if (Jm>Jgs):
			return (r)
		else:
			return (1.-r)


	#################################################################
	#-- Summary Table for ipython notebook				#
	#################################################################

	def summaryTable (self):
		"""Returns a table for use in an ipython notebook which gives summary information about the CGMF simulation (for the fragments, neutrons, and gammas)"""

		table = ListTable()

		table.append(['','All Fragments','Light Fragments','Heavy Fragments','Pre-Fission', 'Total'])

		#-- A, Z, TXE, TKE, J, U, pi
		table.append(['A',"{0:5.2f}".format(np.mean(self.A)),"{0:5.2f}".format(np.mean(self.Al)),"{0:5.2f}".format(np.mean(self.Ah))])
		table.append(['Z',"{0:5.2f}".format(np.mean(self.Z)),"{0:5.2f}".format(np.mean(self.Zl)),"{0:5.2f}".format(np.mean(self.Zh))])
		table.append(['TXE / U (MeV)',"{0:5.2f}".format(np.mean(self.TXE)),"{0:5.2f}".format(np.mean(self.Ul)),"{0:5.2f}".format(np.mean(self.Uh))])
		table.append(['TKE / KE (MeV)',"{0:5.2f}".format(np.mean(self.TKEpre)),"{0:5.2f}".format(np.mean(self.KEl)),"{0:5.2f}".format(np.mean(self.KEh))])
		table.append(['J ($\hbar$)',"{0:5.2f}".format(np.mean(self.J)),"{0:5.2f}".format(np.mean(self.Jl)),"{0:5.2f}".format(np.mean(self.Jh))])
		table.append(['parity',"{0:5.2f}".format(np.mean(self.P)),"{0:5.2f}".format(np.mean(self.Pl)),"{0:5.2f}".format(np.mean(self.Ph))])

		#-- Neutrons
		table.append([r'$\langle \nu\rangle$',"{0:.3f}".format(np.mean(self.nuLF+self.nuHF)),"{0:.3f}".format(np.mean(self.nuLF)),"{0:.3f}".
			format(np.mean(self.nuHF)),"{0:.3f}".format(np.mean(self.preFissionNu)),"{0:.3f}".format(np.mean(self.nutot+self.preFissionNu))])		
		table.append([r'$\langle \epsilon_n^{cm}\rangle$ (MeV)', "{0:.3f}".format(self.meanNeutronEcmFragments()), 
			"{0:.3f}".format(self.meanNeutronEcmLF()),"{0:.3f}".format(self.meanNeutronEcmHF())])
		table.append([r'$\langle E_n^{lab}\rangle$ (MeV)',"{0:.3f}".format(self.meanNeutronElabFragments()),
			"{0:.3f}".format(self.meanNeutronElabLF()),"{0:.3f}".format(self.meanNeutronElabHF()),
			"{0:.3f}".format(self.meanPreFissionNeutronElab()),"{0:.3f}".format(self.meanNeutronElab())])

		#-- Gammas
		table.append([r'$\langle \nu_\gamma\rangle$',"{0:5.2f}".format(np.mean(self.nugtot)),
			"{0:5.2f}".format(np.mean(self.nugLF)),"{0:5.2f}".format(np.mean(self.nugHF))])
		table.append([r'$\langle E_\gamma^{lab}\rangle$ (MeV)',"{0:5.2f}".format(self.meanGammaElab()),
			"{0:5.2f}".format(self.meanGammaElabLF()),"{0:5.2f}".format(self.meanGammaElabHF())])

		return table


class ListTable(list):
	""" Overridden list class which takes a 2-dimensional list of
		the form [[1,2,3],[4,5,6]], and renders an HTML Table in
		IPython Notebook. """
			
	def _repr_html_(self):
		html = ["<table width=60%>"]
		for row in self:
			html.append("<tr>")
		
			for col in row:
				html.append("<td>{0}</td>".format(col))
			
			html.append("</tr>")
		html.append("</center>")
		return ''.join(html)

