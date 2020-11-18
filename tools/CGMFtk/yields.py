##
#
# @mainpage CGMFtk class to read & analyze CGMF initial fission fragment yields
#
# -----------------------------------------------------------------------------
#  CGMF-1.0
#  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
#  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
# -----------------------------------------------------------------------------
#
# @section Description
# The class Yields provides routines to read CGMF yields, as obtained using the 
# 'n' option with a negative number (corresponding to the number of initial
# conditions considered). 
#

# Imports
import numpy as np
import sys

class Yields:

	#################################################################
	#-- define the Yields						#
	#################################################################

	def __init__ (self, filename, nevents=None):
		"""Initializes the Yields class

		filename -- file path/name
		nevents -- never of events to read in
		"""

		# check that the file exists, otherwise exit
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
				sys.exit('nevents must be greater than zero')
		except ValueError:
			sys.exit('nevents must be a number greater than zero')
		except TypeError:
			if (nevents is not None):
				sys.exit('nevents must be a number')
                
		self.yields = self.readYieldFileFromCGMF(filename,nevents)

		self.numberFragments = len(self.yields)
		self.numberEvents = int(self.numberFragments/2)
		
		# print a warning if nevents is greater than the number of events read
		if (nevents is not None and nevents>self.numberEvents):
			print ('WARNING')
			print ('You asked for ',int(nevents),' events and there are only ',self.numberEvents,' in this yield file')
		
		self.A = self.yields[:,0].astype(int)
		self.Z = self.yields[:,1].astype(int)
		self.N = self.A - self.Z
		self.U = self.yields[:,2].astype(np.float64)
		self.J = self.yields[:,3].astype(np.float64)
		self.P = self.yields[:,4].astype(int)
		self.KE = self.yields[:,5].astype(np.float64)
		
		# fragment momenta
		pfx = self.yields[:,6].astype(np.float64)
		pfy = self.yields[:,7].astype(np.float64)
		pfz = self.yields[:,8].astype(np.float64)
		self.pFF = np.dstack((pfx,pfy,pfz))
		self.pFF = self.pFF.reshape((self.numberFragments,3))

		# Total Kinetic Energy
		self.TKE = self.KE[::2] + self.KE[1::2]

		# Total eXcitation Energy
		self.TXE = self.U[::2] + self.U[1::2]

		# Light Fragments
		self.Al = self.A[::2]
		self.Zl = self.Z[::2]
		self.Ul = self.U[::2]
		self.Jl = self.J[::2]
		self.Pl = self.P[::2]
		self.KEl = self.KE[::2]


		# Heavy Fragments
		self.Ah = self.A[1::2]
		self.Zh = self.Z[1::2]
		self.Uh = self.U[1::2]
		self.Jh = self.J[1::2]
		self.Ph = self.P[1::2]
		self.KEh = self.KE[1::2]

	

	#################################################################
	# read yield file from CGMF					#
	#################################################################

	def readYieldFileFromCGMF (self, filename, nevents=None):
		"""Reads the CGMF history file (filename) and returns a list of the yield characteristics"""

		f = open (filename)

		if (nevents is not None):
			maxNumberEvents = 2*int(nevents)
			hasMaxNumber = True
		else:
			hasMaxNumber = False

		A = [] #-- fragment masses
		Z = [] #-- fragment charges
		KE = [] #-- fragment kinetic energy
		U = [] #-- fragment excitation energy
		J = [] #-- fragment spin
		P = [] #-- fragment parity
		preFragmentsX = [] #-- x-momentum vector
		preFragmentsY = [] #-- y-momentum vector
		preFragmentsZ = [] #-- z-momentum vector

		c=0
		while True:
			c+=1
			if (hasMaxNumber and c>maxNumberEvents):
				break
			line = f.readline()
			if (len(line)==0):
				break

			line = line.split()
			Z.append(int(line[0]))
			A.append(int(line[1]))
			KE.append(float(line[2]))
			U.append(float(line[3]))
			J.append(float(line[4]))
			P.append(int(line[5]))
			preFragmentsX.append(float(line[6]))
			preFragmentsY.append(float(line[7]))
			preFragmentsZ.append(float(line[8]))
			
		f.close()

		data = np.dstack((A,Z,U,J,P,KE,preFragmentsX,preFragmentsY,preFragmentsZ))

		data = data[0,:,:]

		return (data)

	#################################################################
	#-- get all of the variables - instead of direct access		#
	#################################################################

	def getFissionYields(self):
		"""Returns a list with the full yields"""
		return (self.yields)

	def getNumberFragments(self):
		"""Returns the total number of fission fragments"""
		return (self.numberFragments)

	def getNumberEvents(self):
		"""Returns the number of simulated events"""
		return (self.numberEvents)

	def getA(self):
		"""Returns the mass (A) for each fragment"""
		return (self.A)

	def getZ(self):
		"""Returns the charge (Z) for each fragment"""
		return (self.Z)

	def getN(self):
		"""Returns the number of neutrons (N) for each fragment"""
		return (self.N)

	def getU(self):
		"""Returns the excitation energy (U) for each fragment"""
		return (self.U)

	def getJ(self):
		"""Returns the spin (J) of each fragment"""
		return (self.J)

	def getP(self):
		"""Returns the parity (P) of each fragment"""
		return (self.P)

	def getKE(self):
		"""Returns the kinetic energy for each fragment"""
		return (self.KE)

	def getTXE(self):
		"""Returns the total excitation energy for each event"""
		return (self.TXE)
	
	def getTKE(self):
		"""Returns the total kinetic energy for each event"""
		return (self.TKE)

	def getFragmentMomentum(self):
		"""Returns an array of momenta vectors for each fission fragment"""
		return (self.pFF)

	# quantities for the light fragments

	def getALF(self):
		"""Returns the masses of the light fragments"""
		return (self.Al)

	def getZLF(self):
		"""Returns the charges of the light fragments"""
		return (self.Zl)
	
	def getULF(self):
		"""Returns the excitation energies of the light fragments"""
		return (self.Ul)

	def getJLF(self):
		"""Returns the spins of the light fragments"""
		return (self.Jl)

	def getPLF(self):
		"""Returns the parities of the light fragments"""
		return (self.Pl)

	def getKELF(self):
		"""Returns the kinetic energies of the light fragments"""
		return (self.KEl)

	# quantities for the heavy fragments
	
	def getAHF(self):
		"""Returns the masses of the heavy fragments"""
		return (self.Ah)

	def getZHF(self):
		"""Returns the charges of the heavy fragments"""
		return (self.Zh)

	def getUHF(self):
		"""Returns the excitation energies of the heavy fragments"""
		return (self.Uh)

	def getJHF(self):
		"""Returns the spins of the heavy fragments"""
		return (self.Jh)

	def getPHF(self):
		"""Returns the parities of the light fragments"""
		return (self.Ph)

	def getKEHF(self):
		"""Returns the kinetic energies of the heavy fragments"""
		return (self.KEh)

	#################################################################
	#-- Summary Table for ipython notebook				#
	#################################################################

	def summaryTable (self):
		"""Returns a table for use in an ipython notebook which gives summary information about the CGMF simulation (for the fragments)"""

		table = ListTable()

		table.append(['','All Fragments','Light Fragments','Heavy Fragments'])

		#-- A, Z, TXE, TKE, J, U, pi
		table.append(['A',"{0:5.2f}".format(np.mean(self.A)),"{0:5.2f}".format(np.mean(self.Al)),"{0:5.2f}".format(np.mean(self.Ah))])
		table.append(['Z',"{0:5.2f}".format(np.mean(self.Z)),"{0:5.2f}".format(np.mean(self.Zl)),"{0:5.2f}".format(np.mean(self.Zh))])
		table.append(['TXE / U (MeV)',"{0:5.2f}".format(np.mean(self.TXE)),"{0:5.2f}".format(np.mean(self.Ul)),"{0:5.2f}".format(np.mean(self.Uh))])
		table.append(['TKE / KE (MeV)',"{0:5.2f}".format(np.mean(self.TKE)),"{0:5.2f}".format(np.mean(self.KEl)),"{0:5.2f}".format(np.mean(self.KEh))])
		table.append(['J ($\hbar$)',"{0:5.2f}".format(np.mean(self.J)),"{0:5.2f}".format(np.mean(self.Jl)),"{0:5.2f}".format(np.mean(self.Jh))])
		table.append(['parity',"{0:5.2f}".format(np.mean(self.P)),"{0:5.2f}".format(np.mean(self.Pl)),"{0:5.2f}".format(np.mean(self.Ph))])

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
			
