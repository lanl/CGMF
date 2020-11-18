import json
import numpy as np
import matplotlib.pylab as plt
import matplotlib.figure as figure

from pylab import *


#lib = ['ENDFB71', 'ENDFB80', 'JEFF3T3', 'JENDL40']
#labels = ["ENDF/B-VII.1", "ENDF/B-VIII.0", "JEFF-3.3T4", "JENDL-4.0"]
#ylabels = {1:'(n,tot)', 2:'(n,el)', 4:'(n,inel)', 16:'(n,2n)', 17:'(n,3n)', 18:'(n,f)', 102:'(n,$\gamma$)', 452:'Total Neutron Multiplicity', 456:'Prompt Neutron Multiplicity'}
#MTnumbers = {1:'n,tot', 2:'n,el', 4:'inelastic', 16:'n,2n', 17:'n,3n', 18:'n,f', 102: 'capture', 452:'nutot', 456:'nubar'}


################################################################################
#-- Read CGMF output history file with initial yields only
################################################################################
def readFragmentYieldsFromCGMF (filename):
	data = np.loadtxt(filename)
	return (data)



################################################################################
# old 'develop' format for now
################################################################################
def readHistoryFileFromCGMF (filename, **keyword_parameters):
	
	if ('fmt' in keyword_parameters):
		fmt = keyword_parameters['fmt']
	else:
		fmt = 'old'

	f = open (filename)

	A  = [] #-- fragment masses
	Z  = [] #-- fragment charges
	U  = [] #-- fragment initial excitation energies
	J  = [] #-- fragment initial angular momentum (hbar)
	P  = [] #-- fragment initial parity

	KEpre =  [] #-- pre-neutron emission fragment kinetic energy (MeV)
	KEpost = [] #-- post-neutron emission fragment kinetic energy (MeV)

	nmult = [] #-- neutron multiplicity
	gmult = [] #-- gamma multiplicity

	cmNeutronDircos    = [] #-- neutron directional cosines in center-of-mass frame of the fragment
	cmNeutronEnergies  = [] #-- neutron energies in CM frame (MeV)
	labNeutronDircos   = [] #-- neutron directional cosines in LAB frame
	labNeutronEnergies = [] #-- neutron energies in LAB frame (MeV)

	cmGammaDircos    = [] #-- gamma directional cosines in CM frame
	cmGammaEnergies  = [] #-- gamma energies in CM frame (MeV)
	labGammaDircos   = [] #-- gamma directional cosines in LAB frame
	labGammaEnergies = [] #-- gamma energies in LAB frame (MeV)

	preFragments  = [] #-- pre-neutron emission fragment momentum vectors
	postFragments = [] #-- post-neutron emission fragment momentum vectors

	preFragmentsX = []
	preFragmentsY = []
	preFragmentsZ = []

	postFragmentsX = []
	postFragmentsY = []
	postFragmentsZ = []

	##############################################################################

	line = f.readline()
	if (line.find("#")!=0):
		sys.exit("ERROR: FIRST LINE OF OUTPUT FILE SHOULD CONTAIN '#'")
	else:
		data = line[1:].split()
		ZAIDc   = int(data[0])
		Zc = int(ZAIDc/1000)
		Ac = int(ZAIDc)-1000*Zc
		Asym = int(float(Ac)/2.0)
		Einc    = float(data[1])
#    nevents = int(data[2])
#    alpha   = float(data[3])

	c=0
	index=0
	nfragments=0
	while True:
		c+=1
		line = f.readline()
		if (len(line)==0):
			break
		
		#		if (c>4):
		#	break
	
		data = line.split()
		nfragments+=1
			
		A.append(int(data[0]))
		Z.append(int(data[1]))
		U.append(float(data[2]))
		J.append(float(data[3]))
		P.append(int(data[4]))
		KEpre.append(float(data[5]))
		KEpost.append(float(data[6]))
		
		nn=int(data[7])
		ng=int(data[8])
		
		nmult.append(nn)
		gmult.append(ng)
		
		#-----------------------------------------------------------
		if (fmt is 'new'):
		#-----------------------------------------------------------

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

			if (ng>0):
				data=f.readline().split()
				for i in range(ng):
					Dg.append(np.array([float(data[0+i*4]),float(data[1+i*4]),float(data[2+i*4])]))
					Eg.append(float(data[3+i*4]))

			#-- store in final list for output
			
			cmNeutronEnergies.append(cmEn)
			cmNeutronDircos.append(cmDn)
			labNeutronEnergies.append(labEn)
			labNeutronDircos.append(labDn)
											 
			cmGammaEnergies.append(Eg)
			cmGammaDircos.append(Dg)
			labGammaEnergies.append(Eg)
			labGammaDircos.append(Dg)

			preFragmentsX.append (pfx)
			preFragmentsY.append (pfy)
			preFragmentsZ.append (pfz)

			postFragmentsX.append (pfx2)
			postFragmentsY.append (pfy2)
			postFragmentsZ.append (pfz2)

#			preFragments.append  (preFF)
#			postFragments.append (postFF)


		#-----------------------------------------------------------
		else:
		#-----------------------------------------------------------
		
			#-- neutrons -----------------------------------------------
			ne_cm  = [] # energies
			#    nv_cm  = [] # vectors
			ne_lab = [] # energies
			#    nv_lab = [] # vectors

			if (nn!=0): # neutrons
				line = f.readline()
				data = line.split()
				for i in range(nn):
					ne_cm.append(float(data[4*i]))
					#nv_cm.append(np.array([float(data[8*i]),float(data[8*i+1]),float(data[8*i+2])]))
					ne_lab.append(float(data[4*i+2]))
					#nv_lab.append(np.array([float(data[8*i+4]),float(data[8*i+5]),float(data[8*i+6])]))

			cmNeutronEnergies.append(ne_cm)
			#    cmNeutronDircos.append(nv_cm)
			labNeutronEnergies.append(ne_lab)
			#    labNeutronDircos.append(nv_lab)
		
			#-- gammas -------------------------------------------------
			ge_cm  = [] # energies
			#    gv_cm  = [] # vectors
			ge_lab = [] # energies
			#    gv_lab = [] # vectors
				
			if (ng!=0):
				line = f.readline()
				data = line.split()
				for i in range(ng):
					ge_cm.append(float(data[4*i]))
					# gv_cm.append(np.array([float(data[8*i]),float(data[8*i+1]),float(data[8*i+2])]))
					ge_lab.append(float(data[4*i+2]))
					# gv_lab.append(np.array([float(data[8*i+4]),float(data[8*i+5]),float(data[8*i+6])]))

			cmGammaEnergies.append(ge_cm)
			#    cmGammaDircos.append(gv_cm)
			labGammaEnergies.append(ge_lab)
			#		labGammaDircos.append(gv_lab)

		#    if (c>=100000):
		#        break

		#-----------------------------------------------------------
		# end If new format
		#-----------------------------------------------------------

	nevents = int(nfragments/2)
	#	print (nevents)

	if (fmt is 'new'):
		data = np.dstack((A,Z,U,J,P,KEpre,nmult,gmult,cmNeutronEnergies,labNeutronEnergies,cmGammaEnergies,labGammaEnergies,preFragmentsX,preFragmentsY,preFragmentsZ,postFragmentsX,postFragmentsY,postFragmentsZ,labNeutronDircos))
	else:
		data = np.dstack((A,Z,U,J,P,KEpre,nmult,gmult,cmNeutronEnergies,labNeutronEnergies,cmGammaEnergies,labGammaEnergies))

	data = data[0,:,:]

	return (data)






################################################################################
#-- Plot yields from CGMF history file
################################################################################
def plotYields (yields, quantity):

	l=len(yields)
	q=list()
	
	xlabels = ('Fission Fragment Mass (u)', 'Fission Fragment Charge')

	if (quantity=="YA"):
		for i in range(l):
			q.append(yields[i][:,1])
		xlabel = "Fission Fragment Mass (u)"
	elif (quantity=="YZ"):
		for i in range(l):
			q.append(yields[i][:,0])
		xlabel = "Fission Fragment Charge"
	elif (quantity=="YUi"):
		for i in range(l):
			q.append(yields[i][:,3])
		xlabel = "Initial Fragment Excitation Energy (MeV)"
	elif (quantity=="YKE"):
		for i in range(l):
			q.append(yields[i][:,2])
		xlabel = "Initial Fragment Kinetic Energy (MeV)"
	elif (quantity=="YJi"):
		for i in range(l):
			q.append(yields[i][:,4])
		xlabel = r"Initial Angular Momentum ($\hbar$)"
	else:
		return -1

	fig=figure(figsize=(14,6))
	plt.subplot(1,2,1)
	for i in range(l):
		qbins = np.arange(min(q[i]),max(q[i]))
		hq, binEdges = np.histogram(q[i],bins=qbins,normed=True)
		binCenters=0.5*(binEdges[1:]+binEdges[:-1])
		plt.step(binCenters,hq)
	#cor.plotExperimentalData(exp,"YA",author="Gook")
	plt.xlabel(xlabel)

	plt.subplot(1,2,2)
	for i in range(l):
		qbins = np.arange(min(q[i]),max(q[i]))
		hq, binEdges = np.histogram(q[i],bins=qbins,normed=True)
		binCenters=0.5*(binEdges[1:]+binEdges[:-1])
		plt.step(binCenters,hq)
	plt.xlabel(xlabel)
	plt.yscale('log')
	
	plt.show()
	
	return





################################################################################
#-- Read JSON-formatted data file
################################################################################
def readJSONDataFile (filename):
	"""Read a JSON-formatted data file
		
		:returns: list of JSON entries
		"""
	with open(filename) as jsonFile:
		jsonData = json.load(jsonFile)
		
		exp = list()
		c=0
		for item in jsonData.get("entries"):
			exp.append(item)
			c=c+1

	return exp

################################################################################
#-- Plot experimental data sets corresponding to a particular quantity
################################################################################
def plotExperimentalData (experiments, quantity, **keyword_parameters):
	"""Plot experimental data on current figure.
		
		:param experiments: list of experiments that can be plotted
		:param quantity:    type of data,e.g., 'YA', 'Pnu'
		:param **keyword_parameters: optional parameters that can be used to modify
		the plot. 'author' can be used to specify a particular set of data corresponding
		to just one author.
		'format' can be used to modify the plotting format. By default, 'ko--' is used.
		'renorm'
		
		"""
	if ('format' in keyword_parameters):
		fmtplot = keyword_parameters['format']
	else:
		fmtplot = 'ko--'

	checkAuthor = False
	if ('author' in keyword_parameters):
		author = keyword_parameters['author']
		checkAuthor = True
	for exp in experiments:
		if (exp['quantity']==quantity):
			if (checkAuthor):
				if (author not in exp['authors']):
					continue
			data = np.asarray(exp['data'])
			if ('renorm' in exp):
				coef = exp['renorm']
			else:
				coef = 1.0
			if ('label' in keyword_parameters):
				labelplot = keyword_parameters['label']
			else:
				labelplot = exp['label']
			
			datatype = exp['fmt']
			if (datatype=='xminxmaxydy'):
				plt.errorbar((data[:,0]+data[:,1])/2.0,data[:,2]*coef,xerr=(data[:,1]-data[:,0])/2.0,yerr=data[:,3]*coef,fmt=fmtplot,alpha=0.5,label=labelplot)
			elif (datatype=='xdxydy'):
				plt.errorbar(data[:,0],data[:,2],xerr=data[:,1],yerr=data[:,3],fmt=fmtplot,alpha=0.5,label=labelplot)
			elif (datatype=='xmid xlow xhigh ymid ylow yhigh'):
				xbot = (data[:,0]-data[:,1])*1e-6
				xtop = (data[:,2]-data[:,0])*1e-6
				ybot = (data[:,3]-data[:,4])
				ytop = (data[:,5]-data[:,3])
				plt.errorbar(data[:,0]*1e-6,data[:,3],xerr=(xbot,xtop),yerr=(ybot,ytop),fmt=fmtplot,alpha=0.5,label=labelplot)
			elif (datatype=='xy'):
				plt.plot(data[:,0],data[:,1],fmtplot,alpha=0.5,label=labelplot)
			elif (datatype=='xydy'):
				x = data[:,0]
				y = data[:,1]
				z = data[:,2]
				if (exp['units-dy']=='percent'):
					plt.errorbar(x,y*coef,yerr=y*z/100.0*coef,fmt=fmtplot,alpha=0.5,label=labelplot)
				else:
					plt.errorbar(x,y*coef,yerr=z*coef,fmt=fmtplot,alpha=0.5,label=labelplot)
			elif (datatype=='xdx1dx2ydy1dy2'):
				x = data[:,0]
				dx = data[:,1]
				y = data[:,3]
				dy = data[:,4]
				plt.errorbar(x,y,xerr=dx,yerr=dy,fmt=fmtplot,alpha=0.5,label=labelplot)


################################################################################
#-- Plot (xs,cov) for a given reaction channel.
################################################################################
def plotReaction(quantity, **keyword_parameters):
	
	exp = readJSONDataFile ("../experiments/94239-exp.json")
	#plotAllEvaluatedCrossSections (evaluations, quantity)
	#	plotExperimentalData (exp, quantity)
	
	return


################################################################################
#-- Lists all experiments read in
################################################################################
def listExperimentalData (experiments, **keyword_parameters):
	if ('quantity' in keyword_parameters):
		q = keyword_parameters['quantity']
	else:
		q = 'any'
	for exp in experiments:
		if (q=='any' or exp['quantity']==q):
			print ("{0:10s} |  {1}, {2}".format(exp['quantity'], exp['authors'], exp['year']))

################################################################################
#-- Upload evaluated libraries as given in user input dictionary
################################################################################
def uploadEvaluations (evaluations):
	for k,v in evaluations.items():
		print ('\n** ',v['label'],' **\n')
		v['xs']  = readAllCrossSectionsFromENDF(v['file'])
		print ('>> MF=3,1,5 , ',v['xs'].keys())
		v['cov'] = readAllCovariancesFromENDF(v['file'])
		print ('>> MF=33,31,35 , ',v['cov'].keys())
	
	return

################################################################################
#-- List all evaluations contained in dictionary 'd'
################################################################################
def listEvaluations (evaluations):
	for k,v in evaluations.items():
		print (k,' : ',v['label'],', ',v['file'])
	return

################################################################################
#-- check if the matrix is symmetric
################################################################################
def checkSymmetric(a,tol=1e-8):
	return np.allclose(a, a.T, atol=tol)


################################################################################
# Plot the cross sections corresponding to MT channel in the given library
################################################################################
def plotEvaluatedCrossSections (evaluations, library, MT, **keyword_parameters):
	
	#-- option to include standard deviations on the plot
	if ('unc' in keyword_parameters):
		unc = keyword_parameters['unc']
	else:
		unc = False
	
	#-- option to change color of plot and uncertainty band (if present)
	if ('fmt' in keyword_parameters):
		fmtplot = keyword_parameters['fmt']
	else:
		fmtplot = "k-"

	#-- order of layering
	if ('zorder' in keyword_parameters):
		order = keyword_parameters['zorder']
	else:
		order=0

	for k,v in evaluations.items():
		if (k==library):
			x = v['xs'][MT][0]*1e-6
			y = v['xs'][MT][1]
			title = v['label']
			title = str.replace(title,"\\$","$")
			title = str.replace(title,"\_","_")
			title = r'{0}'.format(title)
			if (unc):
				xstd=v['cov'][MT][1]*1e-6
				ystd = v['cov'][MT][2]
				ystdint = np.interp(x,xstd,ystd)
				ax=plt.gca()
				base_line, = ax.plot(x,y,fmtplot)
				plt.plot(x,y,label=title,color=base_line.get_color(),zorder=order)
				plt.fill_between(x,y*(1.0-ystdint/100.0),y*(1.0+ystdint/100.0),alpha=0.5,facecolor=base_line.get_color(),zorder=order)
			else:
				plt.plot(x, y, fmtplot, label=title, zorder=order)
	plt.xlabel("Neutron Energy (MeV)")
	if (MT==456 or MT==452):
		plt.ylabel(ylabels[MT]+' (n/f)')
	else:
		plt.ylabel(ylabels[MT]+' Cross Sections (b)')
	return

################################################################################
# Plot the cross sections corresponding to MT channel for all libraries
# contained in 'evaluations
################################################################################
def plotAllEvaluatedCrossSections (evaluations, MT):
	colors=['r','g','b','m','c','k']
	i=0
	for k,v in evaluations.items():
		plotEvaluatedCrossSections (evaluations, v['name'], MT, fmt=colors[i])
		i+=1
	return



################################################################################
# Plot the evaluated PFNS in the 'library' and for the incident neutron energy
# 'Einc' given in MeV.
################################################################################
def plotEvaluatedPFNS (evaluations, library, Einc, **keyword_parameters):
	MT=518
	if (Einc==0.0):
		Einc=1e-11
	for k,v in evaluations.items():
		if (k==library):
			En   = v['xs'][MT][0]
			Eout = v['xs'][MT][1]
			PFNS = v['xs'][MT][2]

			for i in range(size(En)):
				if (En[i]>Einc):
					break

			eout = Eout[i]
			pfns = PFNS[i]

			title = v['label']
			title = str.replace(title,"\\$","$")
			title = str.replace(title,"\_","_")
			title = r'{0}'.format(title)

			plt.xscale('log')
			plt.xlim(1e-3,20)
			plt.xlabel('Outgoing Neutron Energy (MeV)')
			plt.ylabel(r'Prompt Fission Neutron Spectrum ($n$/MeV)')

			plt.plot(eout,pfns,label=title)

	return


################################################################################
# Plot the evaluated PFNS correlation matrix in the 'library' and for the
#	incident neutron energy 'Einc' given in MeV.
################################################################################
def plotEvaluatedPFNSCorrelationMatrix (evaluations, library, Einc, **keyword_parameters):

	for k,v in evaluations.items():
		if (k==library):
			d=readPFNSCovFromENDF (v['file'],Einc=Einc)
			x=d['ebins']
			corr=d['corr']
			if (x.any()!=0.0):
				plt.xscale('log')
				plt.yscale('log')
				imgplot = plt.pcolormesh(x,x,corr,vmin=-1,vmax=1)
				plt.colorbar()
				plt.xlabel("Outgoing Neutron Energy (MeV)")
				plt.ylabel("Outgoing Neutron Energy (MeV)")
				plt.xlim(1e-5,20)
				plt.ylim(1e-5,20)
	return




################################################################################
# Plot evaluated PFNS from a suite of 'evaluations' for a particular incident
# energy 'Einc' given in MeV.
################################################################################
def plotAllEvaluatedPFNS (evaluations, Einc):
	colors=['r','g','b','m','c','k']
	i=0
	for k,v in evaluations.items():
		plotEvaluatedPFNS (evaluations, v['name'], Einc, fmt=colors[i])
		i+=1
	lg=plt.legend()
	return



################################################################################
# Plot (x,dy/y*100) relative standard deviations for all evaluations for the MT
# reaction channel.
################################################################################
def plotEvaluatedUncertainties (evaluations, library, MT, **keyword_parameters):
	#-- option to change color of plot and uncertainty band (if present)
	if ('fmt' in keyword_parameters):
		fmtplot = keyword_parameters['fmt']
	else:
		fmtplot = "k-"
	for k,v in evaluations.items():
		if (k==library):
			energies = v['cov'][MT][0]*1e-6
			emid=(energies[:-1]+energies[1:])/2.0
			title = v['label']
			title = str.replace(title,"\\$","$")
			title = str.replace(title,"\_","_")
			title = r'{0}'.format(title)
			if (size(emid)!=size(v['cov'][MT][2])):
				plt.step(0.0, 0.0, label=title) #-- to keep same color increment scheme throughout
			else:
				plt.step(emid,v['cov'][MT][2],color=fmtplot[0],label=title) #,color=fmtplot)
	plt.xlabel("Neutron Energy (MeV)")
	plt.ylabel(ylabels[MT]+' Uncertainties (%)')
	return

################################################################################
# Plot the relative standard deviations for all evaluations for the MT reaction
# channel.
################################################################################
def plotAllEvaluatedUncertainties (evaluations, MT):
	colors=['r','g','b','m','c','k']
	i=0
	for k,v in evaluations.items():
		plotEvaluatedUncertainties (evaluations, v['name'], MT, fmt=colors[i])
		i+=1
	return


################################################################################
# Plot the correlation matrix corresponding to a particular MT reaction channel
# and a given 'library'.
################################################################################
def plotEvaluatedCorrelationMatrix (evaluations, library, MT):
	for k,v in evaluations.items():
		if (k==library):
			x = v['cov'][MT][0]*1e-6
			z = v['cov'][MT][4]
			if (x.any()!=0.0):
				imgplot = plt.pcolormesh(x,x,z,vmin=-1,vmax=1)
				plt.colorbar()
				plt.xlabel("Neutron Energy (MeV)")
				plt.ylabel("Neutron Energy (MeV)")
				plt.xlim(min(x),min(max(x),20.0))
				plt.ylim(min(x),min(max(x),20.0))
			break
	return


################################################################################
# Create a figure with a summary of plots for a given reaction from a given
# library.
################################################################################
def plotReaction (evaluations, library, experiments, MT, **keyword_parameters):
	"""Produce a suite of plots related to the MT reaction channel in 'library'.
		
		-- Arguments --
		
		evaluations: a dictionary of evaluated libraries
		library:     the name of the library containing the reaction to be plotted
		MT:          the MT number corresponding to the reaction to be plotted.
		\
		"""

	#-- option to produce/save a figure
	if ('fig' in keyword_parameters):
		figname = keyword_parameters['fig']
	else:
		figname = ''

	for k,v in evaluations.items():
		if (k==library):
			en = v['xs'][MT][0]*1e-6
			xs = v['xs'][MT][1]
				
			x = v['cov'][MT][0][:-1]*1e-6
			y = v['cov'][MT][1]
			s = v['cov'][MT][2]
			
			cov  = v['cov'][MT][3]
			corr = v['cov'][MT][4]
			
			if (max(x)>20.0): # truncate plot at 20.0 MeV
				i=np.where(x>=20.0)[0][0]+1
				x=x[:i]
				y=y[:i]
				s=s[:i]
				cov  = cov[:i,:i]
				corr = corr[:i,:i]

			fig=figure(figsize=(14,12))
			
			ax1=plt.subplot(2,2,1)
			plotEvaluatedCrossSections (evaluations, k, MT, unc=True, fmt="r")
			lg=plt.legend(loc=2)
			lg.draw_frame(False)

			ax2=plt.subplot(2,2,2,sharex=ax1,sharey=ax1)
			plotExperimentalData (experiments, MTnumbers[MT])
			plotEvaluatedCrossSections (evaluations, k, MT, unc=True, fmt="r")
			lg=plt.legend(fontsize=16)
			lg.draw_frame(False)
			if (MT==102):
				plt.xscale('log')
				plt.yscale('log')
				plt.xlim(0.01,20.0)
				plt.ylim(0.01,10.0)
			elif (MT==1):
				plt.xlim(0.0,20.0)
				plt.xscale('log')
				plt.ylim(5.0,20.0)
				ax1.yaxis.set_major_locator(MultipleLocator(5))
				ax1.yaxis.set_minor_locator(MultipleLocator(1))
			elif (MT==18):
				plt.xlim(0.001, 20.0)
				plt.ylim(1.0, 3.0)
				ax1.xaxis.set_major_locator(MultipleLocator(5))
				ax1.xaxis.set_minor_locator(MultipleLocator(1))
				ax1.yaxis.set_major_locator(MultipleLocator(0.5))
				ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
			else:
				plt.xlim(0.0, 20.0)
				plt.ylim(0.001, 1.0)
				ax1.xaxis.set_major_locator(MultipleLocator(5))
				ax1.xaxis.set_minor_locator(MultipleLocator(1))
				ax1.yaxis.set_major_locator(MultipleLocator(0.2))
				ax1.yaxis.set_minor_locator(MultipleLocator(0.05))

			ax3=plt.subplot(2,2,3,sharex=ax1)
			plt.xlim(max(min(x),0.01), 20.0)
			plt.ylim(0.5*min(s),2.0*max(s))
			if (MT==18):
				plt.ylim(0.0,3.0)
			plotEvaluatedUncertainties (evaluations, k, MT, fmt="r")
			ax3.yaxis.set_ticks_position('both')
			lg=plt.legend()
			lg.draw_frame(False)
			if (MT==1 or MT==18):
				ax3.yaxis.set_major_locator(MultipleLocator(1))
				ax3.yaxis.set_minor_locator(MultipleLocator(0.25))
			else:
				ax3.yaxis.set_major_locator(MultipleLocator(20))
				ax3.yaxis.set_minor_locator(MultipleLocator(5))
			
			ax4=plt.subplot(2,2,4,aspect='equal')
			plotEvaluatedCorrelationMatrix(evaluations, k, MT)
			plt.xlim(min(x),max(x))
			plt.ylim(min(x),max(x))
			if (MT==102):
				plt.xscale('log')
				plt.yscale('log')
				plt.xlim(0.001,20)
				plt.ylim(0.001,20)
			else:
				ax4.xaxis.set_major_locator(MultipleLocator(5))
				ax4.xaxis.set_minor_locator(MultipleLocator(1))
				ax4.yaxis.set_major_locator(MultipleLocator(5.0))
				ax4.yaxis.set_minor_locator(MultipleLocator(1.0))
			
			fig.tight_layout()
			fig.subplots_adjust(top=0.92)
			fig.suptitle("$^{239}$Pu "+ylabels[MT],fontsize=30)
			
			plt.show()

			if (figname!=''):
				fig.savefig(figname)

################################################################################
################################################################################
def readAllCrossSectionsFromENDF (filename):
	mf=3
	mt=[1,2,4,16,17,18,102]
	d={}
	for i,mti in enumerate(mt):
		x, y= readCrossSectionsFromENDF(filename,3,mti)
		d[mti] = x,y
	mf=1
	mt=[452,456]
	for i,mti in enumerate(mt):
		x, y= readCrossSectionsFromENDF(filename,1,mti)
		d[mti] = x,y
	mf=5
	mt=[518]
	for i,mti in enumerate(mt):
		einc,x,y = readPFNS(filename)
		d[mti] = einc, x, y
	return (d)

################################################################################
################################################################################
def readCrossSectionsFromENDF (filename, mf, mt):
	
	f=open(filename,'r')
	
	data = f.readlines()
	section = []
	
	#-- find (MF,MT) section
	for line in data:
		#		mat = line[66:70]
		mf0 = int(line[70:72])
		mt0 = int(line[72:75])
		if (mf0==mf and mt0==mt):
			section.append(line)

	f.close()
	
	ZA,AWR,d,d,d,d = readHEAD(section[0])
	QM,QI,d,LR,NR,NP,Eint,sigma = readTAB1(section[1:])

	return (Eint,sigma)

################################################################################
# Return a list 'data' that contains the section (MF,MT)
################################################################################
def getSectionMFMT (filename, mf, mt):
	f=open(filename,'r')
	data=f.readlines()
	section=[]
	for line in data:
		#		mat = line[66:70]
		mf0 = int(line[70:72])
		mt0 = int(line[72:75])
		if (mf0==mf and mt0==mt):
			section.append(line)
	f.close()
	return section


################################################################################
################################################################################
def readPnuFromENDF (filename, particleType):
	if (particleType!='neutron' and particleType!='gamma'):
		return -1
	section = getSectionMFMT (filename, 6, 18)
	ZA,AWR,JP,LCT,NK,d = readHEAD(section[0])
	section=section[1:]
	Pnu  = [] # neutrons
	Pnug = [] # gammas
	for i in range(NK):
		ZAP,AWP,LIP,LAW,NR,NP,Eint,yi = readTAB1(section)
		if (float(ZAP)==1.0):
			EincNeutrons = Eint
			Pnu.append(yi)
		elif (float(ZAP)==0.0):
			EincGammas = Eint
			Pnug.append(yi)
		if (NP%3==0):
			nl=NP//3+2
		else:
			nl=NP//3+3
		section = section[nl:]
	if (particleType=='neutron'):
		en, nu = readCrossSectionsFromENDF (filename, 1, 456)
		nubar = np.interp(EincNeutrons,en,nu)
		return EincNeutrons, np.asarray(nubar*Pnu)
	else:
		en, nug = readCrossSectionsFromENDF (filename, 12, 18)
		nubarg = np.interp(EincGammas,en,nug)
	return EincGammas, np.asarray(nubarg*Pnug)

################################################################################
################################################################################
def readAllCovariancesFromENDF (filename):
	mf=33
	mt=[1,4,16,17,18,102]
	n=size(mt)
	d={}
	for i,mti in enumerate(mt):
		energies, ebins, std, cov, corr, covListElements = readCovFromENDF(filename,mf,mti)
		d[mti] = energies, ebins, std, cov, corr, covListElements
	mf=31
	mt=[452,456]
	for i,mti in enumerate(mt):
		energies, ebins, std, cov, corr, covListElements = readCovFromENDF(filename,mf,mti)
		d[mti] = energies, ebins, std, cov, corr, covListElements
	mf=35
	mt=[518]
	for i,mti in enumerate(mt):
		energies, ebins, std, cov, corr, covListElements = readCovFromENDF(filename,mf,18)
		d[mti] = energies, ebins, std, cov, corr, covListElements
	return (d)

################################################################################
################################################################################
def readCovFromENDF(filename, MF, MT):
	
	#	print ("readCovFromENDF: MF,MT= ", MF, MT)
	
	#-- setup default values
	n=1
	covd=np.zeros((n,n))
	corrd=np.zeros((n,n))
	energiesd=np.zeros(n)
	ebinsd=np.zeros(n-1)
	stdd=np.zeros(n)
	matrixElements=[]
	
	#-- also need to read corresponding cross section data
	if (MF==31):
		en, xs = readCrossSectionsFromENDF (filename, 1, MT)
	elif (MF==33):
		en, xs = readCrossSectionsFromENDF (filename, 3, MT)
	elif (MF==35):
		en, xs = readCrossSectionsFromENDF (filename, 5, MT)
	else:
		print ("Not implemented yet! - MF, MT = ", MF, MT)
		return (energiesd,ebinsd,stdd,covd,corrd,matrixElements)

	#-- special case of PFNS MF35, MT18
	if (MF==35 and MT==18):
		d = readPFNSCovFromENDF (filename)
		return (d['energies'], d['ebins'], d['std'], d['cov'], d['corr'], d['matrixElements'])

	f=open(filename,'r')
	
	data = f.readlines()
	section = []

	#-- find (MF,MT) section
	found=False
	for line in data:
		MAT = line[66:70]
		MF0 = int(line[70:72])
		MT0 = int(line[72:75])
		if (MF0==MF and MT0==MT):
			section.append(line)
			found=True

	if (not found):
		print("*** Section (", MF, ",",MT,") not found!")
		return (energiesd,ebinsd,stdd,covd,corrd,matrixElements)
	
	ZA,AWR,d,MTL,d,NL = readHEAD(section[0])

	section=section[1:]
	covList = []
	#-- loop over NL subsections -----------------------
	for iNL in range(NL):

#		print ("Subsection # ", iNL)
		XMF1,XLFS1,MAT1,MT1,NC,NI = readCONT(section[0])

#		print ("MAT1, MT1, NC, NI = ", MAT1, MT1, NC, NI)

		if (NC!=0):
			print ("*** Not yet implemented - NC!=0 ", NC, MF, MT)
			return (energiesd,ebinsd,stdd,covd,corrd,matrixElements)
		else:

			section=section[1:]
			#-- loop over NI sub-subsections -------------------------
			for iNI in range(NI):

#				if (MT==18):
#					print ("- iNI: ", iNI)
#					print (section[0])

#				print (section[0])
				d,d,LS,LB,NT,NE,matrixElements = readLIST (section)
				#print ("-- iNI, LS, LB, NT, NE = ", iNI, LS,LB,NT,NE)

				# number of lines in section
				if (NT%6==0):
					nl=int(NT/6)+1
				else:
					nl=int(NT/6)+2

				if (MT!=MT1 or LS!=1):
					section = section[nl:]
					continue

				# initialize arrays
				cov=np.zeros((NE-1,NE-1))
				corr=np.zeros((NE-1,NE-1))
				ebins=np.zeros(NE-1)
				std=np.zeros(NE-1)

				#-----------------------------------------------------------------------
				if (LB<5):
					print (matrixElements)
					energies = matrixElements[::2]
					print (energies)

					#-- interpolate cross sections of (3,MT) on energy grid of (33,MT)
					xsint = np.interp(energies,en,xs)

					k=0
					for i in range(NE-1):
						for j in range(i,NE-1): #- upper symmetric
							cov[i,j]=matrixElements[k]*xsint[i]*xsint[j]
							k=k+1

				else:


					if (LB==6): #-- cross-reaction correlations (non-square matrix)
						d,d,d,LB,NT,NER,matrixElements = readLIST (section)
					elif (LS!=1 or LB!=5):
						print ("*** Not yet implemented - LS!=1 or 5 - MT=", MT, LS, LB)
						if (iNI==0):
							return (energiesd,ebinsd,stdd,covd,corrd,matrixElements)
						else:
							if (NT%6==0):
								nl=int(NT/6)+1
							else:
								nl=int(NT/6)+2
								section = section[nl:]
								break

					energies = np.asarray(matrixElements[0:NE])
					matrixElements = matrixElements[NE:]
					
					#-- interpolate cross sections of (3,MT) on energy grid of (33,MT)
					xsint = np.interp(energies,en,xs)

					k=0
					for i in range(NE-1):
						for j in range(i,NE-1): #- upper symmetric
							cov[i,j]=matrixElements[k]*xsint[i]*xsint[j]
							k=k+1

					#-----------------------------------------------------------------------

					#-- fill out symmetric lower part
					for i in range(NE-1):
						for j in range(i):
							cov[i,j]  = cov[j,i]

					#-- compute correlation matrix
					for i in range(NE-1):
						for j in range(NE-1):
							if (cov[i,i]==0 or cov[j,j]==0):
								corr[i,j]=0.0
							else:
								corr[i,j]=cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])
	
					if (not checkSymmetric(cov)):
						print ("*** Matrix not symmetric!")

					#-- compute relative standard deviations (%)
					for i in range(NE-1):
						if (xsint[i]==0.0):
							std[i]=0.0
						else:
							std[i]=np.sqrt(cov[i,i])/xsint[i]*100.0
						ebins[i]=0.5*(energies[i+1]+energies[i])

#				print ("last terms: ", matrixElements[-5:])
				#-- add sub-section data to list of covariance sub-sections
#				if (MT1==MT):
					#					print ("MT1=MT, first matrixElements: ", matrixElements[0:5])
					#print ("First terms cov: ", cov[0:5,0:5])
					#print (shape(cov))
					covList.append(
											 {
											 'energies':energies,
											 'ebins':ebins,
											 'std':std,
											 'cov':cov,
											 'corr':corr,
											 'matrixElements':matrixElements}
											 )

				section = section[nl:]
					#print ('covList: ', shape(covList))

	f.close()

#	print (shape(covList))
#	print ('energies: ', covList[0]['energies'])
	ic=0
	if (MT==18 and (("n-094_Pu_239.endf" in filename) or ("n-094_Pu_239-beta5.endf" in filename))):
		ic=1

	energies = covList[ic]['energies']
	ebins    = covList[ic]['ebins']
	std      = covList[ic]['std']
	cov      = covList[ic]['cov']
	corr     = covList[ic]['corr']
	matrixElements = covList[ic]['matrixElements']

	#-- for now, just return the first one (or second one for fission!) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TO CHANGE!
	#	return (covList[0])

	return (energies,ebins,std,cov,corr,matrixElements)



#-- read CONT records ----------------------------------------------------------
def readCONT (line):
	C1 = fix(line[0:11])
	C2 = fix(line[11:22])
	L1 = int(line[22:33])
	L2 = int(line[33:44])
	N1 = int(line[44:55])
	N2 = int(line[55:66])
	return (C1,C2,L1,L2,N1,N2)

#-- read HEAD records ----------------------------------------------------------
def readHEAD (line):
	ZA, AWR, L1, L2, N1, N2 = readCONT(line)
	return (ZA,AWR,L1,L2,N1,N2)

#-- read LIST records ----------------------------------------------------------
def readLIST (lines):
	line=lines[0]
	C1 = fix(line[0:11])
	C2 = fix(line[11:22])
	L1 = int(line[22:33])
	L2 = int(line[33:44])
	NPL = int(line[44:55])
	N2  = int(line[55:66])

	if (NPL%6==0):
		nlines=int(NPL/6)
	else:
		nlines=int(NPL/6)+1

	B=list()
	data=lines[1:]
	for i in range(nlines):
		line=data[i][0:66]
		l=""
		for j in range(6):
			l2=line[j*11:(j+1)*11]
			l2=l2.replace('+','e+'); l2=l2[0]+l2[1:].replace('-','e-')
			l=l+" "+l2
		r = list(map(float,l.split()))
		B.append(r)

	B = np.concatenate(B)
	B = B.flatten()

	return (C1,C2,L1,L2,NPL,N2,B)


#-- read TAB1 records ----------------------------------------------------------
def readTAB1 (lines):
	line=lines[0]
	C1 = fix(line[0:11])
	C2 = fix(line[11:22])
	L1 = int(line[22:33])
	L2 = int(line[33:44])
	NR = int(line[44:55])
	NP = int(line[55:66])
	
	x=np.empty(NP,dtype=float)
	y=np.empty(NP,dtype=float)
	
	
	il=1
	NBT=np.empty(NR,dtype=int)
	INT=np.empty(NR,dtype=int)
	line=lines[il]
	
	for n in range(NR):
		NBT[n] = int(line[n*11:(n+1)*11])
		INT[n] = int(line[(n+1)*11:(n+2)*11])
	
	nlines = NP//3+1
	c=0
	for il in range(2,nlines+2):
		data=lines[il]
		for i in range(0,min(3,NP-c)):
			x[c] = fix(data[2*i*11:(2*i+1)*11])
			y[c] = fix(data[(2*i+1)*11:(2*i+2)*11])
			c=c+1
		if (NP==c):
			break

	return (C1,C2,L1,L2,NR,NP,x,y)


#-- read TAB2 records ----------------------------------------------------------
def readTAB2 (lines):
	line=lines[0]
	C1 = fix(line[0:11])
	C2 = fix(line[11:22])
	L1 = int(line[22:33])
	L2 = int(line[33:44])
	NR = int(line[44:55])
	NZ = int(line[55:66])
	
	NBT=np.empty(NR,dtype=int)
	INT=np.empty(NR,dtype=int)
	line=lines[1]
	for n in range(NR):
		NBT[n] = int(line[n*11:(n+1)*11])
		INT[n] = int(line[(n+1)*11:(n+2)*11])
	
	return (C1,C2,L1,L2,NR,NZ,NBT,INT)


#-- read TEXT records ----------------------------------------------------------
def readTEXT(fh):
	line = fh.readline()
	HL = line[0:66]
	MAT = line[67:71]
	MF = line[71:73]
	MT = line[74:77]
	return (HL,MAT,MF,MT)



#-- fixes a number for missing 'e' in scientific notation
def fix(s):
	i=s.find('+')
	if (i>0):
		s = s[0:i]+"e"+s[i:]
	i=s.rfind('-')
	if (i>0):
		s = s[0:i]+"e"+s[i:]
	return s


#-- transforms a string of numbers to conform to ENDF formatting (not "e")
#-- the inverse of fix(s)
def fixENDF(s):
	s = s.replace("e+0","+")
	s = s.replace("e-0","-")
	s = s.replace("0e-1","-1") #<< need to replace any number before "e" by blank...
	return s


################################################################################
################################################################################
def writeCovToENDF (MAT, MF, MT, energies, xs, cov, filename, seqstart):

	f = open (filename,"w")

	n=size(energies)

	seq = seqstart #-- start sequence number here

	covElements=[]
	k=0
	for i in range(n-1):
		for j in range(i,n-1):
			if (xs[i]==0 or xs[j]==0):
				covElements.append(0.0)
			else:
				covElements.append(cov[i,j]/xs[i]/xs[j])
			k=k+1

	NT=n*(n+1)/2
	s = fixENDF("%13.6e%13.6e%11i%11i%11i%11i%4i%2i%3i%5i\n"%(0.0,0.0,1,5,NT,n,MAT,MF,MT,seq))
	f.write(s)
	seq+=1

	k=0
	c=0
	for i in range(n):
		c=c+1
		if (c>6):
			c=1
			f.write("%4i%2i%3i%5i\n"%(MAT,MF,MT,seq))
			seq+=1
		
		f.write(fixENDF("%13.6e"%(energies[i])))
#f.write('%-40s %6s %10s %2s\n' % (filename, type, size, modified))

	for i in range(n-1):
		for j in range(i,n-1):
			c=c+1
			if (c>6):
				c=1
				f.write("%4i%2i%3i%5i\n"%(MAT,MF,MT,seq))
				seq+=1
			f.write(fixENDF("%13.6e"%covElements[k]))
			k=k+1

	if (c<=6):
		for i in range(6-c):
			f.write("           ")
		f.write("%4i%2i%3i%5i\n"%(MAT,MF,MT,seq))
	f.close()

	return 0



################################################################################
################################################################################
def modifyCovarianceMatrix (en, xs, cov, scaling, Emin, coef, **keyword_parameters):
	n=size(cov[:,0])
	if (size(xs)!=(n+1) or size(en)!=(n+1)):
		print ("Size mismatch in modifyCovarianceMatrix!")
		return -1
	
	adhocCov  = np.zeros((n,n))
	adhocCorr = np.zeros((n,n))
	adhocStd  = np.zeros(n)

	std = np.zeros(n)
	for i in range(n):
		if (xs[i]==0.0):
			std[i]=0.0
		else:
			std[i]=np.sqrt(cov[i,i])/xs[i]*100.0

	#-- copy original matrix
	adhocCov=np.copy(cov)
	
	#-- overall scaling
	adhocCov=cov*scaling**2

	#-- set minimum std. dev. (%)
	if ('stdmin' in keyword_parameters):
		stdmin = keyword_parameters['stdmin']
		for i in range(n):
			if (std[i]<stdmin):
				std[i]=stdmin
				adhocCov[i,i] = (stdmin/100.0*xs[i])**2

	#-- add constant term in quadrature
	if ('quad' in keyword_parameters):
		quad = keyword_parameters['quad']
		adhocCov += (quad*np.ones((n,n)))

	#-- high-energy corrections
	#Emin*=1e6 # in eV
	coef*=1e-6
	for i in range(n):
		for j in range(n):
			if (en[i]>Emin and en[j]>Emin):
				adhocCov[i,j]+=np.sqrt(coef)*xs[i]*xs[j]*np.log(1.0+(0.5*(en[i]+en[j])-Emin)/10.0)

	#-- low-energy corrections
	if ('Emax' in keyword_parameters):
		Emax = keyword_parameters['Emax']
		coef2 = keyword_parameters['coef2']
		for i in range(n):
			for j in range(n):
				if (en[i]<Emax and en[j]<Emax):
					adhocCov[i,j]*=coef2**2

	for i in range(n):
		if (xs[i]==0.0):
			adhocStd[i]=0.0
		else:
			adhocStd[i]=np.sqrt(adhocCov[i,i])/xs[i]*100.0

	#-- correlation matrix
	adhocCorr = correlationMatrix(adhocCov)
	return (adhocStd,adhocCov,adhocCorr)


################################################################################
def correlationMatrix (cov):
	n=size(cov[:,0])
	corr=np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			if (cov[i,i]==0.0 or cov[j,j]==0.0):
				corr[i,j]=1.0
			else:
				corr[i,j]=cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])
	return (corr)





################################################################################
# Algorithm from Higham, 2000.

def _getAplus(A):
	eigval, eigvec = np.linalg.eig(A)
	Q = np.matrix(eigvec)
	xdiag = np.matrix(np.diag(np.maximum(eigval, 0)))
	return Q*xdiag*Q.T

def _getPs(A, W=None):
	W05 = np.matrix(W**.5)
	return  W05.I * _getAplus(W05 * A * W05) * W05.I

def _getPu(A, W=None):
	Aret = np.array(A.copy())
	Aret[W > 0] = np.array(W)[W > 0]
	return np.matrix(Aret)

def nearPD(A, nit=10):
	n = A.shape[0]
	W = np.identity(n)
	# W is the matrix used for the norm (assumed to be Identity matrix here)
	# the algorithm should work for any diagonal W
	deltaS = 0
	Yk = A.copy()
	for k in range(nit):
		Rk = Yk - deltaS
		Xk = _getPs(Rk, W=W)
		deltaS = Xk - Rk
		Yk = _getPu(Xk, W=W)
	return Yk



################################################################################
# Find the nearest positive-definite matrix to the matrix given in input
# Based on Higham's 1988 paper (https://doi.org/10.1016/0024-3795(88)90223-6 )
################################################################################
def nearestPD(A):
	"""Find the nearest positive-definite matrix to input
		
		A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
		credits [2].
		
		[1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
		
		[2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
		matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
		"""
	B = (A + A.T) / 2
	_, s, V = linalg.svd(B)
	H = np.dot(V.T, np.dot(np.diag(s), V))
	A2 = (B + H) / 2
	A3 = (A2 + A2.T) / 2
	
	if isPD(A3):
		return A3

	spacing = np.spacing(linalg.norm(A))
	# The above is different from [1]. It appears that MATLAB's `chol` Cholesky
	# decomposition will accept matrixes with exactly 0-eigenvalue, whereas
	# Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
	# for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
	# will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
	# the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
	# `spacing` will, for Gaussian random matrixes of small dimension, be on
	# othe order of 1e-16. In practice, both ways converge, as the unit test
	# below suggests.
	I = np.eye(A.shape[0])
	k = 1
	while not isPD(A3):
		mineig = np.min(np.real(linalg.eigvals(A3)))
		A3 += I * (-mineig * k**2 + spacing)
		k += 1
	
	return A3

#-- part of Higham's algorithm
def isPD(B):
	"""Returns true when input is positive-definite, via Cholesky"""
	try:
		_ = linalg.cholesky(B)
		return True
	except linalg.LinAlgError:
		return False

	if __name__ == '__main__':
		import numpy as np
		for i in xrange(10):
			for j in xrange(2, 100):
				A = np.random.randn(j, j)
				B = nearestPD(A)
				assert(isPD(B))
		print('unit test passed!')





################################################################################
#-- Read covariance matrix from a DeCE output, after ENDF read-out
################################################################################
def readCovFromDece(filename):
	energy=[]
	std=[]
	covelements=[]
	f=open(filename,'r')
	for line in f:
		if "#" in line:
			next
		else:
			l = line.split()
			energy.append(l[0])
			std.append(l[1])
			covelements.append(l[2:])
	std=np.asarray(std,float)
	NP = size(energy)
	NT=NP*(NP+1)//2
	corr = np.zeros((NP,NP))
	cov  = np.zeros((NP,NP))
	for i in range(NP):
		x = np.array(covelements[i])
		y = x.astype(np.int)
		cov[i,0:size(y)]=y/1000.0*std[i]
	
	for i in range(NP):
		for j in range(i+1,NP):
			cov[i,j]  = cov[j,i]

	for i in range(NP):
		for j in range(NP):
			if (cov[i,i]==0 or cov[j,j]==0):
				corr[i,j]=0.0
			else:
				corr[i,j]=cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])

	if (not checkSymmetric(cov)):
		print ("Matrix not symmetric!")
	
	f.close()
	return (cov,corr)



################################################################################
################################################################################
def plotENDF (filename, mf, mt, withCov=False):
	
	en, xs = readCrossSectionsFromENDF(filename,mf,mt)
	
	if (withCov):
		mfc=int('3'+str(mf))
		enc, std, cov, corr = readCovFromENDF(filename,mfc,mt)
	
	if (not withCov):
		
		plt.plot(en,xs)
		plt.show()
	
	else:
		
		plt.subplot(2,2,1)
		plt.plot(en,xs)
		
		with np.errstate(invalid='ignore'):
			relstd = std/np.interp(enc,en,xs)*100
		
		plt.subplot(2,2,2)
		plt.step(enc,relstd)
		
		plt.subplot(2,2,3)
		imshow(corr,vmin=-1,vmax=1)
		colorbar()
		
		plt.subplot(2,2,4)
		imshow(cov)
		colorbar()
		
		plt.tight_layout()
		plt.show()
	
		return



################################################################################
# Prompt Fission Neutron Spectrum Utilities
################################################################################

#-------------------------------------------------------------------------------
# read PFNS in section MF5, MT18
#-------------------------------------------------------------------------------
def readPFNS (filename):
	
	f=open(filename,'r')
	
	data = f.readlines()
	section = []
	
	#-- find (MF,MT) section
	for line in data:
		mf0 = int(line[70:72])
		mt0 = int(line[72:75])
		if (mf0==5 and mt0==18):
			section.append(line)

	f.close()
	
	ZA,AWR,d,d,d,d = readHEAD(section[0])
	QM,QI,d,LF,NR,NP,Eint,sigma = readTAB1(section[1:])
	
	#-- LF=1 (arbitrary tabulated function)
	if (LF!=1):
		print ("Error: cannot read LF different from 1!")
		return -1

	d,d,d,d,NR,NE,NBT,INT = readTAB2(section[4:])

	Einc = []
	Eout = []
	PFNS = []
	
	data = section[6:]
	for i in range(NE):
		d,incidentEnergy,d,d,NR,NF,X,Y=readTAB1(data)
		Einc.append(float(incidentEnergy)*1e-6)
		Eout.append(X*1e-6)
		PFNS.append(Y*1e6)
		# switch to next section
		nl=NF//3+2
		if (NF%3):
			nl=nl+1
		data=data[nl:]

	return (Einc,Eout,PFNS)


#-------------------------------------------------------------------------------
# read PFNS covariances in section MF35, MT18
#-------------------------------------------------------------------------------
def readPFNSCovFromENDF (filename, **keyword_parameters):
	
	f = open (filename)
	data = f.readlines()
	section = []
	
	MF=35
	MT=18

#	print ("Reading PFNS Covariance!")

	#-- Incident neutron energy (MeV)
	if ('Einc' in keyword_parameters):
		incidentEnergy = keyword_parameters['Einc']
	else:
		incidentEnergy = 0.0 #-- default

	if (incidentEnergy==0.0):
		incidentEnergy=1e-5

	#-- find (MF,MT) section
	found=False
	for line in data:
		MAT = line[66:70]
		MF0 = int(line[70:72])
		MT0 = int(line[72:75])
		if (MF0==MF and MT0==MT):
			section.append(line)
			found=True
		if (found and (MF0!=MF or MT0!=MT)):
			break # exit once the full section is read

	if (not found):
		print("*** Section (", MF, ",",MT,") not found!")
		return (energiesd,ebinsd,stdd,covd,corrd,matrixElements)

	#-- read PFNS values from (5,18)
	einc, eout, pfns = readPFNS (filename)
	for k in range(size(einc)):
		if (einc[k]>incidentEnergy):
			if (k>=1):
				k=k-1
			outgoingEnergies = eout[k]
			if (k==0):
				spectrum = pfns[k]
			else:
				spectrum = (pfns[k]-pfns[k-1])*(incidentEnergy-einc[k-1])/(einc[k]-einc[k-1])+pfns[k-1]
			break

	ZA,AWR,d,d,NK,d = readHEAD(section[0])

	section=section[1:]

	found=False
	#-- loop over NK sections (number of covariance blocks)
	for k in range(NK):

		E1, E2, LS, LB, NT, NE, matrixElements = readLIST (section)
		E1 = float(E1)*1e-6
		E2 = float(E2)*1e-6
		#		print (k, E1, E2, LS, LB, NT, NE)
		
		# number of lines in section
		if (NT%6==0):
			nl=int(NT/6)+1
		else:
			nl=int(NT/6)+2
		section = section[nl:]

		# check if this matrix block corresponds to the incident neutron energy
		# if not, continue search
		if (incidentEnergy<E1 or incidentEnergy>E2):
			continue

		# initialize arrays
		cov=np.zeros((NE-1,NE-1))
		corr=np.zeros((NE-1,NE-1))
		ebins=np.zeros(NE-1)
		std=np.zeros(NE-1)

		if (LB!=7): # absolute covariance matrix
			print ("Cannot read LB!=7 section for PFNS covariance matrix")
				
		energies = np.asarray(matrixElements[0:NE])
		energies *= 1e-6
		interpSpectrum = np.interp(energies,outgoingEnergies,spectrum)
		
		matrixElements = matrixElements[NE:]
		energyBins = np.zeros(NE-1)
		for i in range(NE-1):
			energyBins[i] = energies[i+1]-energies[i]
																							 
		k=0
		for i in range(NE-1):
			for j in range(i,NE-1): #- upper symmetric
				#cov[i,j]=matrixElements[k]/energyBins[i]/energyBins[j]
				if (interpSpectrum[i]==0.0 or interpSpectrum[j]==0.0):
					cov[i,j]=0.0
				else:
					cov[i,j]=matrixElements[k]/interpSpectrum[i]/interpSpectrum[j]/energyBins[i]/energyBins[j]
				#cov[i,j]=matrixElements[k]/spectrum[i]/spectrum[j]/energyBins[i]/energyBins[j]
				k=k+1

		#-- fill out symmetric lower part
		for i in range(NE-1):
			for j in range(i):
				cov[i,j]  = cov[j,i]
			
		#-- compute correlation matrix
		for i in range(NE-1):
			for j in range(NE-1):
				if (cov[i,i]==0 or cov[j,j]==0):
					corr[i,j]=0.0
				else:
					corr[i,j]=cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])

		if (not checkSymmetric(cov)):
			print ("*** Matrix not symmetric!")

		#-- compute relative standard deviations (%)
		for i in range(NE-1):
			std[i]=np.sqrt(cov[i,i])*100.0
			ebins[i]=0.5*(energies[i+1]+energies[i])
				
		#-- add sub-section data to list of covariance sub-sections
		d = {
			'energies':energies,
			'ebins':ebins,
			'std':std,
			'cov':cov,
			'corr':corr,
			'matrixElements':matrixElements
			}

		return (d)

	if (not found):
		d={}
		print ("not found")
		return(d)


#-------------------------------------------------------------------------------
# read PFNS from data extracted with DeCE: 'dece -f 5 -t 18 file.endf'
#-------------------------------------------------------------------------------
def readPFNSfromDece (filename):
	
	f=open(filename,'r')
	ne=15
	neout=643
	
	x = np.zeros((ne,neout))
	y = np.zeros((ne,neout))
	
	c=-1
	while (True):
		line = f.readline()
		if ('#            E' in line):
			E = float(line.split()[2])
			c=c+1
		if ('#           NF' in line):
			NF = int(line.split()[2])
			
			x[c,:]=0.0
			y[c,:]=0.0
			for i in range(NF):
				(x[c,i],y[c,i])=f.readline().split()
			
			x[c,:]*=1e-6
			y[c,:]*=1e6
			
			if (NF<neout):
				x[c,NF:neout]=x[c,NF-1]
				y[c,NF:neout]=0.0
		
		if (c>=ne):
			break

	return (x,y)

#-------------------------------------------------------------------------------
# Perform an eigenvalue decomposition of the covariance matrix associated with
# the mean values given in "spectrum", and returns a number (nsamples) of
# realizations of perturbed vectors.
#-------------------------------------------------------------------------------
def sampleCovMatrix (cov, spectrum, nsamples):
	
	s = linalg.eig(cov) # [1] diagonalize here
	eigenvalues  = s[0]
	eigenvectors = s[1]
	
	# print ne, eigenvalues.size (should be equal)
	ne = eigenvalues.size
	
	# [2] sample eigenvalue distributions using Gaussian distributions with
	# standard deviations of sqrt(eigenvalues)
	ns=nsamples
	
	# [3] sample PFNS distributions
	samples = np.empty([ns,ne])
	eps = np.empty([ne])
	for i in range(ns): # number of desired samples
		for j in range(ne):
			eps[j] = np.random.normal(0.0, sqrt(abs(eigenvalues[j])))
		for j in range(ne): # corresponds to number of points on outgoing energy grid
			samples[i][j]=spectrum[j]+sum(eigenvectors[j]*eps)

	return samples, eigenvalues

#-------------------------------------------------------------------------------
# Check if matrix is positive, definite.
#-------------------------------------------------------------------------------
def isPositiveSemiDefinite (cov):
	s = linalg.eig(cov)
	eigenvalues = s[0]
	if (np.all(eigenvalues>=0)):
		return 1
	else:
#		print (eigenvalues)
		return 0


#-------------------------------------------------------------------------------
# Check if matrix is symmetric
#-------------------------------------------------------------------------------
def isSymmetric (cov):
	if (np.allclose(cov.transpose(),cov,rtol=1e-5)):
		return 1
	else:
		return 0

#-------------------------------------------------------------------------------
# Check if correlation coefficients are all within -1 and 1
#-------------------------------------------------------------------------------
def isBounded (corr):
	if ((corr<=1.0).all() and (corr>=-1.0).all()):
		return 1
	else:
		return 0

#-------------------------------------------------------------------------------
# Calculates the average outgoing neutron energy for a given spectrum (PFNS)
#-------------------------------------------------------------------------------
def PFNS_aveEout (Eout, PFNS):
	return integrate.trapz(Eout*PFNS,Eout)

#-------------------------------------------------------------------------------
# Returns the Maxwellian function at the energy x and temperature T.
#-------------------------------------------------------------------------------
def Maxwellian (x,T):
	return 2.0/sqrt(3.14159*T**3)*sqrt(x)*exp(-x/T)

#-- vectorized version of the Maxwellian function
vMaxwellian = np.vectorize(Maxwellian)








