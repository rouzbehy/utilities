# format in pdg.dat
#	0	pdg_id
#	1	name
#	2	mass
#	3	width
#	4	2*Spin+1
#	5	Baryon Number
#	6	Strange Number
#	7	Charm Number
#	8	Bottom Number
#	9	2*IsoSpin+1
#	10	Electric Charge
#	11	number of decay modes

PDG_Name= {}
PDG_Mass= {}
PDG_Width= {}
PDG_Spin= {}
PDG_IsoSpin= {}
PDG_Baryon= {}
PDG_Strange= {}
PDG_Charm= {}
PDG_Bottom= {}
PDG_Charge= {}
PDG_Pool= []

for aLine in open('pdg.dat'):
	lineinfo = aLine.split()
	if len(lineinfo) ==12 :
		pdgid = int(lineinfo[0])
		Name = lineinfo[1]
		Mass = float(lineinfo[2])
		Width = float(lineinfo[3])
		Spin = int(lineinfo[4])
		IsoSpin = int(lineinfo[9])
		Baryon = int(lineinfo[5])
		Strange = int(lineinfo[6])
		Charm = int(lineinfo[7])
		Bottom = int(lineinfo[8])
		Charge = int(lineinfo[10])
		PDG_Name[pdgid] = Name
		PDG_Mass[pdgid] = Mass
		PDG_Width[pdgid] = Width
		PDG_Spin[pdgid] = Spin
		PDG_IsoSpin[pdgid] = IsoSpin
		PDG_Baryon[pdgid] = Baryon
		PDG_Strange[pdgid] = Strange
		PDG_Charm[pdgid] = Charm
		PDG_Bottom[pdgid] = Bottom
		PDG_Charge[pdgid] = Charge
		PDG_Pool.append(pdgid)
		if Baryon > 0 :
			PDG_Name[-pdgid] = 'anti-%s'%Name
			PDG_Mass[-pdgid] = Mass
			PDG_Width[-pdgid] = Width
			PDG_Spin[-pdgid] = Spin
			PDG_IsoSpin[-pdgid] = IsoSpin
			PDG_Baryon[-pdgid] = -Baryon
			PDG_Strange[-pdgid] = -Strange
			PDG_Charm[-pdgid] = -Charm
			PDG_Bottom[-pdgid] = -Bottom
			PDG_Charge[-pdgid] = -Charge
			PDG_Pool.append(-pdgid)

