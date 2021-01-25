class cutFlow:
    def __init__(self,nPho,nEle_loose,nEle_tight,nMu_loose,nMu_tight,nTau,Werecoil,Wmurecoil,ZeeRecoil,ZmumuRecoil,pfMet,njets,njets_iso,nBjets,nBjets_iso,h_mass,nAK8JetsSR,nAK8JetsSBand,nAK8JetsZCR,ZeeMass,ZmumuMass):

        self.nPho         =  nPho
        self.nEle_loose   =  nEle_loose
        self.nEle_tight   =  nEle_tight
        self.nMu_loose    =  nMu_loose
        self.nMu_tight    =  nMu_tight
        self.nTau         =  nTau
        self.Werecoil     =  Werecoil
        self.Wmurecoil    =  Wmurecoil
	self.ZeeRecoil    =  ZeeRecoil
        self.ZmumuRecoil  =  ZmumuRecoil
        self.met          =  pfMet
        self.njets        =  njets
        self.njets_iso    =  njets_iso
        self.nBjets       =  nBjets
        self.nBjets_iso   =  nBjets_iso
        self.Mbb          =  h_mass
        self.nAK8JetsSR   =  nAK8JetsSR
        self.nAK8JetsSBand=  nAK8JetsSBand
        self.nAK8JetsZCR  =  nAK8JetsZCR
        self.ZeeMass      =  ZeeMass
        self.ZmumuMass    =  ZmumuMass
        self.value = 0
        self.cuts = {}


    def singleLepton_R(self,value,isEleRegion=True,isTop=True):

        if isEleRegion:
            bin1 = self.nEle_loose==1 and self.nEle_tight==1
            bin2 = bin1 and self.nMu_loose==0 and self.nTau==0 and self.nPho==0
            bin3 = bin2 and self.Werecoil > 200.0
        else:
            bin1 = self.nMu_loose==1 and self.nMu_tight==1
            bin2 = bin1 and self.nEle_loose==0 and self.nTau==0
            bin3 = bin2 and self.Wmurecoil > 200.0

        #bin3 = bin2 and self.recoil > 200.0
        #bin4 = bin3 and self.met > 100.0
        if isTop:
	    bin4 = bin3 and self.met > 50.0
	    bin5 = bin4 and self.njets >2
        else:
	    bin4 = bin3 and self.met > 100.0
	    bin5 = bin4 and self.njets ==2

        bin6 = bin5 and self.nBjets==2
        bin7 = bin6 and self.Mbb > 100.0 and self.Mbb < 150.0

        if bin7:     return value,value,value,value,value,value,value
        elif bin6:   return value,value,value,value,value,value,0.0
        elif bin5:   return value,value,value,value,value,0.0,0.0
        elif bin4:   return value,value,value,value,0.0,0.0,0.0
        elif bin3:   return value,value,value,0.0,0.0,0.0,0.0
        elif bin2:   return value,value,0.0,0.0,0.0,0.0,0.0
        elif bin1:   return value,0.0,0.0,0.0,0.0,0.0,0.0
        else:        return 0.0,0.0,0.0,0.0,0.0,0.0,0.0


    def singleLepton_B(self,value,isEleRegion=True,isTop=True):

        if isEleRegion:
            bin1 = self.nEle_loose==1 and self.nEle_tight==1
            bin2 = bin1 and self.nMu_loose==0 and self.nTau==0 and self.nPho==0
            bin3 = bin2 and self.Werecoil > 200.0
        else:
            bin1 = self.nMu_loose==1 and self.nMu_tight==1
            bin2 = bin1 and self.nEle_loose==0 and self.nTau==0
            bin3 = bin2 and self.Wmurecoil > 200.0

        #bin3 = bin2 and self.recoil > 200.0
        bin4 = bin3 and self.met > 50.0
        bin5 = bin4 and self.nAK8JetsSR==1
        if isTop:
            bin6 = bin5 and self.nBjets_iso==1
        else:
            bin6 = bin5 and self.nBjets_iso==0

        bin7 = bin6 #and self.njets_iso > 1

        if bin7:     return value,value,value,value,value,value,value
        elif bin6:   return value,value,value,value,value,value,0.0
        elif bin5:   return value,value,value,value,value,0.0,0.0
        elif bin4:   return value,value,value,value,0.0,0.0,0.0
        elif bin3:   return value,value,value,0.0,0.0,0.0,0.0
        elif bin2:   return value,value,0.0,0.0,0.0,0.0,0.0
        elif bin1:   return value,0.0,0.0,0.0,0.0,0.0,0.0
        else:        return 0.0,0.0,0.0,0.0,0.0,0.0,0.0


    def diLepton_R(self,value,isEleRegion=True):

        if isEleRegion:
	    bin1 = self.nEle_loose==2 and (self.nEle_tight==1 or self.nEle_tight==2)
            bin2 = bin1 and self.nMu_loose==0 and self.nTau==0 and self.nPho==0
            bin3 = bin2 and self.ZeeRecoil > 200.0
            ZMass = self.ZeeMass
        else:
            bin1 = self.nMu_loose==2 and (self.nMu_tight==1 or self.nMu_tight==2)
            bin2 = bin1 and self.nEle_loose==0 and self.nTau==0
            bin3 = bin2 and self.ZmumuRecoil > 200.0
            ZMass = self.ZmumuMass

        bin4 = bin3 and self.met < 100.0
        bin5 = bin4 and self.njets >2
        bin6 = bin5 and self.nBjets==2
        bin7 = bin6 and ZMass > 60.0 and ZMass < 120.0

        if bin7:     return value,value,value,value,value,value,value
        elif bin6:   return value,value,value,value,value,value,0.0
        elif bin5:   return value,value,value,value,value,0.0,0.0
        elif bin4:   return value,value,value,value,0.0,0.0,0.0
        elif bin3:   return value,value,value,0.0,0.0,0.0,0.0
        elif bin2:   return value,value,0.0,0.0,0.0,0.0,0.0
        elif bin1:   return value,0.0,0.0,0.0,0.0,0.0,0.0
        else:        return 0.0,0.0,0.0,0.0,0.0,0.0,0.0


    def diLepton_B(self,value,isEleRegion=True):

        if isEleRegion:
            bin1 = self.nEle_loose==2 and (self.nEle_tight==1 or self.nEle_tight==2)
            bin2 = bin1 and self.nMu_loose==0 and self.nTau==0 and self.nPho==0
            bin3 = bin2 and self.ZeeRecoil > 200.0
            ZMass = self.ZeeMass
        else:
            bin1 = self.nMu_loose==2 and (self.nMu_tight==1 or self.nMu_tight==2)
            bin2 = bin1 and self.nEle_loose==0 and self.nTau==0
            bin3 = bin2 and self.ZmumuRecoil > 200.0
            ZMass = self.ZmumuMass

        bin4 = bin3 and self.met > 0.0
        bin5 = bin4 and self.nAK8JetsZCR==1
        bin6 = bin5 and self.nBjets_iso==0
        bin7 = bin6 and ZMass > 60.0 and ZMass < 120.0

        if bin7:     return value,value,value,value,value,value,value
        elif bin6:   return value,value,value,value,value,value,0.0
        elif bin5:   return value,value,value,value,value,0.0,0.0
        elif bin4:   return value,value,value,value,0.0,0.0,0.0
        elif bin3:   return value,value,value,0.0,0.0,0.0,0.0
        elif bin2:   return value,value,0.0,0.0,0.0,0.0,0.0
        elif bin1:   return value,0.0,0.0,0.0,0.0,0.0,0.0
        else:        return 0.0,0.0,0.0,0.0,0.0,0.0,0.0

    def SBand_R(self,value):

        bin1 = self.nEle_loose==0 and self.nMu_loose==0
        bin2 = bin1 and self.nTau==0
        bin3 = bin2 and self.nPho==0
        bin4 = bin3 and self.met > 200.0
        bin5 = bin4 and self.njets>=0 
        bin6 = bin5 and self.nBjets==2
        bin7 = bin6 and ((self.Mbb > 50 and self.Mbb < 100) or (self.Mbb > 150 and self.Mbb < 350))

        if bin7:     return value,value,value,value,value,value,value
        elif bin6:   return value,value,value,value,value,value,0.0
        elif bin5:   return value,value,value,value,value,0.0,0.0
        elif bin4:   return value,value,value,value,0.0,0.0,0.0
        elif bin3:   return value,value,value,0.0,0.0,0.0,0.0
        elif bin2:   return value,value,0.0,0.0,0.0,0.0,0.0
        elif bin1:   return value,0.0,0.0,0.0,0.0,0.0,0.0
        else:        return 0.0,0.0,0.0,0.0,0.0,0.0,0.0


    def SBand_B(self,value):

        bin1 = self.nEle_loose==0 and self.nMu_loose==0
        bin2 = bin1 and self.nTau==0
        bin3 = bin2 and self.nPho==0
        bin4 = bin3 and self.met > 200.0
        bin5 = bin4 and self.nAK8JetsSBand==1
        bin6 = bin5 and self.nBjets_iso==0
        bin7 = bin6

        if bin7:     return value,value,value,value,value,value,value
        elif bin6:   return value,value,value,value,value,value,0.0
        elif bin5:   return value,value,value,value,value,0.0,0.0
        elif bin4:   return value,value,value,value,0.0,0.0,0.0
        elif bin3:   return value,value,value,0.0,0.0,0.0,0.0
        elif bin2:   return value,value,0.0,0.0,0.0,0.0,0.0
        elif bin1:   return value,0.0,0.0,0.0,0.0,0.0,0.0
        else:        return 0.0,0.0,0.0,0.0,0.0,0.0,0.0


