import os
import sys
sys.path.append('../../ExoPieUtils/commonutils/')

def ZRecoil_Phi_Zmass(nEle, eleCharge_, elepx_, elepy_, elepz_, elee_,met_,metphi_):
    dummy=-9999.0
    ZeeRecoilPt=dummy; ZeerecoilPhi=-10.0; Zee_mass=dummy
    if nEle == 2:
        iele1=0
        iele2=1
        if eleCharge_[iele1]*eleCharge_[iele2]<0:
            Zee_mass = InvMass(elepx_[iele1],elepy_[iele1],elepz_[iele1],elee_[iele1],elepx_[iele2],elepy_[iele2],elepz_[iele2],elee_[iele2])
            zeeRecoilPx = -( met_*math.cos(metphi_) + elepx_[iele1] + elepx_[iele2])
            zeeRecoilPy = -( met_*math.sin(metphi_) + elepy_[iele1] + elepy_[iele2])
            ZeeRecoilPt =  math.sqrt(zeeRecoilPx**2  +  zeeRecoilPy**2)
            ZeerecoilPhi = mathutil.ep_arctan(zeeRecoilPx,zeeRecoilPy)
    return ZeeRecoilPt, ZeerecoilPhi, Zee_mass

def WRecoil_Phi_Wmass(nEle,elept,elephi,elepx_,elepy_,met_,metphi_):
    dummy=-9999.0
    WenuRecoilPt=dummy; WenurecoilPhi=-10.0;  We_mass=dummy;
    if nEle == 1:
        We_mass = MT(elept[0],met_, DeltaPhi(elephi[0],metphi_)) #transverse mass defined as sqrt{2pT*MET*(1-cos(dphi)}
        WenuRecoilPx = -( met_*math.cos(metphi_) + elepx_[0])
        WenuRecoilPy = -( met_*math.sin(metphi_) + elepy_[0])
        WenuRecoilPt = math.sqrt(WenuRecoilPx**2  +  WenuRecoilPy**2)
        WenurecoilPhi = mathutil.ep_arctan(WenuRecoilPx,WenuRecoilPy)
    return WenuRecoilPt, WenurecoilPhi, We_mass
