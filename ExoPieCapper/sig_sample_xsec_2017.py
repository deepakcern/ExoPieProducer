import os,sys

def getSigXsec(sig_name):
    if 'Ma1000_MChi1_MA1200' in sig_name: Xsec = 0.004647
    if 'Ma100_MChi1_MA600' in sig_name: Xsec = 14.850000
    if 'Ma10_MChi1_MA1200' in sig_name: Xsec = 431.400000
    if 'Ma10_MChi1_MA600' in sig_name: Xsec = 433.400000
    if 'Ma150_MChi1_MA1200' in sig_name: Xsec = 5.738000
    if 'Ma150_MChi1_MA600' in sig_name: Xsec = 5.746000
    if 'Ma200_MChi1_MA1200' in sig_name: Xsec = 2.595000
    if 'Ma200_MChi1_MA600' in sig_name: Xsec = 2.639000
    if 'Ma250_MChi1_MA1200' in sig_name: Xsec = 1.314000
    if 'Ma250_MChi1_MA600' in sig_name: Xsec = 1.345000
    if 'Ma300_MChi1_MA600' in sig_name: Xsec = 0.747900
    if 'Ma350_MChi1_MA1200' in sig_name: Xsec = 0.415900
    if 'Ma350_MChi1_MA600' in sig_name: Xsec = 0.448900
    if 'Ma400_MChi1_MA1200' in sig_name: Xsec = 0.251900
    if 'Ma400_MChi1_MA600' in sig_name: Xsec = 0.289900
    if 'Ma450_MChi1_MA1200' in sig_name: Xsec = 0.159000
    if 'Ma450_MChi1_MA600' in sig_name: Xsec = 0.198300
    if 'Ma500_MChi1_MA1200' in sig_name: Xsec = 0.103100
    if 'Ma500_MChi1_MA600' in sig_name: Xsec = 0.139400
    if 'Ma50_MChi1_MA1200' in sig_name: Xsec = 55.390000
    if 'Ma50_MChi1_MA600' in sig_name: Xsec = 55.710000
    if 'Ma700_MChi1_MA1200' in sig_name: Xsec = 0.023600
    if 'Ma750_MChi1_MA1200' in sig_name: Xsec = 0.017120
    return Xsec

def getSigXsec_official(sig_name):
    if 'ma_400_mA_1200' in sig_name: Xsec =  0.03766
    if 'ma_10_mA_1200' in sig_name: Xsec =  0.9053
    if 'ma_100_mA_1200' in sig_name: Xsec =  0.4096
    if 'ma_1000_mA_1200' in sig_name: Xsec =  0.001415
    if 'ma_150_mA_1200' in sig_name: Xsec =  0.2694
    if 'ma_200_mA_1200' in sig_name: Xsec =  0.1733
    if 'ma_250_mA_1200' in sig_name: Xsec =  0.11630
    if 'ma_300_mA_1200' in sig_name: Xsec =  0.076770
    if 'ma_350_mA_1200' in sig_name: Xsec =  0.052710
    if 'ma_400_mA_1200' in sig_name: Xsec =  0.036630
    if 'ma_450_mA_1200' in sig_name: Xsec =  2.66E-02
    if 'ma_50_mA_1200' in sig_name: Xsec =  0.62480
    if 'ma_500_mA_1200' in sig_name: Xsec =  0.018970
    if 'ma_700_mA_1200' in sig_name: Xsec =  0.005750
    if 'ma_750_mA_1200' in sig_name: Xsec =  0.0043760
    if 'ma_10_mA_600' in sig_name: Xsec =  0.91120
    if 'ma_100_mA_600' in sig_name: Xsec =  0.4232
    if 'ma_150_mA_600' in sig_name: Xsec =  0.282800
    if 'ma_200_mA_600' in sig_name: Xsec =  0.1880
    if 'ma_250_mA_600' in sig_name: Xsec =  0.12550
    if 'ma_300_mA_600' in sig_name: Xsec =  0.085840
    if 'ma_350_mA_600' in sig_name: Xsec =  0.061830
    if 'ma_400_mA_600' in sig_name: Xsec =  0.045220
    if 'ma_450_mA_600' in sig_name: Xsec =  0.034210
    if 'ma_50_mA_600' in sig_name: Xsec =  0.62520
    if 'ma_500_mA_600' in sig_name: Xsec =  0.0257900
    return Xsec
