import os, sys

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
    if 'ma_1000_mA_1200'in sig_name: Xsec = 0.001418
    if 'ma_100_mA_1200'in sig_name: Xsec = 0.465000
    if 'ma_10_mA_1200'in sig_name: Xsec = 1.559000
    if 'ma_150_mA_1200'in sig_name: Xsec = 0.294700
    if 'ma_200_mA_1200'in sig_name: Xsec = 0.188900
    if 'ma_250_mA_1200'in sig_name: Xsec = 0.122900
    if 'ma_300_mA_1200'in sig_name: Xsec = 0.081890
    if 'ma_350_mA_1200'in sig_name: Xsec = 0.055240
    if 'ma_400_mA_1200'in sig_name: Xsec = 0.038170
    if 'ma_450_mA_1200'in sig_name: Xsec = 0.026740
    if 'ma_500_mA_1200'in sig_name: Xsec = 0.018900
    if 'ma_50_mA_1200'in sig_name: Xsec = 0.695000
    if 'ma_700_mA_1200'in sig_name: Xsec = 0.005690
    if 'ma_750_mA_1200'in sig_name: Xsec = 0.004344
    if 'ma_100_mA_600'in sig_name: Xsec = 0.468300
    if 'ma_10_mA_600'in sig_name: Xsec = 1.565000
    if 'ma_150_mA_600'in sig_name: Xsec = 0.299700
    if 'ma_200_mA_600'in sig_name: Xsec = 0.197100
    if 'ma_250_mA_600'in sig_name: Xsec = 0.131400
    if 'ma_300_mA_600'in sig_name: Xsec = 0.090020
    if 'ma_350_mA_600'in sig_name: Xsec = 0.063680
    if 'ma_400_mA_600'in sig_name: Xsec = 0.046770
    if 'ma_450_mA_600'in sig_name: Xsec = 0.035910
    if 'ma_500_mA_600'in sig_name: Xsec = 0.027360
    return Xsec