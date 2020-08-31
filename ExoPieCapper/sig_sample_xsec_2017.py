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