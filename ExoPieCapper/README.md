**For MET Channel:**

```python StackPlotter_syst.py -d MET -m -y 2016 -D analysis_histo_v16_06-04-05 -s -S analysis_histo_v16_06-04-05```


**For Electron Channel:**

```python StackPlotter_syst.py -d SE -e -y 2016 -D analysis_histo_v16_06-04-05```


for change in the directory please update in the following line:
https://github.com/tiwariPC/ExoPieProducer/blob/bbDM_withSyst/ExoPieCapper/StackPlotter_syst.py#L93
ANd then use ```-y``` command for same year you have given the directory for.
