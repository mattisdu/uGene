@ECHO OFF
ECHO Start uGene core.
python uGeneCore.py -f Example/example.phyloprofile -t [{'dev_report':0,'y_axis':'geneID','x_axis':'ncbiID','values':['FAS_F','FAS_B'],'jobs':'gene'}]
ECHO Done !
PAUSE