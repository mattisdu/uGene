# uGene

### Project snapshot:
With uGene, phylogenetic profiles can be read and clustered with uMap. 
Groups of genes or taxa of particular interest can be exported.

### Installation:
##### Direct Download:
- Download from: https://github.com/Mattisudu/uGene/
- Extract to your desired location.

##### Via Pip:
1. Open 'packages.txt'.
2. Run: pip install -r requirements.txt

##### Via Conda/Pip:
1. Install Anaconda.
2. Execute: 
   conda create -n uGene
   conda activate uGene
   conda install pip
   pip install -r requirements.txt

##### Enable PhyloProfile Support:
1. Install PhyloProfile. Visit: https://github.com/BIONF/PhyloProfile for install instructions.
2. Use PhyloProfile to open and parse your taxa. Ensuring all taxa are recognized.
3. Ensure R is installed and on your system's PATH.

### Start uGene:
You can use uGene via command-line or the dashboard.
#### Command-line use:
   python uGeneCore.py -h
#### Garfic-Interface use:
   python uGeneGUI.py
