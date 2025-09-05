#Software/program import
module load CCEnv StdEnv/2020 gcc/9.3.0 iq-tree/2.2.2.7

#IQTree model finder and run
iqtree2 -s iqtree/continent_alignment.fasta -m MF
iqtree2 -s iqtree/continent_alignment.fasta -m TPM2u+F+I+R7 -st DNA -nt AUTO --prefix iqtree/continent_tree -B 1000 --redo

#ClonalFrameML tree reconstruction
$HOME/bin/ClonalFrameML/src/ClonalFrameML clonalframeML/cftree.labelled_tree.newick clonalframeML/final_alignment.mfa clonalframeML/cftreebranch -embranch true -embranch_dispersion 0.1 -initial_values "0.12796937 16.035942 0.086674745"
