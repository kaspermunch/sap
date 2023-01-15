
from Bio import Phylo
from Bio.Phylo import draw_ascii
tree = Phylo.read("tree.nw", "newick")
draw_ascii(tree)
