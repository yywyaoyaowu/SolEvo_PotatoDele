# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 8/2/21
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/03_Kmer.tree")
library(ape)

tree.kmer <- read.tree("topology.nwk")
tree.kmer.root=root(tree.kmer,11:12,resolve.root =TRUE)
plot(root(tree.kmer,11:12,resolve.root =TRUE),size=0.3)
write.tree(tree.kmer.root,"topology.root_11.12.nwk")


tree.kmer.root=root(tree.kmer,52:58,resolve.root =TRUE)
plot(root(tree.kmer,52:58,resolve.root =TRUE),size=0.3)
write.tree(tree.kmer.root,"topology.root_52:58.nwk")

tree.kmer_root.read <- read.tree("topology.root.nwk")
str(tree.kmer)
tree.kmer
root(tree.kmer, outgroup, node = NULL, resolve.root = FALSE,
     interactive = FALSE, edgelabel = FALSE, ...)

data(bird.orders)
tr <- root(bird.orders, 1)
plot(root(tr,1,resolve.root =TRUE))

