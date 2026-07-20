library(ape)
library(ggtree)
library(ggtreeExtra)

tree <- read.tree("coyote_dogs_palaeofaeces_aligned.phy_phyml_tree.txt")

map <- c(
  "KM061574.1"="KM061574_A1a1a",
  "KM061588.1"="KM061588_A1a1b",
  "EU789740.1"="EU789740_A1a1c",
  "DQ480496.1"="DQ480496_A1a1d",
  "EU789756.1"="EU789756_A1a2a",
  "KM061501.1"="KM061501_A1a2b",
  "EU789721.1"="EU789721_A1b1a1",
  "EU408280.2"="EU408280_A1b1a2",
  "EU789759.1"="EU789759_A1b1b",
  "KM061506.1"="KM061506_A1b2a1",
  "KM061529.1"="KM061529_A1b2a2",
  "DQ480497.1"="DQ480497_A1b2b",
  "KU290442.1"="KU290442_A1b3",
  "KM061567.1"="KM061567_A1b4",
  "EU789686.1"="EU789686_A1b5",
  "EU789717.1"="EU789717_A1b6",
  "KU290608.1"="KU290608_A1b7",
  "KM061542.1"="KM061542_A2a1a",
  "EU789677.1"="EU789677_A2a1b",
  "KU290617.1"="KU290617_A2b1a",
  "EU789673.1"="EU789673_A2b2",
  "KT168369.1"="KT168369_A2b2a",
  "EU789776.1"="EU789776_A2b3a",
  "EU789681.1"="EU789681_A2b3b",
  "KM113774.1"="KM113774_A3a1",
  "EU789692.1"="EU789692_A3b1",
  "EU789669.1"="EU789669_A4a",
  "EU789666.1"="EU789666_A5a",
  "EU408300.1"="EU408300_A6a",
  "EU789651.1"="EU789651_B2a1",
  "EU789757.1"="EU789757_B2b",
  "EU408268.1"="EU408268_B1a2",
  "EU789646.1"="EU789646_B1a1",
  "KU290510.1"="KU290510_C1a",
  "KJ637137.1"="KJ637137_C1b1",
  "KU290450.1"="KU290450_C1b2",
  "KX379529.1"="KX379529_C2",
  "EU789655.1"="EU789655_D1a",
  "EU789773.1"="EU789773_D2a",
  "EU789662.1"="EU789662_E",
  "AB499817.1"="AB499817_F",
  "NC_008093.1"="NC_008093_coyote",
  "KT168372"="KT168372_kuri",
  "KT168373"="KT168373_kuri",
  "KT168374"="KT168374_kuri",
  "KU215693"="KU215693_kuri",
  "KU215698"="KU215698_kuri",
  "KU215702"="KU215702_kuri",
  "KY798508"="KY798508_kuri",
  "KY798509"="KY798509_kuri",
  "KY798510"="KY798510_kuri",
  "KY798516"="KY798516_kuri",
  "MS11669_dr"="MS11669_palaeofaeces",
  "MS11670_dr"="MS11670_palaeofaeces",
  "MS11673_dr"="MS11673_palaeofaeces",
  "MS11674_dr"="MS11674_palaeofaeces",
  "MS11675_dr"="MS11675_palaeofaeces",
  "MS11676_dr"="MS11676_palaeofaeces",
  "MS11677_dr"="MS11677_palaeofaeces",
  "MS11678_dr"="MS11678_palaeofaeces",
  "MS11679_dr"="MS11679_palaeofaeces",
  "MS11683_dr"="MS11683_palaeofaeces",
  "MS11684_dr"="MS11684_palaeofaeces",
  "MS11686_dr"="MS11686_palaeofaeces",
  "MS11770_dr"="MS11770_palaeofaeces",
  "MS11771_dr"="MS11771_palaeofaeces",
  "MS11774_dr"="MS11774_palaeofaeces",
  "MS11775_dr"="MS11775_palaeofaeces"
)

tree$tip.label <- ifelse(
  tree$tip.label %in% names(map),
  map[tree$tip.label],
  tree$tip.label
)

tree$node.label <- sprintf("%.2f", as.numeric(tree$node.label))

rooted <- root(tree, outgroup="NC_008093_coyote", resolve.root=TRUE)

#Shortend branch length of Coyote for better visualisation
rooted$edge.length[rooted$edge[,2] == which(rooted$tip.label == "NC_008093_coyote")] <- 
  rooted$edge.length[rooted$edge[,2] == which(rooted$tip.label == "NC_008093_coyote")] * 0.15

ggtree(rooted) +
  geom_tiplab(size=5) +
  geom_nodepoint(size=2) +
  geom_nodelab(aes(label = ifelse(as.numeric(label) > 0.8, label, "")),
               size=3, vjust=-0.7, hjust=1.3, colour="red3")

ggsave("mitogenome_tree_pub.png", width = 12, height = 15, bg = "white")

