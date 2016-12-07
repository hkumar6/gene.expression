source("common/env-kh.R")
source("common/generate-sample.R")
source("kknn/kknn-functions.R")

allGenes <- kknnImputeAllGenes()


# ggplot(allGenes, aes(x=reorder(optimalK, -table(optimalK)[optimalK]))) + geom_bar() + coord_flip()
