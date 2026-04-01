

setwd('~/Desktop/Research/ciona/spatialTranscriptomics/myData/')

# slide 1

slide1 <- read.csv('ciona_brain_nc1/spatial/tissue_positions_list.csv', header = FALSE)
slide_1_1 <- slide1
slide_1_2 <- slide1
slide_1_3 <- slide1

slide_1_1[slide_1_1$V6 > 4000,]$V2 <- 0
slide_1_2[(slide_1_2$V6 < 4000) | (slide_1_2$V6 > 7000),]$V2 <- 0
slide_1_3[slide_1_3$V6 < 7000,]$V2 <- 0


write.csv(slide_1_1,file = 'ciona_brain_nc1/spatial/tissue_positions_list_1.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_1_2,file = 'ciona_brain_nc1/spatial/tissue_positions_list_2.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_1_3,file = 'ciona_brain_nc1/spatial/tissue_positions_list_3.csv',row.names = FALSE,quote = FALSE)


## Command line 
xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_A_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/spatial/tissue_positions_list_1.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc1_1

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_A_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/spatial/tissue_positions_list_2.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc1_2

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_A_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/spatial/tissue_positions_list_3.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc1/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc1_3


# slide 2
slide2<- read.csv('ciona_brain_nc2/spatial/tissue_positions_list.csv', header = FALSE)
slide_2_1 <- slide2
slide_2_2 <- slide2
slide_2_3 <- slide2

slide_2_1[slide_2_1$V6 > 5000,]$V2 <- 0
slide_2_2[(slide_2_2$V6 < 5000) | (slide_2_2$V6 > 8000),]$V2 <- 0
slide_2_3[slide_2_3$V6 < 8000,]$V2 <- 0

write.csv(slide_2_1,file = 'ciona_brain_nc2/spatial/tissue_positions_list_1.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_2_2,file = 'ciona_brain_nc2/spatial/tissue_positions_list_2.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_2_3,file = 'ciona_brain_nc2/spatial/tissue_positions_list_3.csv',row.names = FALSE,quote = FALSE)

## Command line 
xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_B_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/spatial/tissue_positions_list_1.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc2_1

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_B_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/spatial/tissue_positions_list_2.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc2_2

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_B_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/spatial/tissue_positions_list_3.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc2/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc2_3

# slide 3
slide3<- read.csv('ciona_brain_nc3/spatial/tissue_positions_list.csv', header = FALSE)
slide_3_1 <- slide3
slide_3_2 <- slide3
slide_3_3 <- slide3

slide_3_1[slide_3_1$V6 > 4000,]$V2 <- 0
slide_3_2[(slide_3_2$V6 < 4000) | (slide_3_2$V6 > 7000),]$V2 <- 0
slide_3_3[slide_3_3$V6 < 7000,]$V2 <- 0

write.csv(slide_3_1,file = 'ciona_brain_nc3/spatial/tissue_positions_list_1.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_3_2,file = 'ciona_brain_nc3/spatial/tissue_positions_list_2.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_3_3,file = 'ciona_brain_nc3/spatial/tissue_positions_list_3.csv',row.names = FALSE,quote = FALSE)

## Command line 
xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_C2_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/spatial/tissue_positions_list_1.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc3_1

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_C2_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/spatial/tissue_positions_list_2.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc3_2

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_C2_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/spatial/tissue_positions_list_3.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc3/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc3_3


# slide 4
slide4<- read.csv('ciona_brain_nc4/spatial/tissue_positions_list.csv', header = FALSE)
slide_4_1 <- slide4
slide_4_2 <- slide4
slide_4_3 <- slide4

slide_4_1[slide_4_1$V6 > 4000,]$V2 <- 0
slide_4_2[(slide_4_2$V6 < 4000) | (slide_4_2$V6 > 7000),]$V2 <- 0
slide_4_3[slide_4_3$V6 < 7000,]$V2 <- 0

write.csv(slide_4_1,file = 'ciona_brain_nc4/spatial/tissue_positions_list_1.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_4_2,file = 'ciona_brain_nc4/spatial/tissue_positions_list_2.csv',row.names = FALSE,quote = FALSE)
write.csv(slide_4_3,file = 'ciona_brain_nc4/spatial/tissue_positions_list_3.csv',row.names = FALSE,quote = FALSE)

## Command line 
xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_D_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/spatial/tissue_positions_list_1.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc4_1

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_D_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/spatial/tissue_positions_list_2.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc4_2

xfuse convert visium \
--bc-matrix /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/raw_feature_bc_matrix.h5 \
--image /home/xzeng/project/Ciona/raw_data/spatial/IMAGE/220216_GE/220216_GE_HE_D_x10.tif \
--tissue-positions /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/spatial/tissue_positions_list_3.csv \
--scale-factors /home/xzeng/project/Ciona/raw_data/spatial/SPACERANGER/ciona_brain_nc4/outs/spatial/scalefactors_json.json \
--scale 0.2 \
--save-path nc4_3
