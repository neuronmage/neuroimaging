atlas.list <- as.factor(c('36','60','125','250','P1_16','P17_28', 'p29_36'))
atlas.count = 0 
for (atlas in levels(atlas.list)) {

hemi.list<- as.factor(c('L','R'))
hemi.count= 0 
for (hemi in levels(hemi.list)) { 

print(atlas)
print(hemi)

original_color_atlas_hemi <- read.delim(paste('~/Desktop/colortable_and_gcs/original_color_',atlas,'_',hemi,'.txt', sep=""), header= FALSE, sep = "\t")
numregions<-nrow(original_color_atlas_hemi)
print(numregions)
rgbcol1<- as.integer(runif(numregions, 20, 225));
rgbcol2<-as.integer(runif(numregions, 20, 225));
rgbcol3<-as.integer(runif(numregions, 20, 225));
new_LUT<- cbind(original_color_atlas_hemi[c("V1", "V2")], rgbcol1, rgbcol2, rgbcol3, original_color_atlas_hemi[c("V6")])
write.table(new_LUT, file=paste('~/Desktop/colortable_and_gcs/new/original_color_',atlas,'_',hemi,'.txt', sep= ""), sep = " ", col.names= FALSE, row.names= FALSE, quote = FALSE)
hemi.count= 0+1 
}
}
