stop("bet you didn't mean to run the whole script")
setwd('~/d/sci/src/fragments')
library(sqldf)
library(dplyr)
library(tidyr)
library(magrittr)
options(stringsAsFactors=FALSE)
source('../miscellanea/r_helper.r')
peaks = read.table('data/processed/19f_peaks.tsv', sep='\t', header=FALSE) 

colnames(peaks) = c('filename','protein','nucleus','condition','row','column','shift','intensity')
peaks$peak_id = paste(peaks$condition, peaks$row, peaks$column, formatC(peaks$shift,format='f',digits=2,flag='0'), sep='_')
length(peaks$peak_id)
length(unique(peaks$peak_id))

protein_peaks = peaks[peaks$condition=='protein',c('peak_id', 'row', 'column', 'shift', 'intensity')]
control_peaks = peaks[peaks$condition=='control',c('peak_id', 'row', 'column', 'shift', 'intensity')]

reference_peaks = read.table('data/f19_library_072815.tsv', sep='\t', header=TRUE, blank.lines.skip=TRUE, quote='', comment.char='', allowEscapes=FALSE)
colnames(reference_peaks) = fix_colnames(colnames(reference_peaks))

reference_peaks$well_id = paste(reference_peaks$pool_plate, reference_peaks$well_number)

reference_peaks = reference_peaks[!is.na(reference_peaks$pool__),]

reference_peaks$pooled_col = as.integer(gsub('[A-Za-z ]*','',reference_peaks$pool_plate))
reference_peaks$pooled_row = gsub('[0-9]*','',reference_peaks$well_number)
reference_peaks$tfa_shift_1 = as.numeric(reference_peaks$chemical_shift_1__tfa_)
reference_peaks$tfa_shift_2 = as.numeric(reference_peaks$chemical_shift_2__tfa_)
reference_peaks$tfa_shift_3 = as.numeric(reference_peaks$chemical_shift_3__tfa_)
reference_peaks$tfa_shift_4 = as.numeric(reference_peaks$chemical_shift_4__tfa_)
reference_peaks$tfa_shift_5 = as.numeric(reference_peaks$chemical_shift_5__tfa_)
reference_peaks$tfa_shift_6 = as.numeric(reference_peaks$chemical_shift_6__tfa_)
reference_peaks$tfa_shift_7 = as.numeric(reference_peaks$chemical_shift_7__tfa_)
reference_peaks$tfa_shift_8 = as.numeric(reference_peaks$chemical_shift_8__tfa_)
reference_peaks$tfa_shift_9 = as.numeric(reference_peaks$chemical_shift_9__tfa_)
reference_peaks$tfa_shift_10 = as.numeric(reference_peaks$chemical_shift_10__tfa_)
reference_peaks$tfa_shift_11 = as.numeric(reference_peaks$chemical_shift_11__tfa_)
reference_peaks$tfa_shift_12 = as.numeric(reference_peaks$chemical_shift_12__tfa_)
reference_peaks = reference_peaks[,-(18:49)]
colnames(reference_peaks)

# tidy up the multiple shift columns
ref = reference_peaks %>% gather(shift_number, shift, tfa_shift_1:tfa_shift_12)
ref = ref[!is.na(ref$shift),]

# check that the peaks line up reasonably well
A1ref = ref$pooled_col==1 & ref$pooled_row=='A'
plot(ref$shift[A1ref],rep(1,sum(A1ref)),type='h',col='#0001CD',lwd=5,ann=FALSE,axes=FALSE)
A1con = control_peaks$column==1 & control_peaks$row=='A'
points(control_peaks$shift[A1con], rep(1,sum(A1con)),type='h',col='#FF9912',lwd=5)
axis(side=1, at=(-20:-1)*10)
legend('topright',c("Alison's reference","my controls"),col=c('#0001CD','#FF9912'),lwd=5,bty='n')
title('pool A1 as an example')

# try to match peaks
hits = read.table('data/19f_hits_manual_review.tsv',header=TRUE,sep='\t',quote='',comment.char='')
hits = hits[hits$valid=='yes',]
hit_table = unique(hits[,c('row','column')])
table(paste(hits$row, hits$column, sep=''))


for (k in 1:nrow(hit_table)) {
  row = hit_table$row[k]
  col = hit_table$col[k]
# to test just one, uncomment these::
# row = 'G'
# col = 7
ref_select = reference_peaks$pooled_row==row & reference_peaks$pooled_col==col
con_select = control_peaks$row==row & control_peaks$col==col
hit_select = hits$row == row & hits$column == col

png(paste('data/pool',row,col,'.png',sep=''),width=1600,height=600,res=150)
xlim = range(c(reference_peaks$tfa_shift_1[ref_select & !is.na(ref_select)], reference_peaks$tfa_shift_2[ref_select & !is.na(ref_select)], control_peaks$shift[con_select]), na.rm=TRUE)
ylim = c(0, max(control_peaks$intensity))
plot(NA, NA, xlim=xlim, ylim=ylim, axes=FALSE, ann=FALSE, yaxs='i')
peak_name = paste(reference_peaks$pool_plate[ref_select], reference_peaks$well_number[ref_select])
for (i in 1:12) {
  abline(v=reference_peaks[ref_select,paste('tfa_shift_',i,sep='')]+1, col='red', lwd=2)
  text(x=reference_peaks[ref_select,paste('tfa_shift_',i,sep='')]+1, y=10, labels=peak_name, pos=2, srt=90, col='red', cex=.7)
}
points(control_peaks$shift[con_select], control_peaks$intensity[con_select], type='h', col='blue', lwd=2)
abline(h=0)
points(hits$cshift[hit_select], hits$cintens[hit_select], pch=17)
axis(side=1, at=(-18:-1)*10, labels=(-18:-1)*10, lwd=0, lwd.ticks=1)
title(paste('pool ',row,col,sep=''))
dev.off()
}
## note how the hit in F1 doesn't line up quite as well
## and hit in G7 seems to have no match
## G10 is also a little weird - the one next to it is orphaned
## H9 is totally orphaned
# evidence of promiscuity?


#match all our control peaks to Alison's peaks

cand_ref_matches = sqldf("
select   r.well_id, r.pooled_row, r.pooled_col, r.source, r.smiles, r.shift rshift,
         c.peak_id cpeak, c.shift cshift, c.intensity cintens, 
         abs(c.shift - r.shift) abs_delta, c.shift - r.shift delta
from     ref r left outer join control_peaks c
on       r.pooled_row = c.row
and      r.pooled_col = c.column
order by 1 asc, 10 asc
;")

best_ref_matches = cand_ref_matches[!duplicated(cand_ref_matches$well_id),]

hist(best_ref_matches$abs_delta, breaks=1000, xlim=c(0,2), col='#555555', border=NA)

hist(best_ref_matches$delta, breaks=1000, xlim=c(0,2), col='#555555', border=NA)



cand_ref_matches = sqldf("
select   r.well_id, r.pooled_row, r.pooled_col, r.source, r.name, r.smiles, r.mw, r.shift rshift,
         c.peak_id cpeak, c.shift cshift, c.intensity cintens, 
         abs(c.shift - r.shift - 1) abs_delta_1
from     ref r left outer join control_peaks c
on       r.pooled_row = c.row
and      r.pooled_col = c.column
order by 1 asc, 12 asc
;")
best_ref_matches = cand_ref_matches[!duplicated(cand_ref_matches$well_id),]
hist(best_ref_matches$abs_delta, breaks=1000, xlim=c(0,2), col='#555555', border=NA)

hit_matches = sqldf("
select   b.well_id, b.pooled_row, b.pooled_col, b.source, b.name, b.smiles, b.mw, b.rshift, b.cshift, b.abs_delta_1 match_delta, h.delta_trunc hit_delta, h.valid, h.comments, h.how_big
from     hits h, best_ref_matches b
where    h.cpeak = b.cpeak
;")

write.table(hit_matches, 'data/19f_hits_matched.tsv', row.names=F, quote=F, col.names=T, sep='\t')

plot(hit_matches$match_delta, abs(hit_matches$hit_delta), pch=20, xlab='Alison vs. control', ylab='Control vs. protein')

cat(paste(hit_matches$smiles[hit_matches$match_delta < .3], collapse='.'))

