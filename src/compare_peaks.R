setwd('~/d/sci/src/fragments')
require(sqldf)
options(stringsAsFactors=FALSE)
peaks = read.table('data/processed/19f_peaks.tsv', sep='\t', header=FALSE) 

colnames(peaks) = c('filename','protein','nucleus','condition','row','column','shift','intensity')
peaks$peak_id = paste(peaks$condition, peaks$row, peaks$column, formatC(peaks$shift,format='f',digits=2,flag='0'), sep='_')
length(peaks$peak_id)
length(unique(peaks$peak_id))

protein_peaks = peaks[peaks$condition=='protein',c('peak_id', 'row', 'column', 'shift', 'intensity')]
control_peaks = peaks[peaks$condition=='control',c('peak_id', 'row', 'column', 'shift', 'intensity')]

head(protein_peaks)
head(control_peaks)

candidate_matches = sqldf("
select   c.peak_id cpeak, c.row, c.column, 
         c.shift cshift, c.intensity cintens, 
         p.peak_id ppeak, p.shift pshift, p.intensity pintens,
         abs(c.shift - p.shift) abs_delta, c.shift - p.shift delta
from     control_peaks c left outer join protein_peaks p
on       c.row = p.row
and      c.column = p.column
order by 2 asc, 3 asc, 9 asc
;")

best_matches = candidate_matches[!duplicated(candidate_matches$cpeak),]

best_matches$abs_delta_trunc = pmin(best_matches$abs_delta, .2)

hist(best_matches$delta, breaks=10000, border=NA, col='#FF9912',  xaxs='i', yaxs='i', xlim=c(-.1,.1),
     xlab='absolute difference in chemical shift (control vs. protein)',
     ylab='number of peaks',
     main='histogram of shift data')


hist(best_matches$abs_delta_trunc, breaks=100, border=NA, col='#FF9912', xaxt='n', xaxs='i', yaxs='i',
     xlab='absolute difference in chemical shift (control vs. protein)',
     ylab='number of peaks',
     main='histogram of shift data')
axis(side=1, at=(0:4)*5/100, labels=c('0', '.05', '.1', '.15', '.20+'))

sum(best_matches$abs_delta_trunc > .1)

plot(best_matches$cintens, best_matches$pintens, pch='.', col='#FF2016')

best_matches$intens_delta = best_matches$cintens - best_matches$pintens
plot(best_matches$abs_delta, abs(best_matches$intens_delta), pch=20, col='#FF9912', log='xy',
     xlab='Absolute change in chemical shift (log scale)', ylab='Absolute change in intensity (log scale)',
     main='Changes in peaks between control and protein conditions')

best_matches$intens_ratio = best_matches$pintens / best_matches$cintens
best_matches$delta_trunc = best_matches$delta
best_matches$delta_trunc[best_matches$delta > .2] = .2
best_matches$delta_trunc[best_matches$delta < -.2] = -.2
plot(best_matches$delta_trunc, best_matches$intens_ratio, pch=20, col='#FF9912', log='y',
     xlab='Change in chemical shift', ylab='Ratio of intensity',
     main='Changes in peaks between control and protein conditions')
abline(h=1,col='red')
abline(v=0,col='red')

# Ben Chittick's suggestion
best_matches$shift_med_z = (best_matches$delta - median(best_matches$delta))/mad(best_matches$delta)
sum(abs(best_matches$shift_med_z) > 5)

best_matches$intens_med_z = (log(best_matches$intens_ratio) - median(log(best_matches$intens_ratio)))/mad(log(best_matches$intens_ratio))
sum(abs(best_matches$intens_med_z) > 5)

best_matches$hit = abs(best_matches$shift_med_z) > 5 | abs(best_matches$intens_med_z) > 5
best_matches$col = '#FF9912'
best_matches$col[best_matches$hit] = '#0000CD'


plot(best_matches$delta_trunc, best_matches$intens_ratio, pch=20, col=best_matches$col, log='y',
     xlab='Change in chemical shift', ylab='Ratio of intensity',
     main='Changes in peaks between control and protein conditions')
abline(h=1,col='red')
abline(v=0,col='red')
sum(best_matches$hit)

# print list of hits
best_matches$cpeak[best_matches$hit]

table(paste(best_matches$row[best_matches$hit],best_matches$column[best_matches$hit],sep=''))




