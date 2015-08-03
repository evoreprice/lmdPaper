
# real results
tpm <- readRDS('finalAnalysis/tpm/tpm.Rds')
gtfLength <- readRDS('finalAnalysis/tpm/gtfLength.Rds')

# intergenic shuffle
tpmInt <- readRDS('finalAnalysis/dnaTpm/tpm.Rds')
gtfLengthInt <- readRDS('finalAnalysis/dnaTpm/gtfLength.Rds')

# compare length distributions
g <- data.table(gtfLength, keep.rownames = TRUE)
setkey(g, 'rn')
gI <- data.table(gtfLengthInt, keep.rownames = TRUE)
setkey(gI, 'rn')
plotData <- merge(g, gI, all = TRUE)
setnames(plotData, c('gene_id', 'GTF length', 'Shuffled length'))
pd.melt <- reshape2::melt(plotData, id.vars = 'gene_id', measure.vars = c('GTF length', 'Shuffled length'))
ggplot(pd.melt, aes(x = log(value), colour = variable, fill = variable), alpha = 0.5) +
  theme_minimal() +
  scale_fill_brewer(palette = 'Set1') +
  scale_colour_brewer(palette = 'Set1') +
  geom_line(stat = 'density') +
  geom_density(colour = NA, alpha = 0.5)
# distributions totally different

# 95th percentile of the tpm calculated from intergenic regions
q <- quantile(tpmInt, 0.95, type = 7)

sapply(c(1:9), function(i)
  quantile(tpmInt, 0.9, type = i))

# find intersection between densities
upper.limit <- log(50)
lower.limit <- 0
tpm.density <- density(log(tpm[,'N1R3'] + 0.5), from = lower.limit, to = upper.limit, n = 2^10)
tpmInt.density <- density(log(tpmInt[,'N1R3'] + 0.5), from = lower.limit, to = upper.limit, n = 2^10)

density.difference <- tpm.density$y - tpmInt.density$y
intersection.points <- tpm.density$x[which(diff(density.difference > 0) != 0) + 1]

exp(intersection.points)


# test view double-distribution
pd <- data.table(tpm, keep.rownames = TRUE)
pd[,colnames(pd)[!colnames(pd) %in% c('rn', 'N1R1')] := NULL]
setnames(pd, "N1R1", "Genic TPM")

pd2 <- data.table(tpmInt, keep.rownames = TRUE)
pd2[,colnames(pd2)[!colnames(pd2) %in% c('rn', 'N1R1')] := NULL]
setnames(pd2, "N1R1", "Intergenic TPM")

setkey(pd, 'rn')
setkey(pd2, 'rn')

pd.melt <- reshape2::melt(merge(pd, pd2, all = TRUE), id.vars = 'rn')

cols <- RColorBrewer::brewer.pal(9, 'Set1')

# log
ggplot(pd.melt, aes(x = log(value + 0.5), fill = variable, colour = variable)) +
  theme_minimal() +
  #  xlim(log(0.25), log(1000)) +
  scale_fill_brewer(palette = 'Set1') +
  scale_colour_brewer(palette = 'Set1') +
  geom_line(stat = 'density') +
  geom_density(colour = NA, alpha = 0.2) +
  #geom_vline(x = log(q + 0.5), colour = cols[3]) +
  #geom_vline(x = log(8 + 0.5), colour = cols[4]) +
  geom_vline(x = intersection.points, colour = cols[5])

# not log
ggplot(pd.melt, aes(x = value, fill = variable, colour = variable)) +
  theme_minimal() +
  xlim(-1, 10) +
  scale_fill_brewer(palette = 'Set1') +
  scale_colour_brewer(palette = 'Set1') +
  geom_line(stat = 'density') +
  geom_density(colour = NA, alpha = 0.2) +
  geom_vline(x = q, colour = cols[3]) +
  geom_vline(x = 4.4, colour = cols[4])

# number of expressed genes (above q)
apply(tpm, 2, function(x)
  table(x > q))

apply(tpm, 2, function(x)
  table(x > exp(intersection.points[3])))

apply(tpmInt, 2, function(x)
  table(x > q))

# high "intergenic" "genes"
head(
  tpmInt[rev(order(tpmInt[,'N1R1'])),]
)


ggplot(data.frame(log2(tpm + 0.5)), aes(x = N1R1)) +
  geom_line(stat = 'density') +
  xlim(0, log2(100)) +
  geom_vline(x = log2(q + 00.5))

ggplot(data.frame(log2(tpmInt + 0.5)), aes(x = N1R1)) +
  xlim(log2(0.25), log2(50)) +
  geom_density()

hist(tpmInt[,'N1R1'])
which.max(tpmInt[,'N1R1'])

counts(dds)['LOC_Os08g44560',]


# monster tpm table for all comparisons
# first test it

tpm.test <- data.frame(tpm)
tpm.test$ID <- rownames(tpm.test)
tpmInt.test <- data.frame(tpmInt)
tpmInt.test$ID <- rownames(tpmInt.test)

tpm.test.melt <- data.table(reshape2::melt(tpm.test, id.vars = 'ID'))
setkey(tpm.test.melt, 'ID', 'variable')

tpmInt.test.melt <- data.table(reshape2::melt(tpmInt.test, id.vars = 'ID'))
setkey(tpmInt.test.melt, 'ID', 'variable')

# only compare genes in shuffled gtf, otherwise use all = TRUE
tpm.merged <- merge(tpm.test.melt, tpmInt.test.melt)
setnames(tpm.merged, old = c('variable', 'value.x', 'value.y'), new = c('Sample', 'Genic TPM', 'Intergenic TPM'))

tpm.merged.melt <- reshape2::melt(tpm.merged, id.vars = c('ID', 'Sample'), value.name = 'TPM')

# calculate the cutoff for whole datasets at once
upper.limit <- log(50)
lower.limit <- 0
tpm.density <- density(log(tpm.merged$`Genic TPM` + 0.5), from = lower.limit, to = upper.limit, n = 2^10)
tpmInt.density <- density(log(tpm.merged$`Intergenic TPM` + 0.5),
                          from = lower.limit, to = upper.limit, n = 2^10, na.rm = TRUE)

density.difference <- tpm.density$y - tpmInt.density$y
intersection.points <- tpm.density$x[which(diff(density.difference > 0) != 0) + 1]

cutoff <- exp(intersection.points) - 0.5
q <- quantile(tpm.merged$`Intergenic TPM`, 0.95, na.rm = TRUE)

cols <- RColorBrewer::brewer.pal(9, 'Set1')

ggplot(tpm.merged.melt, aes(x = log(TPM + 0.5), colour = variable, fill = variable)) +
  theme_minimal() +
  #facet_wrap(~Sample, ncol = 3, scales = 'free_y') +
  xlab(expression(log[italic(e)](TPM+0.5))) +
  ylab('Density') +
  xlim(log(0.25), 4) +
  scale_fill_brewer(palette = 'Set1', guide = FALSE) +
  scale_colour_brewer(palette = 'Set1', guide = guide_legend(title = NULL)) +
  geom_line(stat = 'density') +
  geom_density(colour = NA, alpha = 0.2) +
  geom_vline(x = intersection.points, colour = cols[3]) +
  annotate(x = intersection.points, y = 1, geom = 'text', hjust = -0.1, colour = cols[3],
           label = paste0('TPM = ', round(cutoff, 1)))
#geom_vline(x = log(q + 0.5), colour = cols[4]) +
#annotate(x = log(q + 0.5), y = 1, geom = 'text', hjust = -0.1, colour = cols[4],
#           label = paste0('95th percentile\nTPM = ',
#                         round(q), 1))

saveRDS(cutoff[length(cutoff)], 'finalAnalysis/dnaTpm/cutoff.Rds')

apply(tpm, 2, function(x)
  table(x > cutoff[length(cutoff)]))

apply(tpmInt, 2, function(x)
  table(x > cutoff[length(cutoff)]))


