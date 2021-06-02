#
#    This file is part of Lep-MAP3.
#
#    Lep-MAP3 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Lep-MAP3 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Lep-MAP3.  If not, see <http://www.gnu.org/licenses/>.
#
#        Copyright (C) 2021 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki


#QTL LOD plot

library(scales)

p=read.table("qtlphenotypes.txt", sep="\t")

png("qtl.png",600,600)

mh = list()
mhc = list()

mh2 = list()

for (chr in 1:1) {#edit to work for multiple chromosomes, e.g. 1:21
	d=read.table(paste0("qtldata",chr,".12"), sep="\t")
	n=nrow(d)
	m=ncol(d)

	#different families
	fam=unique(t(d[1,5:m]))

	gen=d[7:n,5:m]

	data1 = 2  * (gen == "1 1") + 2 * (gen == "1 2") - 1
	data2 = 2  * (gen == "1 1") + 2 * (gen == "2 1") - 1
	data3 = 2  * (gen == "1 1") + 2 * (gen == "2 2") - 1

	numMarkers = n - 6
	numFamilies = length(fam)
	print(paste0(numMarkers, " markers and ", numFamilies, " families"))

	sex = p$V1 - 1
	weight = p$V2

	ret = matrix(0, nrow = numMarkers, ncol = numFamilies)

	famIndex = 1
	for (f in fam) {
		fi = (d[1,5:m] == f)

		d1 = data1[,fi]
		d2 = data2[,fi]
		d3 = data3[,fi]
		s  = sex[fi]
		w = weight[fi]

		m0 = glm(w ~ s)
		l0 = logLik(m0)
	
		#cat(s, "\n")

		m1 = glm(w ~ s + d1[1,] + d2[1,] + d3[1,])
		l1 = logLik(m1)
		ret[1, famIndex] = l1 - l0
		for (i in 2:numMarkers) {
			if (!all(d1[i - 1,] == d1[i,]) | !all(d2[i - 1,] == d2[i,])) {
				m1 = glm(w ~ s + d1[i,] + d2[i,] + d3[i,])
				l1 = logLik(m1)
			}
			ret[i, famIndex] = l1 - l0
		}
		famIndex = famIndex + 1
	}
	rs = rowSums(ret)
#	plot(rs, ylim=c(20, max(max(rs), qchisq(.95, df=numFamilies * 3))))
	mh[[chr]] = rs
	mhc[[chr]] = rep((chr %% 2) + 1, length(rs))

	mh2[[chr]] = ret

}
y=Reduce(c, mh)
plot(y, col=alpha(Reduce(c, mhc), 0.5), pch=19, xlab="Marker", ylab="LOD", cex.lab=1.5, ylim=c(0,max(y)))

#plot significant lines, could be loaded from a file
#abline(h=103.10) #5%, from qtlPerm.R
#abline(h=109.76) #1%, from qtlPerm.R
dev.off()
warnings()
