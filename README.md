# Linkage

Scripts to do linkage with Merlin.

Extract SNPs for linkage:
```
Rscript scripts/Linkage/snps_for_linkage.R --cm.step 0.1
```

This will get SNPs at an interval of 0.1cm and run Merlin with non parametric linkage, while skipping all WGS and WES individuals:
```
Rscript scripts/Linkage/prepare_linkage.R --chr 22 --linkage.markers linkage_markers.csv --base.dir ~/People/PetraLiskova/ --skip J2,J3,J5,J6,J8,J9,J10,J11,J12,J13,J14,JW11,JW15,JW3,JW5 --trim PIV17,PIV19 --cm.step 1
```
```
Rscript scripts/Linkage/prepare_linkage.R --chr 22 --linkage.markers linkage_markers.csv --base.dir ~/People/PetraLiskova/ --skip J2,J3,J5,J6,J8,J9,J10,J11,J12,J13,J14,JW11,JW15,JW3,JW5 --trim PIV17,PIV19 --cm.step 1
```
Add individual 397:
```
Rscript scripts/Linkage/prepare_linkage.R --chr 8 --linkage.markers linkage_markers.csv --base.dir ~/People/PetraLiskova/ --trim PIV17 --skip J6,PIII5,PIV22,PV19,JW15 --cm.step .1
```

The WGS samples were added in to see if we could boost the LOD score:
```
Rscript scripts/Linkage/prepare_linkage.R --chr 22 --linkage.markers linkage_markers.csv --base.dir ~/People/PetraLiskova/ --skip JW11,JW15,JW3 --cm.skip 0.1
```




