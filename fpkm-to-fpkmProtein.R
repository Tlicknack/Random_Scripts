Files came from /N/dc2/projects/paramecium/Bigproject/Expression_Data ---> go into species and look for isoform fpkm

Bi:
sed 's/BIA.V1_4.1.G/BIA.V1_4.1.P/' pbi-fpkm.tab > tmp
sed 's/BIA.V1_4.1.G/BIA.V1_4.1.P/' tmp > tmp2
mv tmp2 pbi-fpkm_prot.tab

Dec:
sed 's/DEC.223.1.G/DEC.223.1.P/' pdec-fpkm.tab > tmp
sed 's/DEC.223.1.G/DEC.223.1.P/' tmp > tmp2
mv tmp2 pdec-fpkm_prot.tab

Dodec:
sed 's/DODEC.274.1.G/DODEC.274.1.P/' pdodec-fpkm.tab > tmp
sed 's/DODEC.274.1.G/DODEC.274.1.P/' tmp > tmp2
mv tmp2 pdodec-fpkm_prot.tab

Nov:
sed 's/NOV.TE.1.G/NOV.TE.M.1.P/' pnov-fpkm.tab > tmp
sed 's/NOV.TE.1.G/NOV.TE.M.1.P/' tmp > tmp2
mv tmp2 pnov-fpkm_prot.tab

Oct:
sed 's/OCT.K8.1.G/OCT.K8.1.P/' poct-fpkm.tab > tmp
sed 's/OCT.K8.1.G/OCT.K8.1.P/' tmp > tmp2
mv tmp2 poct-fpkm_prot.tab

Nov:
sed 's/SEPT.38.1.G/SEPT.38.1.P/' psept-fpkm.tab > tmp
sed 's/SEPT.38.1.G/SEPT.38.1.P/' tmp > tmp2
mv tmp2 psept-fpkm_prot.tab

Sex:
sed 's/SEX.AZ8_4.1.G/SEX.AZ8_4.1.P/' psex-fpkm.tab > tmp
sed 's/SEX.AZ8_4.1.G/SEX.AZ8_4.1.P/' tmp > tmp2
mv tmp2 psex-fpkm_prot.tab

Tet:
sed 's/TET.51.1.G/TET.51.1.P/' ptet-fpkm.tab > tmp
sed 's/TET.51.1.G/TET.51.1.P/' tmp > tmp2
mv tmp2 ptet-fpkm_prot.tab  
