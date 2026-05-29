#!/usr/bin/python
# Prompts for parameters and creates in.hybrid for LAMMPS
# LAMMPS script obtains desired system and generates position
# dump files (minimize.lammpstrj, deform.lammpstrj, dump1.lammpstrj)
# as well as stress data (str2.str)

import sys
import os
import math
import datetime
import shutil


def inputdata():
  sys.path.insert(1,'.')
  try: 
    with open('param.py'): pass
  except IOError:
    sys.exit(1)

  from param import pbsname, rseed, diameter, grftdens
  from param import ltether, attraction, matrixattraction
  from param import systemdens, wtprcnt, runtime
  from param import lammpsreplicate, lmatrix

  if attraction > matrixattraction:
    matrixattraction = attraction
  if diameter == 0:
    print ("Pure polymer system requested")
    radius = mass = grftdens = ltether = ntether = wtprcnt = 0 
    diffmatrix = 0
    attraction = 0
    nmatrix = 25
    lammpsreplicate = '1 1 1'
  else:
    radius = diameter/2
    mass = math.pow(diameter,3)
  repli = [int(n) for n in lammpsreplicate.split()]
  replicates = repli[0]*repli[1]*repli[2]

  lammpsinput = 'lammps.input'
  o = open(lammpsinput, 'w') 
  
  if diameter == 0 and matrixattraction == 0:
    type1 = 1
    type2 = 1
    type3 = 1
    type4 = 1
    maxtype = 1
    NPexist = 'No'
  elif diameter == 0 and matrixattraction > 0:
    type1 = 1
    type2 = 1
    type3 = 1
    type4 = 1
    type5 = 1
    type6 = 2
    maxtype = 2
    NPexist = 'No'
    WCAtable = '/home/dan/Desktop/wca_rep_08_11.table'
    shutil.copy(WCAtable,'.')
  elif grftdens == 0 and matrixattraction > 0:
    type1 = 1
    type2 = 2
    type3 = 1
    type4 = 2
    type5 = 1
    type6 = 3
    maxtype = 3
    NPexist = 'Yes'
    WCAtable = '/home/dan/Desktop/wca_rep_08_11.table'
    shutil.copy(WCAtable,'.')
  elif grftdens != 0 and attraction == 0 and matrixattraction ==0: 
    type1 = 1 # NP Shell
    type2 = 2 # NP Grafted Chains
    type3 = 3 # NP Center
    type4 = 4 # Matrix Chains
    type5 = 5
    maxtype = 5
    NPexist = 'Yes'
  elif attraction > 0:
    type1 = 1 # NP Shell
    type2 = 2 # NP Grafted Chains
    type3 = 3 # NP Center
    type4 = 4 # Matrix Chains
    type5 = 5 # Tether Attraction Unit
    type6 = 6 # Matrix Attraction Unit
    maxtype = 6
    NPexist = 'Yes'
    WCAtable = '/home/dan/Desktop/wca_rep_08_11.table'
    shutil.copy(WCAtable,'.')
  elif grftdens != 0 and matrixattraction > 0:
    type1 = 1
    type2 = 2
    type3 = 3
    type4 = 4
    type5 = 1
    type6 = 5
    maxtype = 5
    NPexist = 'Yes'
    #WCAtable = '/home/dan/Desktop/wca_rep_08_11.table'
    #shutil.copy(WCAtable,'.')
  else:
    type1 = 1 # NP Shell
    type2 = 2 # NP Grafted Chains
    type3 = 1 # NP Center
    type4 = 2 # Matrix Chains
    maxtype = 2
    NPexist = 'Yes'
  

  o.write("LAMMPS FENE input data file\n\n")
  if NPexist == 'Yes':
    ntether = int(round(math.pi*diameter*diameter*grftdens))
    nmatrix = int(math.ceil((math.pow(diameter,3)/(wtprcnt/100)-ntether*ltether-math.pow(diameter,3))/lmatrix))
    expectedwtprcnt = replicates*math.pow(diameter,3)/(replicates*nmatrix*lmatrix+replicates*ntether*ltether+replicates*math.pow(diameter,3))
    testwtprcnt = replicates*math.pow(diameter,3)/(replicates*(nmatrix-1)*lmatrix+replicates*ntether*ltether+replicates*math.pow(diameter,3))
    truenmatrix = int(round((replicates*math.pow(diameter,3)/(wtprcnt/100)-replicates*ntether*ltether-replicates*math.pow(diameter,3))/float(lmatrix)))
    diffmatrix = nmatrix*replicates - truenmatrix
    if diffmatrix > nmatrix:                      #??? 
      numdeletereps = int(math.ceil(diffmatrix / float(nmatrix)))
      extramat = diffmatrix - nmatrix
      diffmatrix = nmatrix
      newexpectedwtprcnt = replicates*math.pow(diameter,3)/(replicates*nmatrix*lmatrix-diffmatrix*lmatrix+replicates*ntether*ltether+replicates*math.pow(diameter,3))
      if math.fabs(newexpectedwtprcnt-(wtprcnt/100)) > math.fabs(testwtprcnt-(wtprcnt/100)):
        diffmatrix = 0
        extramat = 0
        numdeletereps = 0
        nmatrix -= 1
        newexpectedwtprcnt = replicates*math.pow(diameter,3)/(replicates*nmatrix*lmatrix-diffmatrix*lmatrix+replicates*ntether*ltether+replicates*math.pow(diameter,3))
        print ("Not deleting any polymers. Weight percentage may be higher than requested.")
        print ("Expected Nanoparticle Weight Percentage: %f%%\n" %(newexpectedwtprcnt*100))
      else:
        print ("\nWARNING: Not enough matrix molecules to groom competely\n\tAdapt code for futher precision\n\t(ie.remove matrix polymers before last NP)\n\tThis may require adaptation of LAMMPS code to remap bonds\n")
        print ("Expected Nanoparticle Weight Percentage: %f%%\n" %(newexpectedwtprcnt*100))
    elif diffmatrix < 0:
      print ("WARNING: Need more matrix polymers to acheive desired NP Weight %")
    elif diffmatrix == 0:
      print ("Not deleting any polymers. Weight percentage may be higher than requested.")
      print ("Expected Nanoparticle Weight Percentage: %f%%\n" %(expectedwtprcnt*100))
      extramat = 0
      numdeletereps = 0
      
    print ("Number of chain of 1st NP are %d" %(ntether))
    print ("Number of matrix chain are %d" %(nmatrix))
 
  natom = nmatrix*lmatrix+ntether*(ltether+1)+1
  if diameter == 0:
    natom -= 1
  nbond = ntether+ntether*(ltether-1)+nmatrix*(lmatrix-1)
  volume = ((replicates*nmatrix-diffmatrix) * lmatrix + replicates*ntether*ltether)/systemdens+replicates*math.pi/float(6)*math.pow(diameter,3)
  lsimbox = math.pow(volume,(1/3.))
  xlo = -100
  xhi = 100
  ylo = -50
  yhi = 50
  zlo = -100
  zhi = 100
  
  XLO = xlo
  XHI = xlo+((xhi-xlo)*repli[0])
  YLO = ylo
  YHI = ylo+((yhi-ylo)*repli[1])
  ZLO = zlo 
  ZHI = zlo+((zhi-zlo)*repli[2])

  xxo = 0
  yyo = 0
  zzo = 0

  o.write(str(natom)+" atoms\n")
  o.write(str(nbond)+" bonds\n\n")
  o.write(str(maxtype)+" atom types\n1 bond types\n\n")
  o.write(str(xlo)+' '+str(xhi)+ " xlo xhi\n")
  o.write(str(ylo)+' '+str(yhi)+ " ylo yhi\n")
  o.write(str(zlo)+' '+str(zhi)+ " zlo zhi\n\n")
  o.write("Masses\n\n")
  if diameter == 0 and matrixattraction == 0:
    o.write("1  1\n\n")
    i = 1
    molid = 1
    bead1 = 1
    bead2 = 2
  elif diameter == 0 and matrixattraction > 0:
    o.write("1  1\n")
    o.write("2  1\n\n")
    i = 1
    molid = 1
    bead1 = 1
    bead2 = 2
  elif attraction == 0 and grftdens !=0 and matrixattraction ==0:
    o.write("1  1\n2  1\n3  "+str(mass)+"\n4  1\n5  1\n\n")
  elif attraction > 0 and grftdens != 0:
    o.write("1  1\n2  1\n3  "+str(mass)+"\n4  1\n5  1\n6  1\n\n")
  elif grftdens == 0 and matrixattraction == 0:
    o.write("1  "+str(mass)+"\n2  1\n\n")
  elif grftdens == 0 and matrixattraction > 0:
    o.write("1  "+str(mass)+"\n2  1\n3  1\n\n")
  elif grftdens != 0 and matrixattraction > 0:
    o.write("1  1\n2  1\n3  "+str(mass)+"\n4  1\n5  1\n\n")
  o.write("Atoms\n\n")
  if diameter != 0:
    o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (1,1,type3,xxo,yyo,zzo))
    i = 2
    molid = 2
    bead1 = 2
    bead2 = 3
 
  NPlammpslist = [1]
  for k in range(1,ntether+1):
    h = -1 + 2*(k-1)/float(ntether-1)
    theta = math.acos(h)
    if k in (1,ntether):
      phi = 0
    else:
      phi = (phi+3.6/float(math.pow(ntether*(1-h*h),.5))) % 6.285714
    xx = math.cos(phi)*math.sin(theta) * radius
    yy = math.sin(phi)*math.sin(theta) * radius
    zz = math.cos(theta) * radius
    o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,2,type1,xx,yy,zz))
    NPlammpslist.append(i)     # add i to NPlammpslist
    i += 1
    molid += 1
    for kk in range(1,ltether+1):
      factor = radius/float(radius+kk)
      xxb = xx/factor
      yyb = yy/factor
      zzb = zz/factor
      if kk > (ltether-attraction):
        o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,molid,type5,xxb,yyb,zzb))
      else:
        o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,molid,type2,xxb,yyb,zzb))
      i += 1
  
  lastNPtet = i-1
  if grftdens == 0:
    molid = 2
  else:
    molid += 1
  xx = 25
  width = int(round(math.pow(nmatrix,.5)))
  height = int(math.pow(nmatrix,.5))
  remainder = nmatrix - width*height
 
  for n in range(width):
    yy = -25
    for m in range(height):
      zz = 1
      for j in range(lmatrix):
        if (lmatrix/2-matrixattraction) <= j < (lmatrix/2) or j >= (lmatrix-matrixattraction):   
          o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,molid,type6,xx,yy,zz))
        else:  
          o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,molid,type4,xx,yy,zz))
        zz += 1
        i += 1
      molid += 1
      yy += 1
    xx += 1
  for n in range(remainder):
    zz = 1                          # remainder x does not change but lay up on y
    for j in range(lmatrix):
      if (lmatrix/2-attraction) <= j < (lmatrix/2) or j >= (lmatrix-attraction):   
        o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,molid,type6,xx,yy,zz))
      else:
        o.write("%8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,molid,type4,xx,yy,zz))
      zz += 1
      i += 1
    molid += 1
    yy +=1

  lastmatid = (i-1)*replicates
  o.write("\nBonds\n\n")
  bondid = 1
  bondtype = 1
  for i in range(ntether):
    for j in range(ltether):
      o.write("%8i %8i %8i %8i\n" % (bondid, bondtype, bead1, bead2))
      bead1 = bead2
      bead2 += 1
      bondid += 1
    bead1 += 1
    bead2 += 1

  for i in range(nmatrix):
    for j in range(lmatrix-1): 
      o.write("%8i %8i %8i %8i\n" % (bondid, bondtype, bead1, bead2))
      bead1 = bead2
      bead2 += 1
      bondid += 1
    bead1 += 1
    bead2 += 1
  
  o.close()
  print ('END OF PROGRAMMING \n') 

def main():
  if len(sys.argv) != 1:
    print ('usage: No Flags')
    sys.exit(1)
  
  inputdata()

if __name__ == '__main__':
  main()
