#!/usr/bin/python3

import sys
import math

def CxG1Search(m):
 if m < 0 or m > 5:
  print("Bullet mach out of bound 0-5")
  return(0)

 MCxG1=[[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.9,0.925,0.95,0.975,1,1.025,1.05,1.075,1.1,1.125,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.2,4.4,4.6,4.8,5],[0.2629,0.2558,0.2487,0.2413,0.2344,0.2278,0.2214,0.2155,0.2104,0.2061,0.2032,0.2020,0.2034,0.2165,0.2230,0.2313,0.2417,0.2546,0.2706,0.2901,0.3415,0.3734,0.4084,0.4448,0.4805,0.5135,0.5427,0.5677,0.5883,0.6053,0.6191,0.6393,0.6518,0.6589,0.6621,0.6625,0.6607,0.6573,0.6528,0.6474,0.6413,0.6347,0.628,0.6210,0.6141,0.6072,0.6003,0.5934,0.5867,0.5804,0.5743,0.5685,0.5630,0.5577,0.5527,0.5481,0.5438,0.5397,0.5325,0.5264,0.5211,0.5168,0.5133,0.5105,0.5084,0.5067,0.5054,0.5040,0.5030,0.5022,0.5016,0.5010,0.5006,0.4998,0.495,0.4992,0.4990,0.4988]]

 l=len(MCxG1)
 c=len(MCxG1[0])

# print("Nombre de couple Mac Cx G1 : ", c)
# print("G1 Cx by Mach number ")
# for i in range(c):
#  print("Mach :",MCxG1[0][i]," Cx G1 :",MCxG1[1][i])

# return(-1)

# print("Searching for CX G1 for bullet with Mach number = ",m)

 i=0
 End=0
 m_inf=MCxG1[0][0]
 Cx_inf=MCxG1[1][0]
 m_sup=MCxG1[0][0]
 Cx_sup=MCxG1[1][0]
 CxG1=0

 while (i<c and End!=1):
  if i>0:
   if m >= MCxG1[0][i-1] and m <= MCxG1[0][i] :
    m_inf=MCxG1[0][i-1]
    Cx_inf=MCxG1[1][i-1]
    m_sup=MCxG1[0][i]
    Cx_sup=MCxG1[1][i]
    End=1
   else:
    i=i+1
  else:
   i=i+1

# print("i =",i)
# print("End =",End)
# print("m_inf =",m_inf)
# print("Cx_inf =",Cx_inf)
# print("m_sup =",m_sup)
# print("Cx_sup =",Cx_sup)

 if m == m_inf:
  CxG1=float(Cx_inf)
#  print("Cx inf :",CxG1)
 elif m == m_sup:
   CxG1=float(Cx_sup)
#   print("Cx sup :",CxG1)
 else:
   a=float((Cx_sup-Cx_inf)/(m_sup-m_inf))
#   print("a :",a)
   b=float(Cx_inf-a*m_inf)
#   print("b :",b)
   CxG1=float(a*m+b)
#   print("Cx G1 AL :",CxG1)

 return CxG1

def CxG7Search(m):
 if m < 0 or m > 5:
  print("Bullet mach out of bound 0-5")
  return(0)

 MCxG7=[[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.725,0.75,0.775,0.8,0.825,0.875,0.9,0.925,0.95,0.975,1,1.025,1.05,1.075,1.1,1.125,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.2,4.4,4.6,4.8,5],[0.1198,0.1197,0.1196,0.1194,0.1193,0.1194,0.1194,0.1194,0.1193,0.1193,0.1194,0.1193,0.1194,0.1197,0.1202,0.1207,0.1215,0.1226,0.1242,0.1266,0.1368,0.1464,0.166,0.2054,0.2993,0.3803,0.4015,0.4043,0.4034,0.4014,0.3987,0.3955,0.3884,0.381,0.3732,0.3657,0.358,0.344,0.3376,0.3315,0.326,0.3209,0.316,0.3117,0.3078,0.3042,0.301,0.298,0.2951,0.2922,0.2892,0.2864,0.2835,0.2807,0.2779,0.2752,0.2725,0.2697,0.267,0.2643,0.2615,0.2588,0.2561,0.2533,0.2506,0.2479,0.2451,0.2424,0.2368,0.2313,0.2258,0.2205,0.2154,0.2106,0.206,0.2017,0.1975,0.1935,0.1861,0.1793,0.173,0.1672,0.1618]]

 l=len(MCxG7)
 c=len(MCxG7[0])

# print("Nombre de couple Mac Cx G7 : ", c)
# print("G7 Cx by Mach number ")
# for i in range(c):
#  print("Mach :",MCxG7[0][i]," Cx G7 :",MCxG7[1][i])

# return(-1)

# print("Searching for CX G7 for bullet with Mach number = ",m)

 i=0
 End=0
 m_inf=MCxG7[0][0]
 Cx_inf=MCxG7[1][0]
 m_sup=MCxG7[0][0]
 Cx_sup=MCxG7[1][0]
 CxG7=0

 while (i<c and End!=1):
  if i>0:
   if m >= MCxG7[0][i-1] and m <= MCxG7[0][i] :
    m_inf=MCxG7[0][i-1]
    Cx_inf=MCxG7[1][i-1]
    m_sup=MCxG7[0][i]
    Cx_sup=MCxG7[1][i]
    End=1
   else:
    i=i+1
  else:
   i=i+1

# print("i =",i)
# print("End =",End)
# print("m_inf =",m_inf)
# print("Cx_inf =",Cx_inf)
# print("m_sup =",m_sup)
# print("Cx_sup =",Cx_sup)

 if m == m_inf:
  CxG7=float(Cx_inf)
#  print("Cx inf :",CxG7)
 elif m == m_sup:
   CxG7=float(Cx_sup)
#   print("Cx sup :",CxG7)
 else:
   a=float((Cx_sup-Cx_inf)/(m_sup-m_inf))
#   print("a :",a)
   b=float(Cx_inf-a*m_inf)
#   print("b :",b)
   CxG7=float(a*m+b)
#   print("Cx G7 AL :",CxG7)

 return CxG7


BCG1=float(0.421)
BCG7=float(0.327)
m = float(0)
CxG1 = float(-1)
CxG7 = float(-7)

m = float(sys.argv[1])

print("Bullet Mach : ",m)

CxG1=CxG1Search(m)
#print("Bullet Cx G1 : ",CxG1)

CxG7=CxG7Search(m)
#print("Bullet Cx G7 : ",CxG7)


print("Bullet Cx G1 : %0.4f" % CxG1)
print("Bullet Cx G7 : %.04f" % CxG7)



