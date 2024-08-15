#!/usr/bin/python3

import sys
import math

print("==================================================")
print(" PBS stand for Python Ballistic Solver")
print(" PBS is an Open Source Balistic Software ")
print(" Written in Python by Fabien FIGUERAS (he/him)")
print(" v1.00 was released in 2024 ")
print(" Current Version is v1.05 2024")
print(" Call example Python3 ./PBS-v100.py 999.99 0.0001 30.1")
print(" Where param1 is the caliber (inch) ")
print(" Where param2 is the bullet mass (gr) ")
print(" Where param3 is the Muzzle Speed (m/s) ")
print(" Where param4 is the Shooting distance ")
print(" Where param5 is the Azimut (shooting angle relative to the North ")
print(" Where param6 is the Coriolis Option (Y/N) ")
print(" Where param7 is the Average Wind Speed ")
print(" Where param8 is the Wind Speed related to shooting direction (hour) ")
print(" Where param9 is the Spind Drift Option (Y/N) ")
print(" Where param10 is the time increment (s) ")
print(" Where param11 is Zeroing (Y/N) ")

print(" Sources available in GitHub : Link to be added ")
print("==================================================")


SightHeight=float(0.065)
Rifle="Tikka T3X"
Twist=float(11)
Twist_C="Right"

d_i=float(sys.argv[1])
T_cal=float(Twist/d_i)
m_gr=float(sys.argv[2])
bl_cm=float(3.437)
bl_i=float(bl_cm/2.54)
bl_cal=float(bl_i/d_i)
MuzzleSpeed=float(sys.argv[3])
alpha_d=float(0.083871)
D_Tir=float(sys.argv[4])
Az_d=float(sys.argv[5])
Co=sys.argv[6]
Ws=float(sys.argv[7])
Wa_o=float(sys.argv[8])
SD_p=sys.argv[9]
h=float(sys.argv[10])
Zero_C=sys.argv[11]

E_c=float(0)
W_c=float(0)

tph=float(0)
Vxtph=float(0)
Vytph=float(0)
Vztph=float(0)
SD=float(0)
Sg=float(0)
WDx=float(0)
WDy=float(0)
Wsx=float(0)
Wsz=float(0)
V0=float(MuzzleSpeed)
Lat_d=float(46)
Lat_m=float(22)
Lat_s=float(25)
Lat_D=float(Lat_d+(Lat_m/60)+(Lat_s/(60*60)))
NS="North"
Az=float((Az_d*math.pi/180))
Vmod=float(0)
Xtph=float(0)
Ytph=float(0)
Ztph=float(0)
Gamma=float((1+(1/2.48)))

Alt_ICAO=float(0)
P_ICAO=float(101325)
T_C_ICAO=float(15)
T_K_ICAO=float(T_C_ICAO+273.15)
RH_ICAO=0
Rho_ICAO=float(0.003483*((1 + RH_ICAO) / (1 + 1.6078 * RH_ICAO))* ( P_ICAO/ (T_C_ICAO + 273.15) ))
PVS_ICAO=float(6.1078*math.pow(10,((7.5*T_K_ICAO-2048.625)/(T_K_ICAO-35.85))))
PV_ICAO=float(RH_ICAO*PVS_ICAO)
TV_ICAO=float(T_K_ICAO/(1-0.3785*(PV_ICAO/P_ICAO)))

SoundSpeed_ICAO=float(math.sqrt(Gamma*P_ICAO/Rho_ICAO))

Sg_ICAO=float((30*m_gr)/(math.pow(T_cal,2)*math.pow(d_i,3)*bl_cal*(1+math.pow(bl_cal,2))*(math.pow((V0/853.4),(1/3))*((T_K_ICAO*101325)/(288.15*P_ICAO)))))

Alt_Zero=float(440)
P_Zero=float(96330)
T_C_Zero=float(26)
T_K_Zero=float(T_C_Zero+273.15)
RH_Zero=0.66

Alt_Zero=float(0)
P_Zero=float(101325)
T_C_Zero=float(15)
T_K_Zero=float(T_C_Zero+273.15)
RH_Zero=0

Rho_Zero=float(0.003483*((1 + RH_Zero) / (1 + 1.6078 * RH_Zero))* ( P_Zero/ (T_C_Zero + 273.15) ))
PVS_Zero=float(6.1078*math.pow(10,((7.5*T_K_Zero-2048.625)/(T_K_Zero-35.85))))
PV_Zero=float(RH_Zero*PVS_Zero)
TV_Zero=float((T_K_Zero/(1-0.3785*(PV_Zero/P_Zero))))

Alt=float(440)
P=float(96330)
T_C=float(26)
T_K=float(T_C+273.15)
RH=0.66
Rho=float(0.003483*((1 + RH) / (1 + 1.6078 * RH))* ( P/ (T_C + 273.15) ))
PVS=float(6.1078*math.pow(10,((7.5*T_K-2048.625)/(T_K-35.85))))
PV=float(RH*PVS)
TV=float((T_K/(1-0.3785*(PV/P))))

Omega=float(2*math.pi/((23+(56/60)+(4.09/(60*60)))*(60*60)))

D_zero=float(100)
Prec=0.001

#Test d'un angle pour une distance donnée
def TestAngle(Alpha,D_zero,Prec):
 global Ytph,h,MuzzleSpeed,NS,Az_d,Az,Co,Ws,Wa_o,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle,d_i,bl_cal,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO,Alt_Zero,P_Zero,T_C_Zero,T_K_Zero,RH_Zero,Rho_Zero,PVS_Zero,PV_Zero,TV_Zero

 ta=0

 RK(D_zero,h,MuzzleSpeed,NS,Az_d,Az,Co,Ws,Wa_o,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle,d_i,bl_cal,Alpha,Alt_Zero,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P_Zero,T_C_Zero,T_K_Zero,RH_Zero,Rho_Zero,PVS_Zero,PV_Zero,TV_Zero,Omega,SoundSpeed_ICAO)

 if (Ytph>Prec):
  ta=-1
 if (Ytph<0):
  ta=+1

 return(ta)

#Calcul de l'angle Alpha pour que Y=0 pour une distance D_zero
def Zeroing(D_zero,Prec):
 Alpha=(100/3)
# print("Initial Angle (deg): ",math.degrees(lza))
 ta=1
 while (ta!=0):
  ta=int(TestAngle(Alpha,D_zero,Prec))
#  print(ta, " returned")
#  te=input('Type Enter to resume...')
  if (ta<0):
   Alpha=Alpha-(Alpha/2)
#   print("-1 Modified Angle (deg): ",math.degrees(lza))
  if (ta>0):
   Alpha=Alpha+(Alpha/2)
#   print("+1 Modified Angle (deg): ",math.degrees(lza))
# print("Zeroing Angle (deg): ",math.degrees(lza))
 return(Alpha)

def SpinDrift(Sg,ToF):
 Dg=float(3.175*(Sg+1.2)*math.pow(ToF,1.83))
 return(Dg)

def CdG1Search(m):
 if m < 0 or m > 5:
  print("Bullet mach out of bound 0-5")
  return(0)

 MCxG1=[[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.9,0.925,0.95,0.975,1,1.025,1.05,1.075,1.1,1.125,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.2,4.4,4.6,4.8,5],[0.2629,0.2558,0.2487,0.2413,0.2344,0.2278,0.2214,0.2155,0.2104,0.2061,0.2032,0.2020,0.2034,0.2165,0.2230,0.2313,0.2417,0.2546,0.2706,0.2901,0.3415,0.3734,0.4084,0.4448,0.4805,0.5135,0.5427,0.5677,0.5883,0.6053,0.6191,0.6393,0.6518,0.6589,0.6621,0.6625,0.6607,0.6573,0.6528,0.6474,0.6413,0.6347,0.628,0.6210,0.6141,0.6072,0.6003,0.5934,0.5867,0.5804,0.5743,0.5685,0.5630,0.5577,0.5527,0.5481,0.5438,0.5397,0.5325,0.5264,0.5211,0.5168,0.5133,0.5105,0.5084,0.5067,0.5054,0.5040,0.5030,0.5022,0.5016,0.5010,0.5006,0.4998,0.495,0.4992,0.4990,0.4988]]

 l=len(MCxG1)
 c=len(MCxG1[0])

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

 if m == m_inf:
  CxG1=float(Cx_inf)

 elif m == m_sup:
   CxG1=float(Cx_sup)

 else:
   a=float((Cx_sup-Cx_inf)/(m_sup-m_inf))

   b=float(Cx_inf-a*m_inf)

   CxG1=float(a*m+b)

 return(CxG1)

def CdG7Search(m):
 if m < 0 or m > 5:
  print("Bullet mach out of bound 0-5")
  return(0)

 MCxG7=[[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.725,0.75,0.775,0.8,0.825,0.875,0.9,0.925,0.95,0.975,1,1.025,1.05,1.075,1.1,1.125,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.2,4.4,4.6,4.8,5],[0.1198,0.1197,0.1196,0.1194,0.1193,0.1194,0.1194,0.1194,0.1193,0.1193,0.1194,0.1193,0.1194,0.1197,0.1202,0.1207,0.1215,0.1226,0.1242,0.1266,0.1368,0.1464,0.166,0.2054,0.2993,0.3803,0.4015,0.4043,0.4034,0.4014,0.3987,0.3955,0.3884,0.381,0.3732,0.3657,0.358,0.344,0.3376,0.3315,0.326,0.3209,0.316,0.3117,0.3078,0.3042,0.301,0.298,0.2951,0.2922,0.2892,0.2864,0.2835,0.2807,0.2779,0.2752,0.2725,0.2697,0.267,0.2643,0.2615,0.2588,0.2561,0.2533,0.2506,0.2479,0.2451,0.2424,0.2368,0.2313,0.2258,0.2205,0.2154,0.2106,0.206,0.2017,0.1975,0.1935,0.1861,0.1793,0.173,0.1672,0.1618]]

 l=len(MCxG7)
 c=len(MCxG7[0])

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

 if m == m_inf:
  CxG7=float(Cx_inf)

 elif m == m_sup:
   CxG7=float(Cx_sup)

 else:
   a=float((Cx_sup-Cx_inf)/(m_sup-m_inf))

   b=float(Cx_inf-a*m_inf)

   CxG7=float(a*m+b)


 return(CxG7)


def Fx(t,Vx,Vy,Vz,Ox,Oy,Oz,LCo,k,m,g):

    if LCo == "N":
        Fx=float((-k/m)*abs(Vx)*Vx)
    else:
        Fx=float((-k/m)*abs(Vx)*Vx-2*(Oy*Vz-Oz*Vy))

    return(Fx)

def Fy(t,Vx,Vy,Vz,Ox,Oy,Oz,LCo,k,m,g):
    if LCo == "N":
        Fy=float((-k/m)*abs(Vy)*Vy-g)
    else:
        Fy=float((-k/m)*abs(Vy)*Vy-g-2*(-(Ox*Vz-Oz*Vx)))

    return(Fy)

def Fz(t,Vx,Vy,Vz,Ox,Oy,Oz,LCo,k,m,g):

    if LCo == "N":
        Fz=float((-k/m)*abs(Vz)*Vz)

    else:
        Fz=float((-k/m)*abs(Vz)*Vz-2*(Ox*Vy-Oy*Vx))

    return(Fz)

def CdSearch(m,BC_G1,Sd):
    Cd_G1=CdG1Search(m)
    Cd=float(Sd*(Cd_G1/BC_G1))
    return(Cd)

def RK(D_Tir,h,MuzzleSpeed,NS,Az_d,Az,Co,Ws,Wa_o,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle,d_i,bl_cal,alpha_d,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO):

    global E_c,W_c,tph,Vxtph,Vytph,Vztph,Sg,WDx,WDy,Wsx,Wsz,V0,Lat_D,Vmod,Xtph,Ytph,Ztph

    d_mm=float(d_i*25.4)
    d=float(d_mm/1000)
    m_g=float(m_gr*0.0647989)
    m=float(m_g/1000)

    A=float(math.pi*(d*d)/4)

    Sd=float((m_gr/7000)/(d_i*d_i))

#    print("bullet mass (g) : ",round(m_g,2))
#    print("bullet diameter (mm) : ",round(d_mm,2))
#    print("bullet diameter (m) : ",round(d,5))
#    print("bullet length (cm) : ",round(bl_cm,3)," (inch) : ",round(bl_i,3))
#    print("bullet mass (kg) : ",round(m,5))
#    print("Sd (pound/inch) :",round(Sd,3))

    Lat=float(Lat_D*math.pi/180)

    g=float(9.780327*(1+5.3024*(0.001)*math.sin(Lat)*math.sin(Lat)-5.8*(0.000001)*math.sin(2*Lat)*math.sin(2*Lat-3.086*(0.0000001)*Alt)))

#    print("Lattitude (°) :",round(Lat_D,2),NS)
#    print("Lattitude (rad) :",round(Lat,3))
#    print("Earth Acceleration (N/m) :",round(g,2))
#    print("Gamma : ",round(Gamma,3))

    alpha=float(alpha_d*math.pi/180)
#    print("Initial Horizontal Shooting Angle - Alpha 0 (°) :",alpha_d," (rad) : ",round(alpha,3))

    Ox=float(Omega*math.cos(Lat)*math.cos(Az))
    Oy=float(Omega*math.sin(Lat))
    Oz=float(-Omega*math.cos(Lat)*math.sin(Az))
    ModOmega=math.sqrt(Ox*Ox+Oy*Oy+Oz*Oz)

#    print("Ox : ",Ox," Oy : ",Oy," Oz : ",Oz)
#    print("Omega Module (rad/s) : ",ModOmega)

#    print("======================")
#    print("ICAO BC_G1 value based on current speed ")
#    print("Source: Bullet Builder - Sierra")
#    print("======================")

    V0=MuzzleSpeed
    Vx=float(V0*math.cos(alpha))
    Vy=float(V0*math.sin(alpha))
    Vz=float(0)
    V=float(math.sqrt(Vx*Vx+Vy*Vy+Vz*Vz))
#    print("Initial Speed module V (m/s) :",round(V,2))

    if V>=525:
        BC_G1_ICAO=float(0.525)
    if V>=689:
        BC_G1_ICAO=float(0.533)
    if V<525:
        BC_G1_ICAO=float(0.515)

    Mach_ICAO=float(V/SoundSpeed_ICAO)

    Cd_G1_ICAO=CdG1Search(Mach_ICAO)
    Cd_G7_ICAO=CdG7Search(Mach_ICAO)
    BC_G7_ICAO=float(BC_G1_ICAO/(Cd_G1_ICAO/Cd_G7_ICAO))

    Cd_ICAO=float(Sd*(Cd_G7_ICAO/BC_G7_ICAO))

#    print("BC_G1 ICAO (no unit) :",round(BC_G1_ICAO,3))
#    print("ICAO Mach Number :",round(Mach_ICAO,2))
#    print("Cd_G1 ICAO (no unit) :",round(Cd_G1_ICAO,3))
#    print("Cd_G7 ICAO (no unit) :",round(Cd_G7_ICAO,3))
#    print("BC_G7 ICAO (no unit) :",round(BC_G7_ICAO,3))
#    print("Cd = bullet Drag Coefficient ICAO (no unit) : ",round(Cd_ICAO,3))

    k_ICAO=float(0.5*Rho_ICAO*Cd_ICAO*A)
#    print("k coefficient ICAO : ",k_ICAO)

#    print("======================")
#    print("Drag Coefficient (Cd) Determination ")
#    print("======================")

    SoundSpeed=float(math.sqrt(Gamma*P/Rho))
#    print("Speed of sound (m/s) :",round(SoundSpeed,2))

#    print("======================")
#    print("BC_G1 value based on current speed & Atmo Conditions ")
#    print("Source: Bullet Builder - Sierra")
#    print("======================")

    J=float((P_ICAO/P)*math.sqrt(TV/TV_ICAO))
#    print("Impedance Coef (no unit) : ",round(J,5))

    V0=MuzzleSpeed
    Vx=float(V0*math.cos(alpha))
    Vy=float(V0*math.sin(alpha))
    Vz=float(0)
    V=float(math.sqrt(Vx*Vx+Vy*Vy+Vz*Vz))
#    print("Initial Speed module V (m/s) :",round(V,2))

    if V>=525:
        BC_G1=float(BC_G1_ICAO*J)
    if V>=689:
        BC_G1=float(BC_G1_ICAO*J)
    if V<525:
        BC_G1=float(BC_G1_ICAO*J)

    Mach=float(V/SoundSpeed)

    Cd_G1=CdG1Search(Mach)
    Cd_G7=CdG7Search(Mach)
    BC_G7=float(BC_G1/(Cd_G1/Cd_G7))

    Cd=float(Sd*(Cd_G7/BC_G7))

#    print("BC_G1 (no unit) :",round(BC_G1,3))
#    print("Mach Number :",round(Mach,2))
#    print("Cd_G1 (no unit) :",round(Cd_G1,3))
#    print("Cd_G7 (no unit) :",round(Cd_G7,3))
#    print("BC_G7 (no unit) :",round(BC_G7,3))
#    print("Cd = bullet Drag Coefficient (no unit) : ",round(Cd,3))

    Sg=float((30*m_gr)/(math.pow(T_cal,2)*math.pow(d_i,3)*bl_cal*(1+math.pow(bl_cal,2))*(math.pow((V0/853.4),(1/3))*((T_K*101325)/(288.15*P)))))

#    print("Bullet Stability Factor ",round(Sg,2))

    if Sg<1:
        print("Bullet is NOT stable")
    if Sg<1.5:
        print("Bullet Marginaly Stable")
    if Sg >=1.5:
        print("Stable Bullet")

    k=float(0.5*Rho*Cd*A)
#    print("k coefficient : ",k)

    Wa=float((Wa_o*(math.pi/2))/3)
    Wsx=-Ws*math.cos(Wa)
    Wsy=-float(0)
    Wsz=-Ws*math.sin(Wa)

#    print("wind Angle relative to shooting direction (RAD) :",round(Wa,5))
#    print("wind speed along X (m/s) :",round(Wsx,2))
#    print("wind speed along Y (m/s) :",round(Wsy,2))
#    print("wind speed along Z (m/s) :",round(Wsz,2))

#    print("===========================")
#    print("Variables Initialization")
#    print("===========================")

    t=float(0)
    tph=float(t)
    X=float(0)
    Y=float(-SightHeight)
    Z=float(0)
    V0=MuzzleSpeed
    Vx=float(V0*math.cos(alpha))
    Vy=float(V0*math.sin(alpha))
    Vz=float(0)
    tph=float(tph+h)
    V=float(math.sqrt(Vx*Vx+Vy*Vy+Vz*Vz))
    tph=float(tph+h)

#    print("time t (s) :",t)
#    print("Initial X position (m) :",X)
#    print("Initial Y position (m) :",Y)
#    print("Initial Z position (m) :",Z)
#    print("Initial Muzzle Speed V(0) (m/s) :",round(V0,2))
#    print("Initial Speed along X axis Vx (m/s) :",round(Vx,2))
#    print("Initial Speed along Y axis Vy (m/s) :",round(Vy,2))
#    print("Initial Speed along Z axis Vz (m/s) :",round(Vz,2))
#    print("Elapsed time (s) :",tph)
#    print("Initial Speed module V (m/s) :",round(V,2))

    Cd_org=Cd
    BC_G1_B=BC_G1

#    print("===========================")
#    print("Ruge-Kutta Method")
#    print("======== STEP 1 ===========")

    Kx1=float(Fx(t,Vx,Vy,Vz,Ox,Oy,Oz,Co,k,m,g))
    Ky1=float(Fy(t,Vx,Vy,Vz,Ox,Oy,Oz,Co,k,m,g))
    Kz1=float(Fz(t,Vx,Vy,Vz,Ox,Oy,Oz,Co,k,m,g))

#    print("Calculated K1 along X axis Kx1 (m/s) :",Kx1)
#    print("Calculated K1 along Y axis Ky1 (m/s) :",Ky1)
#    print("Calculated K1 along Z axis Kz1 (m/s) :",Kz1)

    Kx2=float(Fx(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,Ox,Oy,Oz,Co,k,m,g))
    Ky2=float(Fy(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,Ox,Oy,Oz,Co,k,m,g))
    Kz2=float(Fz(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,Ox,Oy,Oz,Co,k,m,g))

#    print("Calculated K2 along X axis Kx2 (m/s) :",Kx2)
#    print("Calculated K2 along Y axis Ky2 (m/s) :",Ky2)
#    print("Calculated K2 along Z axis Kz2 (m/s) :",Kz2)

    Kx3=float(Fx(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,Ox,Oy,Oz,Co,k,m,g))
    Ky3=float(Fy(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,Ox,Oy,Oz,Co,k,m,g))
    Kz3=float(Fz(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,Ox,Oy,Oz,Co,k,m,g))

#    print("Calculated K3 along X axis Kx3 (m/s) :",Kx3)
#    print("Calculated K3 along Y axis Ky3 (m/s) :",Ky3)
#    print("Calculated K3 along Z axis Kz3 (m/s) :",Kz3)

    Kx4=float(Fx(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,Ox,Oy,Oz,Co,k,m,g))
    Ky4=float(Fy(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,Ox,Oy,Oz,Co,k,m,g))
    Kz4=float(Fz(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,Ox,Oy,Oz,Co,k,m,g))

#    print("Calculated K4 along X axis Kx4 (m/s) :",Kx4)
#    print("Calculated K4 along Y axis Ky4 (m/s) :",Ky4)
#    print("Calculated K4 along Z axis Kz4 (m/s) :",Kz4)

    Vxtph=float(Vx+(h/6)*(Kx1+2*Kx2+2*Kx3+Kx4))
    Vytph=float(Vy+(h/6)*(Ky1+2*Ky2+2*Ky3+Ky4))
    Vztph=float(Vz+(h/6)*(Kz1+2*Kz2+2*Kz3+Kz4))

#    print("Calculated Speed along X axis Vx(t+h) (m/s) :",round(Vxtph,2))
#    print("Calculated Speed along Y axis Vy(t+h) (m/s) :",round(Vytph,2))
#    print("Calculated Speed along Z axis Vz(t+h) (m/s) :",round(Vztph,2))

    Vxm=float((Vx+Vxtph)/2)
    Vym=float((Vy+Vytph)/2)
    Vzm=float((Vz+Vztph)/2)

#    print("Calculated Average Speed along X axis (Vx(t)+Vx(t+h))/2 (m/s) :",round(Vxm,2))
#    print("Calculated Average Speed along Y axis (Vy(t)+Vy(t+h))/2 (m/s) :",round(Vym,2))
#    print("Calculated Average Speed along Z axis (Vz(t)+Vz(t+h))/2 (m/s) :",round(Vzm,2))

    Xtph=float(X+h*Vxm)
    Ytph=float(Y+h*Vym)
    Ztph=float(Z+h*Vzm)

#    print("Calculated X position (m) :",round(Xtph,3))
#    print("Calculated Y position (m) :",round(Ytph,3))
#    print("Calculated Z position (m) :",round(Ztph,3))

    while (X<D_Tir) :
        t=float(tph)

        Vx=float(Vxtph)
        Vy=float(Vytph)
        Vz=float(Vztph)
        Vz=float(Vztph)
        V=float(math.sqrt(Vx*Vx+Vy*Vy+Vz*Vz))

        if V>=525:
            BC_G1=float(BC_G1_ICAO*J)
        if V>=689:
            BC_G1=float(BC_G1_ICAO*J)
        if V<525:
            BC_G1=float(BC_G1_ICAO*J)

        if BC_G1_B!=BC_G1:

            BC_G1_B=BC_G1

            Mach=float(V/SoundSpeed)

            Cd_G1=CdG1Search(Mach)

            Cd_G7=CdG7Search(Mach)

            BC_G7=float(BC_G1/(Cd_G1/Cd_G7))

            Cd=float(Sd*(Cd_G7/BC_G7))

            if Cd_org != Cd :

                Cd_org=Cd

            k=float(0.5*Rho*Cd*A)

        Kx1=float(Fx(t,Vx,Vy,Vz,Ox,Oy,Oz,Co,k,m,g))
        Ky1=float(Fy(t,Vx,Vy,Vz,Ox,Oy,Oz,Co,k,m,g))
        Kz1=float(Fz(t,Vx,Vy,Vz,Ox,Oy,Oz,Co,k,m,g))

        Kx2=float(Fx(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,Ox,Oy,Oz,Co,k,m,g))
        Ky2=float(Fy(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,Ox,Oy,Oz,Co,k,m,g))
        Kz2=float(Fz(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,Ox,Oy,Oz,Co,k,m,g))

        Kx3=float(Fx(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,Ox,Oy,Oz,Co,k,m,g))
        Ky3=float(Fy(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,Ox,Oy,Oz,Co,k,m,g))
        Kz3=float(Fz(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,Ox,Oy,Oz,Co,k,m,g))

        Kx4=float(Fx(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,Ox,Oy,Oz,Co,k,m,g))
        Ky4=float(Fy(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,Ox,Oy,Oz,Co,k,m,g))
        Kz4=float(Fz(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,Ox,Oy,Oz,Co,k,m,g))

        Vxtph=float(Vx+(h/6)*(Kx1+2*Kx2+2*Kx3+Kx4))
        Vytph=float(Vy+(h/6)*(Ky1+2*Ky2+2*Ky3+Ky4))
        Vztph=float(Vz+(h/6)*(Kz1+2*Kz2+2*Kz3+Kz4))

        Vxm=float((Vx+Vxtph)/2)
        Vym=float((Vy+Vytph)/2)
        Vzm=float((Vz+Vztph)/2)

        X=float(Xtph)
        Xtph=float(X+h*Vxm)
        Y=float(Ytph)
        Ytph=float(Y+h*Vym)
        Z=float(Ztph)
        Ztph=float(Z+h*Vzm)

        tph=float(tph+h)

print("======================")
print("Rifle variables ")
print("======================")

print("Rifle Brand and Model : ",Rifle)
print("Barrel Twist (inch) 1 :",int(Twist),Twist_C)
print("Distance between Sight and Barrel (m) :",SightHeight," (cm) :",(SightHeight*100))

print("======================")
print("Bullet ")
print("======================")

print("bullet diameter (inch) : ",d_i)
print("bullet mass (gr) : ",int(m_gr))
print("Muzzle Speed (m/s) :",round(MuzzleSpeed,0))

print("======================")
print("Localization ")
print("======================")

print("Latitude ",Lat_d," ° ",Lat_m," min ",Lat_s," s")

print("======================")
print("ICAO Standard Athmosphere ")
print("======================")

print("Altitude (m) :",int(Alt_ICAO))
print("Air pressure (Pa) : ",P_ICAO)
print("Air Temperature (°C) : ",T_C_ICAO,"(°K) : ",T_K_ICAO)
print("Air Relative Humidity (%) : ",RH_ICAO)
print("Wet Air Volumic Mass ICAO (kg/m3) : ",round(Rho_ICAO,3))
print("Saturated Vapor Pressure ICAO (Pa) : ",round(PVS_ICAO,5))
print("Vapor Pressure ICAO (Pa) : ",round(PV_ICAO,2))
print("Virtual Temperature ICAO (K) : ",round(TV_ICAO,2))

print("======================")
print("Zeroing Athmosphere ")
print("======================")

print("Altitude (m) :",int(Alt_Zero))
# values provided by Kestrel
print("Air pressure (Pa) : ",P_Zero)
print("Air Temperature (°C) : ",T_C_Zero,"(°K) : ",T_K_Zero)
print("Air Relative Humidity (%) : ",(RH_Zero*100))
print("Wet Air Volumic Mass (kg/m3) : ",round(Rho_Zero,3))
print("Saturated Vapor Pressure (Pa) : ",round(PVS_Zero,5))
print("Vapor Pressure (Pa) : ",round(PV_Zero,2))
print("Virtual Temperature (K) : ",round(TV_Zero,2))

print("======================")
print("Shooting Athmosphere ")
print("======================")

print("Altitude (m) :",int(Alt))
# values provided by Kestrel
print("Air pressure (Pa) : ",P)
print("Air Temperature (°C) : ",T_C,"(°K) : ",T_K)
print("Air Relative Humidity (%) : ",(RH*100))
print("Wet Air Volumic Mass (kg/m3) : ",round(Rho,3))
print("Saturated Vapor Pressure (Pa) : ",round(PVS,5))
print("Vapor Pressure (Pa) : ",round(PV,2))
print("Virtual Temperature (K) : ",round(TV,2))

print("======================")
print("Shooting Conditions ")
print("======================")

print("Shooting Distance : ",int(D_Tir))
print("Time increment (s) : ",h)
print("Shooting Angle relative to North (Azimut °) : ",Az_d," (RAD) : ",Az)

print("======================")
print("Coriolis Data ")
print("======================")

print("Earth Angular Speed - Omega (rad/s) :",Omega)

print("======================")
print("ICAO Drag Coefficient (Cd) Determination")
print("======================")

print("Speed of sound ICAO (m/s) :",round(SoundSpeed_ICAO,2))
print("Bullet Stability Factor ICAO ",round(Sg_ICAO,2))

if Sg_ICAO<1:
    print("ICAO Bullet is NOT stable")
if Sg_ICAO<1.5:
    print("ICAO Bullet Marginaly Stable")
if Sg_ICAO >=1.5:
    print("ICAO Stable Bullet")

print("======================")
print("Wind Speed and Direction ")
print("======================")

print("wind speed (m/s) :",Ws)
print("wind Angle relative to shooting direction (hour) :",int(Wa_o))

print("======================")
print(" Calculus Options ")
print("======================")
print("Corriolis : ",Co)

if Co == "N":
    print("Corriolis = N Effect NOT taken into account")
else:
    print("Corriolis = Y Effect taken into account")


if Zero_C == "Y":
    print("========== ZEROING ============")
    print("Zeroing at : ",D_zero," (m)")
    alpha_d=Zeroing(D_zero,Prec)
    print("======== ZEROING OVER =========")
else:
    print("======================")
    print("No Zeroing")
    print("======================")
print("Alpha(0) (deg) : ",round(alpha_d,6))

print("=======================================")
print("Ruge-Kutta Numerical Method Runing ....")
print("=======================================")

RK(D_Tir,h,MuzzleSpeed,NS,Az_d,Az,Co,Ws,Wa_o,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle,d_i,bl_cal,alpha_d,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO)

tph=float(tph-h)

Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))

SD=SpinDrift(Sg,tph)

WDx=float(0)
WDy=float(0)
WDy=float(Wsx*(tph+(round(D_Tir,0)/V0)))
WDz=float(Wsz*(tph+(round(D_Tir,0)/V0)))

print("===================================")
print("Ruge-Kutta Numerical Method Results")
print("===================================")
print("Lattitude (° N/S) :",round(Lat_D,2),NS)
print("Shooting Direction (Azimut Angle) : ",Az_d," °",Az,"RAD")
if Co == "N":
    print("Corriolis parameter = ",Co,"Effects NOT taken into account")
else:
    print("Corriolis parameter = ", Co, " Effects taken into account")

if SD_p == "N":
    print("Spin Drift parameter = ",SD_p,"Effects NOT taken into account")
else:
    print("Spin Drift parameter = ",SD_p,"Effects taken into account")
print("Goal Distance (m) :",round(D_Tir,0))
print("wind speed (m/s) :",Ws)
print("wind Angle relative to shooting direction (hour) :",Wa_o)
print("Wind Drift Along Y (m) :",WDy)
print("Wind Drift Along Z (m) :",WDz)
print("Time increment (s) : ",h," (ms) :",(h*1000))
print("Time of Flight (s) :",round(tph,4))
print("Bullet Stability Factor",round(Sg,2))
print("Spin Drift (m) : ",round((SD/100),3),"(cm): ",round(SD,2))
print("Calculated Speed Module |V(t+h)| (m/s) :",round(Vmod,2))
print("Calculated X position (m) :",round(Xtph,2))
print("Calculated Y position (m) :",round(Ytph,2)," (cm) :",round(Ytph*100,2))
Y_a=float((Ytph*1000)/D_Tir)
print("Calculated Y Angle (mRAD) :",round(Y_a,2))
E_c=float(-Y_a/0.1)

print("Calculated Y position includind Wind Drift (m) :",round((WDy+Ytph),2)," (cm) :",round((WDy+Ytph)*100,2))

Y_a=float(((WDy+Ytph)*1000)/D_Tir)
print("Calculated Y Angle includind Wind Drift (mRAD) :",round(Y_a,2))
E_c=float(-Y_a/0.1)

print("Calculated Z position (m) :",round(Ztph,2))
Z_a=float((Ztph*1000)/D_Tir)
print("Calculated Z Angle (mRAD) :",round(Z_a,2))

if SD_p == "N":
    SD=float(0)
print("Calculated Z position including Spin Drift (m) :",round(((SD/100)+Ztph),2))
Z_a=float(((SD/100)+Ztph)*1000/D_Tir)
print("Calculated Z Angle including Spin Drift (mRAD) :",round(Z_a,2))

print("Calculated Z position including Spin Drift and Wind Drift (m) :",round((WDz+(SD/100)+Ztph),2)," (cm) :",round((WDz*100+SD+Ztph*100),2))
Z_a=float((WDz+(SD/100)+Ztph)*1000/D_Tir)
print("Calculated Z Angle including Spin Drift and Wind Drift (mRAD) :",round(Z_a,2))
W_c=float(-Z_a/0.1)

print("========== CORRECTIONS TO BE APPLIED ================")

print("Elevation to be applied (clicks) +=>Up -=>Down:",round(E_c,1))

print("Windage to be applied (clicks) +=>Rigt -=>Left:",round(W_c,1))
print("=====================================================")
