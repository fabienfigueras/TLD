#!/usr/bin/python3

import sys
import math
import time
from datetime import date


# Class Definition

class Montage:
    def __init__(self,F_SightHeight=float(0.070),F_FixedAngle_d=float(0)):
        self.SightHeight=F_SightHeight
        self.FixedAngle_d=F_FixedAngle_d
    def show(self):
        print("SightHeight (mm): ",self.SightHeight*1000,"Fixed Angle (mRAD):",round(self.FixedAngle_d,2))

#FNlunette= "./CSV/lunette.csv"

class Rifle:
    def __init__(self,F_Brand="Rifle Brand",F_Model="Rifle Model",F_d_i=float(220),F_Twist=float(9),F_Twist_C="Left"):
        self.Rifle_Brand=F_Brand
        self.Rifle_Model=F_Model
        self.d_i=F_d_i
        self.Twist=F_Twist
        self.Twist_C=F_Twist_C
        self.T_cal=float(self.Twist/self.d_i)
    def show(self):
        print("Rifle Brand : ",self.Rifle_Brand,"Rifle Model :",self.Rifle_Model,"Rifle Caliber (inch) :",self.d_i,"Barrel Twist (inch) 1:",int(self.Twist),"Barrel Twist (R/L) :",self.Twist_C,"Barrel Twist in Caliber",self.T_cal)
    def update(self):
        self.T_cal=float(self.Twist/self.d_i)

class Amo:
    def __init__(self,F_Amo_Brand="Home Made Amo",F_Amo_Model="Home Made Amo Model",F_BC_G1_a=float(-1),F_BC_G7_a=float(-7),F_MS=float(800)):
        self.Amo_Brand=F_Amo_Brand
        self.Amo_Model=F_Amo_Model
        self.BC_G1_a=F_BC_G1_a
        self.BC_G7_a=F_BC_G7_a
        self.MuzzleSpeed=F_MS
    def show(self):
        print("Amo Brand : ",self.Amo_Brand,"Amo Model :",self.Amo_Model," BC G1 :",self.BC_G1_a," BC G7 : ",self.BC_G7_a," Muzzle Speed (m/s) ",self.MuzzleSpeed)

class Bullet:
    def __init__(self,F_Bullet_Brand="Home Made Bullet",F_Bullet_Model="Home Made Bullet Model",F_d_i=float(0.308),F_m_gr=float(190),F_bl_cm=float(3.437),F_MuzzleSpeed=float(780),F_BC_G1=float(-1),F_BC_G7=float(-7)):
        self.Bullet_Brand=F_Bullet_Brand
        self.Bullet_Model=F_Bullet_Model
        self.d_i=F_d_i
        self.m_gr=F_m_gr
        self.bl_cm=F_bl_cm
        self.MuzzleSpeed=F_MuzzleSpeed
        self.BC_G1=float(F_BC_G1)
        self.BC_G7=float(F_BC_G7)
        self.bl_i=float(self.bl_cm/2.54)
        self.bl_cal=float(self.bl_i/self.d_i)
        self.V0=float(self.MuzzleSpeed)
    def show(self):
        print("Bullet Brand : ",self.Bullet_Brand," Bullet Model :",self.Bullet_Model,"Bullet Diameter (inch) : ",self.d_i," Bullet Mass (gr) :",self.m_gr," Bullet Length (cm) :",self.bl_cm," Muzzle Speed (m/s) :",self.MuzzleSpeed,"BC_G1 :",self.BC_G1,"BC_G7 :",self.BC_G7,)
    def update(self):
        self.bl_i=float(self.bl_cm/2.54)
        self.bl_cal=float(self.bl_i/self.d_i)

class Bullet_BC:
    def __init__(self,F_Vmin_b=0,F_BC_min_b=0,F_BC_int_b=0,F_Vmax_b=0,F_BC_max_b=0):
        self.Vmin_b=F_Vmin_b
        self.BC_min_b=F_BC_min_b
        self.BC_int_b=F_BC_int_b
        self.Vmax_b=F_Vmax_b
        self.BC_max_b=F_BC_max_b
    def show(self):
        print(" BC = ", self.BC_min_b, " if speed (m/s) <= ", self.Vmin_b)
        print(" BC = ", self.BC_int_b,   " if speed is > (m/s) ",self.Vmin_b, " and <= ",self.Vmax_b)
        print(" BC = ", self.BC_max_b,   " if speed is > (m/s) ",self.Vmax_b)

#FNzero= "./CSV/zero.csv"

class Zero:
    def __init__(self,F_D_zero=100,F_Prec=0.001,F_alpha_d=0,F_Wind_cm=0):
        self.D_zero = float(F_D_zero)
        self.Prec = float(F_Prec)
        self.atm = self.ZAtm()
        self.alpha_d = float(F_alpha_d)
        self.wind_cm = float(F_Wind_cm)
    def show(self):
        print("Zero Distance (m)",self.D_zero, " Error tolerance (m) ", self.Prec)
        print("Zero Atmosphere Data ")
        self.atm.show()
        print("Zero Angle (deg) ",self.alpha_d,"Windage (cm) - =Left + =Right ",self.wind_cm)

    class ZAtm:
        def __init__(self,F_Alt=0,F_P=101325,F_T_C=15,F_RH=0):
            self.Alt=float(F_Alt)
            self.P=float(F_P)
            self.T_C=float(F_T_C)
            self.T_K=float(self.T_C+273.15)
            self.RH=float(F_RH)
            self.Rho=float(0.003483*((1 + self.RH) / (1 + 1.6078 * self.RH))* ( self.P/ (self.T_C + 273.15) ))
            self.PVS=float(6.1078*math.pow(10,((7.5*self.T_K-2048.625)/(self.T_K-35.85))))
            self.PV=float(self.RH*self.PVS)
            self.TV=float((self.T_K/(1-0.3785*(self.PV/self.P))))
        def show(self):
            print("Altitude (m)", self.Alt, "Absolute Pressure (Pa) ", self.P)
            print("Air Temperature (°C)", self.T_C, "Air Temperature (°K)", self.T_K)
            print("Air Relative Humidity (%) ",self.RH*100)
            print("Wet Air Volumic Mass (kg/m3)", round(self.Rho,3))
            print("Saturated Vapor Pressure (Pa) : ",round(self.PVS,5))
            print("Vapor Pressure (Pa) : ",round(self.PV,2))
            print("Virtual Temperature (K) : ",round(self.TV,2))
        def update(self):
            self.T_K=float(self.T_C+273.15)
            self.Rho=float(0.003483*((1 + self.RH) / (1 + 1.6078 * self.RH))* ( self.P/ (self.T_C + 273.15) ))
            self.PVS=float(6.1078*math.pow(10,((7.5*self.T_K-2048.625)/(self.T_K-35.85))))
            self.PV=float(self.RH*self.PVS)
            self.TV=float((self.T_K/(1-0.3785*(self.PV/self.P))))

class Atm:
    def __init__(self,F_Alt=0,F_P=101325,F_T_C=15,F_RH=0.0):
        self.Alt=float(F_Alt)
        self.P=float(F_P)
        self.T_C=float(F_T_C)
        self.T_K=float(self.T_C+273.15)
        self.RH=float(F_RH)
        self.Rho=float(0.003483*((1 + self.RH) / (1 + 1.6078 * self.RH))* ( self.P/ (self.T_C + 273.15) ))
        self.PVS=float(6.1078*math.pow(10,((7.5*self.T_K-2048.625)/(self.T_K-35.85))))
        self.PV=float(self.RH*self.PVS)
        self.TV=float((self.T_K/(1-0.3785*(self.PV/self.P))))
    def show(self):
        print("Altitude (m)", self.Alt, "Absolute Pressure (Pa) ", self.P)
        print("Air Temperature (°C)", self.T_C, "Air Temperature (°K)", self.T_K)
        print("Air Relative Humidity (%) ",self.RH*100)
        print("Wet Air Volumic Mass (kg/m3)", round(self.Rho,3))
        print("Saturated Vapor Pressure (Pa) : ",round(self.PVS,5))
        print("Vapor Pressure (Pa) : ",round(self.PV,2))
        print("Virtual Temperature (K) : ",round(self.TV,2))
    def update(self):
        self.T_K=float(self.T_C+273.15)
        self.Rho=float(0.003483*((1 + self.RH) / (1 + 1.6078 * self.RH))* ( self.P/ (self.T_C + 273.15) ))
        self.PVS=float(6.1078*math.pow(10,((7.5*self.T_K-2048.625)/(self.T_K-35.85))))
        self.PV=float(self.RH*self.PVS)
        self.TV=float((self.T_K/(1-0.3785*(self.PV/self.P))))

#FNenv= "./CSV/env.csv"

class Wind:
    def __init__(self,F_Ws=0,F_Wa_o=3):
        self.Ws=float(F_Ws)
        self.Wa_o=float(F_Wa_o)
        self.Wa=float((self.Wa_o*(math.pi/2))/3)
        self.Wsx=-self.Ws*math.cos(self.Wa)
        self.Wsy=-float(0)
        self.Wsx=-self.Ws*math.sin(self.Wa)
        self.WDx=float(0)
        self.WDy=float(0)
        self.WDz=float(0)
    def show(self):
        print("Average Wind intensity (m/s)", self.Ws,"Heading from (hour)",self.Wa_o," related to shooting direction")
        print("Heading Angle in RAD :",round(self.Wa,5))
        print("Resulting Wind Speed on X axis (m/s)", round(self.Wsx,2))
        print("Resulting Wind Speed on Y axis (m/s)", round(self.Wsy,2))
        print("Resulting Wind Speed on Z axis (m/s)", round(self.Wsz,2))
        print("Resulting Deviation on X direction (m)", round(self.WDx,3))
        print("Resulting Deviation on Y direction (m)", round(self.WDy,3))
        print("Resulting Deviation on Z direction (m)", round(self.WDz,3))
    def update(self):
        self.Wa=float((self.Wa_o*(math.pi/2))/3)
        self.Wsx=-self.Ws*math.cos(self.Wa)
        self.Wsy=-float(0)
        self.Wsz=-self.Ws*math.sin(self.Wa)

class Impact:
    def __init__(self,iX=1,iY=0,iZ=0,iVx=0,iVy=0,iVz=0):
        self.X=float(iX)
        self.Y=float(iY)
        self.Z=float(iZ)
        self.Vx=float(iVx)
        self.Vy=float(iVy)
        self.Vz=float(iVz)
        self.V=float(math.sqrt(self.Vx*self.Vx+self.Vy*self.Vy+self.Vz*self.Vz))
        self.E_a=float(0)
        self.W_a=float(0)
        self.Vxz=float(self.V*math.cos(self.W_a))
    def show(self):
        print("Impact point Data ")
        print("X coordinate (distance from shooting point) (m)", round(self.X,4))
        print("Y coordinate (Lateral drift from shooting point) (m)", round(self.Y,4))
        print("Z coordinate (Vertical Drop from shooting point) (m)", round(self.Z,4))
        print("Y coordinate (Lateral drift from shooting point) (cm)", round(self.Y*100,1))
        print("Z coordinate (Vertical Drop from shooting point) (cm)", round(self.Z*100,1))
        print("Speed coordinate on X axis (m/s)", round(self.Vx,5))
        print("Speed coordinate on Y axis (m/s)", round(self.Vy,5))
        print("Speed coordinate on Z axis (m/s)", round(self.Vz,5))
        print("Speed Module (m/s)", round(self.V,5))
        print("Speed Module on XY plan (m/s)", round(self.Vxz,5))
        print("Elevation Angle (RAD)", round(self.E_a,5))
        print("Windage Angle (RAD)", round(self.W_a,5))
    def update(self):
        self.V=float(math.sqrt((self.Vx*self.Vx)+(self.Vy*self.Vy)+(self.Vz*self.Vz)))
        self.Vxy=float(self.V*math.cos(self.W_a))
        self.Vz=float(self.V*math.sin(self.W_a))
        self.Vx=float(self.Vxy*math.cos(self.E_a))
        self.Vy=float(self.Vxy*math.sin(self.E_a))

# Software introduction for users

PBS_version="1.28"
PBS_year="2025"

print("==================================================")
print(" PBS stands for Python Ballistic Solver")
print(" PBS is an Open Source Balistic Software ")
print(" Written in Object Oriented Python3 by Fabien FIGUERAS (he/him)")
print(" v1.00 was released in 2024 ")
print(" Current Version is v",PBS_version,PBS_year)
print(" Call example Python3 ./PBS-vxyz.py to get this message")
print("")
print(" Next 3 Parameters will be overwritten by Files values ")
print(" Where param1 is the caliber [inch] ")
print(" Where param2 is the bullet mass [gr] ")
print(" Where param3 is the Muzzle Speed [m/s] ")
print(" Where param4 is the Shooting distance [m]")
print(" Where param5 is the Azimut (shooting angle relative to the North) [deg]")
print(" Where param6 is the Coriolis Option [Y/N] ")
print(" Where param7 is the Average Wind Speed [m/s]")
print(" Where param8 is the Wind Speed direction related to shooting direction [hour] ")
print(" Where param9 is the Spind Drift Option [Y/N] ")
print(" Where param10 is the time increment for numerical solution [s] ")
print(" Where param11 is Zeroing the sight ? [Y/N] ")
print(" Where param12 is the Shooting Angle (relative to the Horizontal plan) required for Coriolis option [deg] ")
print(" Where param13 is Aerodynamic Jump Option ? [Y/N] ")
print(" Where param14 is BC_Gx type ? [G1/G7]")
print(" Next Option could force BC_Gx to be overwritten by Files values ")
print(" Where param15 is BC_Gx value ? [0 constant, 1 Speed related]")
print(" Where param16 is the option to allow calculation of Card or Abacus or Nothing  [C/Y/N] ")
print("")
print(" Sources available in GitHub : https://github.com/fabienfigueras/TLD")
print("==================================================")

# Objects creation

todays_date = date.today()
montage = Montage()
rifle = Rifle()
amo = Amo()
bullet = Bullet()
bullet_BC = Bullet_BC()
zero = Zero()
atm_zero = zero.atm
ICAO_Atm = Atm()
Shoot_Atm = Atm()
Shoot_Wind = Wind()
impact = Impact()
impact_HA = Impact()
impact_HA_Co = Impact()
impact_Co = Impact()
impact_NoHa_NoCo = Impact()
impact_NoC = Impact()

# Objects and variables First initialization

Wind_Range=0

bullet.d_i=float(sys.argv[1])
bullet.m_gr=float(sys.argv[2])
bullet.bl_cm=float(3.437)
bullet.bl_i=float(bullet.bl_cm/2.54)
bullet.bl_cal=float(bullet.bl_i/bullet.d_i)
bullet.MuzzleSpeed=float(sys.argv[3])

#Shooting distance (m))
D_Tir=float(sys.argv[4])
D_Tir_c=float(D_Tir)
#Shooting angle relative to the North [deg]
Az_d=float(sys.argv[5])
Az=float((Az_d*math.pi/180))
#Coriolis Option [Y/N]
Co=sys.argv[6]
#Average Wind Speed [m/s]
Shoot_Wind.Ws=float(sys.argv[7])
#Wind Speed direction related to shooting direction [hour]
Shoot_Wind.Wa_o=float(sys.argv[8])
#Spind Drift Option [Y/N]
SD_p=sys.argv[9]
#Time increment for numerical solution [s]
h=float(sys.argv[10])
#Zeroing the sight ? [Y/N]
Zero_C=sys.argv[11]
#Shooting Angle (relative to the Horizontal plan) required for Coriolis option [deg]
H_Angle_d=float(sys.argv[12])
H_Angle=float((H_Angle_d*math.pi/180))
#Shooting distance corrected for Horizontal Angle (m))
D_Tir_HA=float(D_Tir*math.cos(H_Angle))
Ytph_HA=float(0)
Ytph_HA_Co=float(0)
Ytph_Co=float(0)
Ytph_NoHa_NoCo=float(0)
Y_Drop=float(0)
Y_a_NoHa_NoCo=float(0)
Y_a_NoHa_Co=float(0)
Y_a_NoC=float(0)
Y_a_Ha_Co=float(0)
ToF=float(0)

#Aerodynamic Jump Option ? [Y/N]
AJ=sys.argv[13]
#BC_Gx type ? [G1/G7]
BC_Gx=sys.argv[14]
#BC_Gx value ? [0 constant, 1 Speed related]
BC_Go=float(sys.argv[15])

#Variables used for numerical resolution
#Time
tph=float(0)
#Speed along Axes X,Y,Z depending on time
Vxtph=float(0)
Vytph=float(0)
Vztph=float(0)
#Spin Drift
SD=float(0)
#Dynamic Stability
Sg=float(0)
#Aerodynamic Jump variables
AJ_MOA=float(0)
AJ_mRAD=float(0)
AJD_z=float(0)
AJE=float(0)

#Wind Speed along Axis X,Y,Z (m/s) done when object was created
#Wind Drift along Axis X,Y,Z (m)) done when object was created
#update for others variables
Shoot_Wind.update()

#initial bullet speed
bullet.V0=float(bullet.MuzzleSpeed)
#Latitude
Lat_d=float(46)
Lat_m=float(22)
Lat_s=float(25)
Lat_D=float(Lat_d+(Lat_m/60)+(Lat_s/(60*60)))
NS="North"

#Average Speed
Vmod=float(0)
#Bullet position depending on time
Xtph=float(0)
Ytph=float(0)
Ztph=float(0)

#Elevation correction in click
E_c=float(0)
#Windage correction in click
W_c=float(0)

#Variables used to calculate the maximum height of the bullet
Y_max=float(0)
XY_max=float(0)
TY_max=float(0)
Arrow=0

#Physical constants
Gamma=float((1+(1/2.48)))
Omega=float(2*math.pi/((23+(56/60)+(4.09/(60*60)))*(60*60)))

#Settings ICAO Standard Atmosphere values
ICAO_Atm.Alt=float(0)
ICAO_Atm.P=float(101325)
ICAO_Atm.T_C=float(15)
ICAO_Atm.RH=float(0)
ICAO_Atm.update()

SoundSpeed_ICAO=float(math.sqrt(Gamma*ICAO_Atm.P/ICAO_Atm.Rho))

Sg_ICAO=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((ICAO_Atm.T_K*101325)/(288.15*ICAO_Atm.P)))))

#Default values, will be overwritten by those from files
Shoot_Atm.Alt=float(440)
Shoot_Atm.P=float(96330)
Shoot_Atm.T_C=float(26)
Shoot_Atm.RH=0.66
Shoot_Atm.update()

c_mRAD=float(0.05)

print("===============================================")
print("Gathering and printing Data from Files ")
print("File parameters overcome some Command line parameters")
print("===============================================")

# EXTERNAL DATA

FNmontage= "./CSV/montage.csv"
FNlunette= "./CSV/lunette.csv"
FNrifle= "./CSV/rifle.csv"
FNamo= "./CSV/amo.csv"
FNbullet= "./CSV/bullet.csv"
FNbullet_BC= "./CSV/bullet_BC.csv"
FNzero= "./CSV/zero.csv"
FNenv= "./CSV/env.csv"
FNabacus= "./Abacus.csv"

# CHARGEMENT DES PARAMETRES SELON FICHIERS

def InitFromFile():
    global FNmontage,montage,FNlunette,FNrifle,rifle,FNamo,amo,FNbullet,bullet,FNbullet_BC,bullet_BC,FNzero,zero,atm_zero,FNenv,Lat_d,Lat_m,Lat_s,Lat_D,NS,Shoot_Atm,Sg_ICAO,c_mRAD

    fh = open(FNmontage,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    montage.SightHeight=float(DL[0])/1000
    montage.FixedAngle_d=float(DL[1])
#    montage.show()
    fh.close()

    fh = open(FNlunette,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    c_mRAD=float(DL[1])
    fh.close()

    fh = open(FNrifle,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    rifle.Rifle_Brand=DL[0]
    rifle.Rifle_Model=DL[1]
    rifle.d_i=float(DL[2])
    rifle.Twist=float(DL[3])
    rifle.Twist_C=DL[4]
    rifle.T_cal=float(rifle.Twist/rifle.d_i)
#    rifle.show()
    fh.close()

    fh = open(FNamo,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    amo.Amo_Brand=DL[0]
    amo.Amo_Model=DL[1]
    amo.BC_G1_a=float(DL[2])
    amo.BC_G7_a=float(DL[3])
    amo.MuzzleSpeed=float(DL[4])
#    amo.show()
    fh.close()

    fh = open(FNbullet,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    bullet.Bullet_Brand=DL[0]
    bullet.Bullet_Model=DL[1]
    bullet.d_i=float(DL[2])
    bullet.m_gr=float(DL[3])
    bullet.bl_cm=float(DL[4])
    bullet.MuzzleSpeed=float(DL[5])
    bullet.BC_G1=float(DL[6])
    bullet.BC_G7=float(DL[7])
    bullet.bl_i=float(bullet.bl_cm/2.54)
    bullet.bl_cal=float(bullet.bl_i/bullet.d_i)
    bullet.V0=float(bullet.MuzzleSpeed)
#    bullet.show()
    fh.close()

    fh = open(FNbullet_BC,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    bullet_BC.Vmin_b=float(DL[0])
    bullet_BC.BC_min_b=float(DL[1])
    bullet_BC.BC_int_b=float(DL[2])
    bullet_BC.Vmax_b=float(DL[3])
    bullet_BC.BC_max_b=float(DL[4])
#    bullet_BC.show()
    fh.close()

    fh = open(FNzero,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    zero.D_zero=float(DL[0])
    zero.Prec=float(DL[1])
    atm_zero.Alt=float(DL[2])
    atm_zero.P=float(DL[3])
    atm_zero.T_C=float(DL[4])
    atm_zero.RH=float(DL[5])
    zero.alpha_d=float(DL[6])
    zero.wind_cm=float(DL[7])
    atm_zero.update()
#    zero.show()
#    atm_zero.show()
    fh.close()

    fh = open(FNenv,"r")
    header=fh.readline()
    data=fh.readline()
    data=data.replace("\n","")
    DL=data.split(",")
    #Ws=float(DL[4])
    #Wa_o=float(DL[5])
    Lat_d=float(DL[6])
    Lat_m=float(DL[7])
    Lat_s=float(DL[8])
    Lat_D=float(Lat_d+(Lat_m/60)+(Lat_s/(60*60)))
    NS=DL[9]
    Shoot_Atm.Alt=float(DL[0])
    Shoot_Atm.P=float(DL[1])
    Shoot_Atm.T_C=float(DL[2])
    Shoot_Atm.RH=float(DL[3])
    Shoot_Atm.update()
#    Shoot_Atm.show()

    fh.close()

    Sg_ICAO=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((ICAO_Atm.T_K*101325)/(288.15*ICAO_Atm.P)))))

# Definitions of all functions

def ResetData():
    global bullet,D_Tir,Az_d,Co,Ws,Wa_o,SD_p,h,Zero_C,H_Angle_d,H_Angle,AJ,BC_Gx,BC_Go,E_c,W_c,tph,Vxtph,Vytph,Vztph,SD,Sg,AJ_MOA,AJ_mRAD,AJD_z,AJE,WDx,WDy,WDz,Wsx,Wsy,Wsz,Lat_d,Lat_m,Lat_s,Lat_D,NS,Az,Vmod,Xtph,Ytph,Ztph,Gamma,D_Tir_HA,Ytph_HA,Y_max,XY_max,TY_max,Arrow,ICAO_Atm,SoundSpeed_ICAO,Shoot_Atm,Omega,Shoot_Wind,Wind_Range

    bullet.d_i=float(sys.argv[1])
    bullet.m_gr=float(sys.argv[2])
    bullet.bl_cm=float(3.437)
    bullet.bl_i=float(bullet.bl_cm/2.54)
    bullet.bl_cal=float(bullet.bl_i/bullet.d_i)
    bullet.MuzzleSpeed=float(sys.argv[3])

    D_Tir=float(sys.argv[4])
    D_Tir_c=float(D_Tir)
    Az_d=float(sys.argv[5])
    Co=sys.argv[6]
    Shoot_Wind.Ws=float(sys.argv[7])
    Shoot_Wind.Wa_o=float(sys.argv[8])
    SD_p=sys.argv[9]
    h=float(sys.argv[10])
    Zero_C=sys.argv[11]
    H_Angle_d=float(sys.argv[12])
    H_Angle=float((H_Angle_d*math.pi/180))
    AJ=sys.argv[13]
    BC_Gx=sys.argv[14]
    BC_Go=float(sys.argv[15])

    E_c=float(0)
    W_c=float(0)

    tph=float(0)
    Vxtph=float(0)
    Vytph=float(0)
    Vztph=float(0)
    SD=float(0)
    Sg=float(0)
    AJ_MOA=float(0)
    AJ_mRAD=float(0)
    AJD_z=float(0)
    AJE=float(0)

    Shoot_Wind.WDx=float(0)
    Shoot_Wind.WDy=float(0)
    Shoot_Wind.WDz=float(0)
    Shoot_Wind.update()

    bullet.V0=float(bullet.MuzzleSpeed)
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
    D_Tir_HA=float(0)
    Ytph_HA=float(0)
    Ytph_HA_Co=float(0)
    Ytph_Co=float(0)
    Ytph_NoHa_NoCo=float(0)
    Y_Drop=float(0)
    Y_a_NoHa_NoCo=float(0)
    Y_a_NoHa_Co=float(0)
    Y_a_NoC=float(0)
    Y_a_Ha_Co=float(0)
    ToF=float(0)

    Y_max=float(0)
    XY_max=float(0)
    TY_max=float(0)
    Arrow=0

    InitFromFile()

# BC_Gx calculation
def BC_Search(V):
    global BC_Gx,BC_Go,amo,bullet_BC

    # amo.show()
    match BC_Gx:
        case "G1":
            match BC_Go:
                case 0:
                    BC_G=amo.BC_G1_a
                case 1:
                    if V>=bullet_BC.Vmin_b:
                        BC_G=bullet_BC.BC_int_b
                    if V>=bullet_BC.Vmax_b:
                        BC_G=bullet_BC.BC_max_b
                    if V<bullet_BC.Vmin_b:
                        BC_G=bullet_BC.BC_min_b
                case _:
                    BC_G=-1
        case "G7":
            match BC_Go:
                case 0:
                    BC_G=amo.BC_G7_a
                case 1:
                    if V>=bullet_BC.Vmin_b:
                        BC_G=bullet_BC.BC_int_b
                    if V>=bullet_BC.Vmax_b:
                        BC_G=bullet_BC.BC_max_b
                    if V<bullet_BC.Vmin_b:
                        BC_G=bullet_BC.BC_min_b
                case _:
                    BC_G=-7
        case _:
            BC_G=float(0)
    return(BC_G)

# Search Departure Angle to get Y=0 for a given distance

def Zeroing(D_zero,Prec):
    global h,bullet,NS,Az_d,Az,Co,montage,rifle,Shoot_Atm,Gamma,ICAO_Atm,Omega,SoundSpeed_ICAO,tph,Sg,zero,Impact,c_mRAD

    Delta_Y=float(1)
    cpt=1

#    print("Zero distance (m):",D_zero," Precision required (m): ",Prec)

    Alpha_d=float(0)

    while Delta_Y>Prec:
#        print("Zeroing Phase ", cpt)

        PRS_Solver(D_zero,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,Alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)

#        impact.show()

        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))

        tph=float(tph-h)

        SD=SpinDrift(Sg,tph)

        zero.wind_cm=SD

        Alpha_RAD_i=float(-math.atan(impact.Y/D_zero))

#        print("Passe ",cpt," : Ecart en Y (m) ", (impact.Y))
#        print("Passe ",cpt," : Corrections en Elevation pour le Zero (click) ", (Alpha_RAD_i*1000/c_mRAD))

        Alpha_d+=float(Alpha_RAD_i*180/math.pi)

        Delta_Y=abs(impact.Y)
        cpt+=1

    return(Alpha_d)

# Spin Drift calculation

def SpinDrift(Sg,ToF):
 Dg=float(3.175*(Sg+1.2)*math.pow(ToF,1.83))
 return(Dg)

# Retreive Cd_G1 for a given Mach number

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

# Retreive Cd_G7 for a given Mach number

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

# Retreive Drag Force on All axis

# Retreive Drag Force on X axis

def Fx(t,Vx,Vy,Vz,V,Ox,Oy,Oz,LCo,k,m,g):
    global Wind_Range

    if LCo == "N":
#
        Fx=float((-k/m)*abs(Vx)*Vx)
#        Fx=float((-k/m)*abs(V)*Vx)
    else:
#
        Fx=float((-k/m)*abs(Vx)*Vx-2*(Oy*Vz-Oz*Vy))
#        Fx=float((-k/m)*abs(V)*Vx-2*(Oy*Vz-Oz*Vy))
# G Klimi        Fx=float((-k/m)*abs(Vx)*Vx-2*(-Oz*Vy))
#        Fx=float((-k/m)*abs(V)*Vx-2*(-Oz*Vy))
    return(Fx)

# Retreive Drag Force on Y axis

def Fy(t,Vx,Vy,Vz,V,Ox,Oy,Oz,LCo,k,m,g):
    global Wind_Range

    if LCo == "N":
#
        Fy=float((-k/m)*abs(Vy)*Vy-g)
#        Fy=float((-k/m)*abs(V)*Vy-g)
    else:
#
        Fy=float((-k/m)*abs(Vy)*Vy-g-2*(-(Ox*Vz-Oz*Vx)))
#        Fy=float((-k/m)*abs(V)*Vy-g-2*(-(Ox*Vz-Oz*Vx)))
# G Klimi        Fy=float((-k/m)*abs(Vy)*Vy-g-2*(-(-Oz*Vx)))
#        Fy=float((-k/m)*abs(V)*Vy-g-2*(-(-Oz*Vx)))
    return(Fy)

# Retreive Drag Force on Z axis

def Fz(t,Vx,Vy,Vz,V,Ox,Oy,Oz,LCo,k,m,g):
    global Wind_Range

    if LCo == "N":
#
        Fz=float((-k/m)*abs(Vz)*Vz)
#        Fz=float((-k/m)*abs(V)*Vz)
    else:
#
        Fz=float((-k/m)*abs(Vz)*Vz-2*(Ox*Vy-Oy*Vx))
#        Fz=float((-k/m)*abs(V)*Vz-2*(Ox*Vy-Oy*Vx))

    return(Fz)

# Retreive Cd knowing mass, BC_G1 and Stability Factor

def CdSearch(m,BC_G1,Sd):
    Cd_G1=CdG1Search(m)
    Cd=float(Sd*(Cd_G1/BC_G1))
    return(Cd)

# Runge Kutta Numrical Method used to solve Differential Equation of Movement

def RK(D_Tir,h,MuzzleSpeed,NS,Az_d,Az,Co,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle_Brand,Rifle_Model,d_i,bl_cal,alpha_d,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO):

    global E_c,W_c,tph,Vxtph,Vytph,Vztph,Sg,Shoot_Wind,bullet,Lat_D,Vmod,Xtph,Ytph,Ztph,Arrow,Y_max,XY_max,TY_max,BC_Gx,BC_Go,Shoot_Wind,impact,D_Tir_c,amo

    #Convert from Imperial to Metric units
    #Bullet diameter
    d_mm=float(d_i*25.4)
    d=float(d_mm/1000)
    #Bullet mass https://fr.wikipedia.org/wiki/Grain_(unit%C3%A9)
    m_g=float(m_gr*0.06479891)
    m=float(m_g/1000)

    #Bullet surface
    A=float(math.pi*(d*d)/4)

    #Form factor calculation
    Sd=float((m_gr/7000)/(d_i*d_i))

    #Latitude from degre to RAD
    Lat=float(Lat_D*math.pi/180)

    #Acceleration of Gravity calculation https://en.wikipedia.org/wiki/Gravity_of_Earth

    g=float(9.780327*(1+0.0053024*math.sin(Lat)*math.sin(Lat)-0.0000058*math.sin(2*Lat)*math.sin(2*Lat)))

    #Shooting angle converted from degre to RAD
    alpha=float(alpha_d*(math.pi/180))

    #Setting variables for Coriolis effect
    Ox=float(Omega*math.cos(Lat)*math.cos(Az))
    Oy=float(Omega*math.sin(Lat))
    Oz=float(-Omega*math.cos(Lat)*math.sin(Az))
    ModOmega=math.sqrt(Ox*Ox+Oy*Oy+Oz*Oz)

    #Setting variables for Speed
    bullet.V0=MuzzleSpeed

    Vx=float(bullet.V0*math.cos(alpha))+Shoot_Wind.Wsx
    Vy=float(bullet.V0*math.sin(alpha))
    Vz=float(0)
    V=float(math.sqrt(Vx*Vx+Vy*Vy+Vz*Vz))

    #Mach number calculation in ICAO conditions
#    Mach_ICAO=float(V/SoundSpeed_ICAO)
    Mach_ICAO=float(amo.MuzzleSpeed/SoundSpeed_ICAO)

    #Cd_Gx calculation in ICAO conditions
    Cd_G1_ICAO=CdG1Search(Mach_ICAO)
    Cd_G7_ICAO=CdG7Search(Mach_ICAO)

    match BC_Gx:
        case "G1":
            BC_G1_ICAO=BC_Search(V)
            BC_G7_ICAO=float(BC_G1_ICAO/(Cd_G1_ICAO/Cd_G7_ICAO))
        case "G7":
            BC_G7_ICAO=BC_Search(V)
            BC_G1_ICAO=float(BC_G7_ICAO/(Cd_G7_ICAO/Cd_G1_ICAO))

    #Cd calculation in ICAO conditions
    Cd_ICAO=float(Sd*(Cd_G7_ICAO/BC_G7_ICAO))

    #k calculation in ICAO conditions
    k_ICAO=float(0.5*Rho_ICAO*Cd_ICAO*A)

    #SoundSpeed calculation in Shooting conditions
    SoundSpeed=float(math.sqrt(Gamma*P/Rho))

    #Impedance Coefficient calculation with reference to ICAO conditions
    J=float((P_ICAO/P)*math.sqrt(TV/TV_ICAO))

    #Mach number calculation in Shooting conditions
    Mach=float(V/SoundSpeed)

    #Cd_Gx calculation in Shooting conditions
    Cd_G1=CdG1Search(Mach)
    Cd_G7=CdG7Search(Mach)

    #BC_Gx Calculation using Impedance Coefficient (J)
    BC_G1=float(BC_G1_ICAO*J)
    BC_G7=float(BC_G7_ICAO*J)

    #Cd calculation in shooting conditions
    Cd=float(Sd*(Cd_G7/BC_G7))

    #k calculation in Shooting conditions
    k=float(0.5*Rho*Cd*A)

    #initiating variables
    t=float(0)
    tph=float(t)
    X=float(0)
    Y=float(-SightHeight)
    Z=float(0)
    tph=float(tph+h)
    Y_max=float(0)
    XY_max=float(0)
    TY_max=float(0)

    #initiating impact
    impact.X=X
    impact.Y=Y
    impact.Z=Z
    impact.Vx=Vx
    impact.Vy=Vy
    impact.Vz=Vz
    impact.E_a=alpha
    impact.W_a=float(0)
    impact.update()
#    impact.show()

    Cd_org=Cd
    BC_G1_B=BC_G1
    BC_G7_B=BC_G7

    #Calculating Runge-Kutta coefficients

    Kx1=float(Fx(t,Vx,Vy,Vz,V,Ox,Oy,Oz,Co,k,m,g))
    Ky1=float(Fy(t,Vx,Vy,Vz,V,Ox,Oy,Oz,Co,k,m,g))
    Kz1=float(Fz(t,Vx,Vy,Vz,V,Ox,Oy,Oz,Co,k,m,g))

    Kx2=float(Fx(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,V,Ox,Oy,Oz,Co,k,m,g))
    Ky2=float(Fy(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,V,Ox,Oy,Oz,Co,k,m,g))
    Kz2=float(Fz(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,V,Ox,Oy,Oz,Co,k,m,g))

    Kx3=float(Fx(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,V,Ox,Oy,Oz,Co,k,m,g))
    Ky3=float(Fy(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,V,Ox,Oy,Oz,Co,k,m,g))
    Kz3=float(Fz(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,V,Ox,Oy,Oz,Co,k,m,g))

    Kx4=float(Fx(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,V,Ox,Oy,Oz,Co,k,m,g))
    Ky4=float(Fy(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,V,Ox,Oy,Oz,Co,k,m,g))
    Kz4=float(Fz(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,V,Ox,Oy,Oz,Co,k,m,g))

    #Calculating Next step Speed
    Vxtph=float(Vx+(h/6)*(Kx1+2*Kx2+2*Kx3+Kx4))
    Vytph=float(Vy+(h/6)*(Ky1+2*Ky2+2*Ky3+Ky4))
    Vztph=float(Vz+(h/6)*(Kz1+2*Kz2+2*Kz3+Kz4))

    #Caclulating average speed between t and t+h
    Vxm=float((Vx+Vxtph)/2)
    Vym=float((Vy+Vytph)/2)
    Vzm=float((Vz+Vztph)/2)

    #Calculating the new bullet's coordinates
    X=float(X+h*Vxm)
    Y=float(Y+h*Vym)
    Z=float(Z+h*Vzm)

    #repeat the same process until the requested shooting distance is reached
    while (X<D_Tir) :
        #prepare for next step
        t=float(tph)

        Xtph=float(X+h*Vxm)
        Ytph=float(Y+h*Vym)
        Ztph=float(Z+h*Vzm)

        #If requested calculated the highest altitude of the bullet during its flight
        if Arrow!=0:
            if Ytph>=Y_max:
                Y_max=Ytph
                XY_max=Xtph
                TY_max=tph

        V=float(math.sqrt((Vx*Vx)+(Vy*Vy)+(Vz*Vz)))

        #update Mach,BC_Gx,Cd,k
        match BC_Gx:
            case "G1":
                BC_G1_ICAO=BC_Search(V)
                BC_G7_ICAO=float(BC_G1_ICAO/(Cd_G1_ICAO/Cd_G7_ICAO))
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
            case "G7":
                BC_G7_ICAO=BC_Search(V)
                BC_G1_ICAO=float(BC_G7_ICAO/(Cd_G7_ICAO/Cd_G1_ICAO))
                BC_G1=float(BC_G1_ICAO*J)
                if BC_G7_B!=BC_G7:

                    BC_G7_B=BC_G7

                    Mach=float(V/SoundSpeed)
                    Cd_G1=CdG1Search(Mach)
                    Cd_G7=CdG7Search(Mach)

                    BC_G1=float(BC_G7/(Cd_G1/Cd_G7))
                    Cd=float(Sd*(Cd_G7/BC_G7))

                    if Cd_org != Cd :
                        Cd_org=Cd

                    k=float(0.5*Rho*Cd*A)

        BC_G1=float(BC_G1_ICAO*J)

        #Calculating Runge-Kutta coefficients
        Kx1=float(Fx(t,Vx,Vy,Vz,V,Ox,Oy,Oz,Co,k,m,g))
        Ky1=float(Fy(t,Vx,Vy,Vz,V,Ox,Oy,Oz,Co,k,m,g))
        Kz1=float(Fz(t,Vx,Vy,Vz,V,Ox,Oy,Oz,Co,k,m,g))

        Kx2=float(Fx(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,V,Ox,Oy,Oz,Co,k,m,g))
        Ky2=float(Fy(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,V,Ox,Oy,Oz,Co,k,m,g))
        Kz2=float(Fz(t+(h/2),Vx+(h/2)*Kx1,Vy+(h/2)*Ky1,Vz+(h/2)*Kz1,V,Ox,Oy,Oz,Co,k,m,g))

        Kx3=float(Fx(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,V,Ox,Oy,Oz,Co,k,m,g))
        Ky3=float(Fy(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,V,Ox,Oy,Oz,Co,k,m,g))
        Kz3=float(Fz(t+(h/2),Vx+(h/2)*Kx2,Vy+(h/2)*Ky2,Vz+(h/2)*Kz2,V,Ox,Oy,Oz,Co,k,m,g))

        Kx4=float(Fx(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,V,Ox,Oy,Oz,Co,k,m,g))
        Ky4=float(Fy(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,V,Ox,Oy,Oz,Co,k,m,g))
        Kz4=float(Fz(t+(h/2),Vx+h*Kx3,Vy+h*Ky3,Vz+h*Kz3,V,Ox,Oy,Oz,Co,k,m,g))

        #Calculating Next step Speed
        Vxtph=float(Vx+(h/6)*(Kx1+2*Kx2+2*Kx3+Kx4))
        Vytph=float(Vy+(h/6)*(Ky1+2*Ky2+2*Ky3+Ky4))
        Vztph=float(Vz+(h/6)*(Kz1+2*Kz2+2*Kz3+Kz4))

        #Caclulating average speed between t and t+h
        Vxm=float((Vx+Vxtph)/2)
        Vym=float((Vy+Vytph)/2)
        Vzm=float((Vz+Vztph)/2)

        #Calculating the new bullet's coordinates
        X=float(Xtph)
        Y=float(Ytph)
        Z=float(Ztph)

        #Next time step
        tph=float(tph+h)

        Vx=float(Vxtph)
        Vy=float(Vytph)
        Vz=float(Vztph)

        #updating impact
        impact.X=X
        impact.Y=Y
        impact.Z=Z
        impact.Vx=Vx
        impact.Vy=Vy
        impact.Vz=Vz

        impact.E_a=math.atan(Vy/Vx)
        impact.W_a=math.atan(Vz/Vx)
        impact.update()
#        impact.show()
        D_Tir_c=impact.X

#    print("End of Solving ")

def PRS_Solver(D_Tir,h,MuzzleSpeed,NS,Az_d,Az,Co,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle_Brand,Rifle_Model,d_i,bl_cal,alpha_d,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO):

    global H_Angle_d,H_Angle,Ytph,tph,Vmod,SD,Sg,Shoot_Wind,bullet,Ytph_HA,Ytph_HA_Co,Ytph_Co,Ytph_NoHa_NoCo,impact,impact_HA,impact_HA_Co,impact_Co,impact_NoHa_NoCo,D_Tir_HA,Arrow,zero,D_Tir_c

    Ytph_HA=float(0)
    D_Tir_HA=float(0)

    if H_Angle_d == 0:

        RK(D_Tir,h,MuzzleSpeed,NS,Az_d,Az,Co,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle_Brand,Rifle_Model,d_i,bl_cal,alpha_d,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO)

        if Co == "Y":
#            print("Impact No Ha with Coriolis")
            impact_Co.X=impact.X
            impact_Co.Y=impact.Y
            impact_Co.Z=impact.Z
            impact_Co.Vx=impact.Vx
            impact_Co.Vy=impact.Vy
            impact_Co.Vz=impact.Vz
            impact_Co.E_a=impact.E_a
            impact_Co.W_a=impact.W_a
            impact_Co.update()
#            impact_Co.show()
            Ytph_Co=float(impact_Co.Y)
        if Co == "N":
#        print("Impact No Ha and No Coriolis")
            impact_NoHa_NoCo.X=impact.X
            impact_NoHa_NoCo.Y=impact.Y
            impact_NoHa_NoCo.Z=impact.Z
            impact_NoHa_NoCo.Vx=impact.Vx
            impact_NoHa_NoCo.Vy=impact.Vy
            impact_NoHa_NoCo.Vz=impact.Vz
            impact_NoHa_NoCo.E_a=impact.E_a
            impact_NoHa_NoCo.W_a=impact.W_a
            impact_NoHa_NoCo.update()
#            impact_NoHa_NoCo.show()
            Ytph_NoHa_NoCo=float(impact_NoHa_NoCo.Y)
    if H_Angle_d > 0:
        D_Tir_HA=D_Tir*math.cos(H_Angle)
        print("PRS_Solver : shoot with Horizontal Angle (deg) :",H_Angle_d)
        print("PRS_Solver : Solving Balistic for corrected distance (m) :",round(D_Tir_HA,2))

        RK(D_Tir_HA,h,MuzzleSpeed,NS,Az_d,Az,Co,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle_Brand,Rifle_Model,d_i,bl_cal,alpha_d,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO)

#        print("Solving ended")
#        impact.show()

        if Co == "N":
#            print("Impact with HA No Coriolis")
            impact_HA.X=impact.X
            impact_HA.Y=impact.Y
            impact_HA.Z=impact.Z
            impact_HA.Vx=impact.Vx
            impact_HA.Vy=impact.Vy
            impact_HA.Vz=impact.Vz
            impact_HA.E_a=impact.E_a
            impact_HA.W_a=impact.W_a
            impact_HA.update()
#            impact_HA.show()
            Ytph_HA=float(impact_HA.Y)
        if Co == "Y":
#               print("Impact with HA and Coriolis")
            impact_HA_Co.X=impact.X
            impact_HA_Co.Y=impact.Y
            impact_HA_Co.Z=impact.Z
            impact_HA_Co.Vx=impact.Vx
            impact_HA_Co.Vy=impact.Vy
            impact_HA_Co.Vz=impact.Vz
            impact_HA_Co.E_a=impact.E_a
            impact_HA_Co.W_a=impact.W_a
            impact_HA_Co.update()
#            impact_HA_Co.show()
            Ytph_HA_Co=float(impact_HA_Co.Y)

        RK(D_Tir,h,MuzzleSpeed,NS,Az_d,Az,Co,m_gr,SightHeight,Twist,T_cal,Twist_C,Rifle_Brand,Rifle_Model,d_i,bl_cal,alpha_d,Alt,Gamma,Alt_ICAO,P_ICAO,T_C_ICAO,T_K_ICAO,RH_ICAO,Rho_ICAO,PVS_ICAO,PV_ICAO,TV_ICAO,P,T_C,T_K,RH,Rho,PVS,PV,TV,Omega,SoundSpeed_ICAO)


def MaxY(LE_c,Prec):
    global Arrow,zero,Y_max,XY_max,TY_max,c_mRAD,D_Tir,h,bullet,NS,Az_d,Az,Co,montage,rifle,Shoot_Atm,Gamma,ICAO_Atm,Omega,SoundSpeed_ICAO,impact

    cpt=1
    Delta_Y=float(1)

    Arrow=1

    alpha_d_backup=zero.alpha_d

    zero.alpha_d=float(zero.alpha_d+((LE_c*c_mRAD)/1000)*(180/math.pi))

    PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)

    zero.alpha_d=alpha_d_backup
    Arrow=0

# Reset Data and import from Files
ResetData()

# Print all Data

print("==============================")
print("Rifle and Scope related parameters ")
print("==============================")

rifle.show()
montage.show()

print("=========================")
print("Bullet related parameters")
print("=========================")

bullet.show()
print("bullet length (inch) : ",round(bullet.bl_i,3))

print("==============================")
print("Earth Localization ")
print("==============================")

print("Latitude ",Lat_d," ° ",Lat_m," min ",Lat_s," s")
print("Latitude degree",round(Lat_D,3))

print("==============================")
print("ICAO Standard Atmosphere ")
print("Hard coded values ")
print("==============================")

ICAO_Atm.show()

print("==============================")
print("Zeroing Atmosphere ")
print("==============================")

# Altitude provided by Smart Phone Compass
# Atmospheric values provided by Kestrel

zero.show()

print("==============================")
print("Shooting Atmosphere ")
print("==============================")

# Altitude provided by Smart Phone Compass
# Atmospheric values provided by Kestrel

Shoot_Atm.show()

print("============================")
print("Shot related parameters ")
print("============================")

print("Shooting Distance : ",int(D_Tir))
print("Time increment (s) : ",h)
print("Shooting Angle relative to Horizontal plan (deg) : ",H_Angle_d," (RAD) : ",round(H_Angle,6))
print("Shooting Angle relative to North (Azimut °) : ",Az_d," (RAD) : ",round(Az,6))

print("======================")
print("Coriolis Data ")
print("======================")

print("Earth Angular Speed - Omega (rad/s) :",Omega)

print("========================================")
print("ICAO Drag Coefficient (Cd) Determination")
print("========================================")

print("Speed of sound ICAO (m/s) :",round(SoundSpeed_ICAO,2))
print("Bullet Stability Factor ICAO ",round(Sg_ICAO,2))

if Sg_ICAO<1:
    print("ICAO Bullet is NOT stable")
if Sg_ICAO<1.5:
    print("ICAO Bullet Marginaly Stable")
if Sg_ICAO >=1.5:
    print("ICAO Stable Bullet")

print("=========================")
print("Wind Speed and Direction ")
print("=========================")

print("wind speed (m/s) :",Shoot_Wind.Ws)
print("wind Angle relative to shooting direction (hour) :",int(Shoot_Wind.Wa_o))

print("======================")
print(" Options choice ")
print("======================")
print("Spin Drift : ",SD_p)
print("Aerodynamic Jump : ",AJ)
print("Corriolis : ",Co)
print("Zeroing : ",Zero_C)
print("Calculate Abacus :",sys.argv[16])

if Zero_C == "Y":
    print("========== ZEROING ============")
    print("Zeroing at : ",zero.D_zero," (m)")
    print("error size :",zero.Prec)
    print("======== ZEROING in Progress... =========")
    zero.alpha_d=Zeroing(zero.D_zero,zero.Prec)
else:
    print("======================")
    print("No Zeroing requested")
    print("======================")

print("Alpha(0) used (deg) : ",round(zero.alpha_d,8),"Windage correction used (cm)",round(zero.wind_cm,8))

print("====================================================================================")
print("Ballistic differential equations being solved numerically using Ruge-Kutta Method...")
print("====================================================================================")

# Simulation without Coriollis and Ha >=0

Co_B=Co
Co="N"

print("Doing a Simulation without Coriollis ")

PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)

#
print("")
#
impact.show()
#
print("")

#updating impact No Coriolis
impact_NoC.X=impact.X
impact_NoC.Y=impact.Y
impact_NoC.Z=impact.Z

impact_NoC.Vx=impact.Vx
impact_NoC.Vy=impact.Vy
impact_NoC.Vz=impact.Vz

impact_NoC.E_a=impact.E_a
impact_NoC.W_a=impact.W_a
impact_NoC.update()

Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))

tph=float(tph-h)

Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))

SD=SpinDrift(Sg,tph)

Shoot_Wind.WDx=float(0)
Shoot_Wind.WDy=float(0)
Shoot_Wind.WDz=float(0)

#Using Dedion Formula
Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))

Co=Co_B

# Simulation with Coriollis if requested
if Co=="Y":
#
    print("Simulation with Coriollis due to chosen option")

    PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)

#
    print("")
#
    print("Coriolis Results")
#
    impact.show()
#
    print("")

    impact_Delta = Impact()

    #updating impact Delta No HA No Coriolis - No HA Coriolis
    impact_Delta.X=impact_NoHa_NoCo.X-impact_Co.X
    impact_Delta.Y=impact_NoHa_NoCo.Y-impact_Co.Y
    impact_Delta.Z=impact_NoHa_NoCo.Z-impact_Co.Z
    impact_Delta.Vx=impact_NoHa_NoCo.Vx-impact_Co.Vx
    impact_Delta.Vy=impact_NoHa_NoCo.Vy-impact_Co.Vy
    impact_Delta.Vz=impact_NoHa_NoCo.Vz-impact_Co.Vz
    impact_Delta.E_a=impact_NoHa_NoCo.E_a-impact_Co.E_a
    impact_Delta.W_a=impact_NoHa_NoCo.W_a-impact_Co.W_a
    impact_Delta.update()
#
    print("Delta No Co - Co")
#
    impact_Delta.show()
#
    print("")

    Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))

    tph=float(tph-h)

    Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))

    SD=SpinDrift(Sg,tph)

    Shoot_Wind.WDx=float(0)
    Shoot_Wind.WDy=float(0)
    Shoot_Wind.WDz=float(0)

    #Using Dedion Formula
    Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
    Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))

print("===================================")
print("Printing All Results")
print("===================================")

print("===================================")
print("Shot Parameters ")
print("===================================")

AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
AJD_y=-AJ_mRAD*Shoot_Wind.Wsz
AJE=(D_Tir*math.tan(AJD_y/1000))

print("Lattitude (° N/S) :",round(Lat_D,2),NS)
print("Shooting Direction (Azimut Angle) related to North (deg) : ",Az_d,"RAD",round(Az,5))
print("Shooting Direction (Horizontal Angle) related to vertical (deg) : ",H_Angle_d)
print("Goal Distance (m) :",round(D_Tir,0))
print("wind speed (m/s) :",Shoot_Wind.Ws)
print("wind Angle relative to shooting direction (hour) :",Shoot_Wind.Wa_o)
print("Time increment (s) : ",h," (ms) :",(h*1000))

print("===========================================")
print("Calculated values not linked to any options")
print("===========================================")

Shoot_Wind.show()

print("Wind Drift Along X (m) :",round(Shoot_Wind.WDx,4)," (cm) :",round(Shoot_Wind.WDx*100,2))
print("Wind Drift Along Y (m) :",round(Shoot_Wind.WDy,4)," (cm) :",round(Shoot_Wind.WDy*100,2))
print("Wind Drift Along Z (m) :",round(Shoot_Wind.WDz,4)," (cm) :",round(Shoot_Wind.WDz*100,2))

print("Calculated Z shift due to Wind Drift (m) :",round((Shoot_Wind.WDz),2)," (cm) :",round((Shoot_Wind.WDz)*100,2))

print("Calculated Z Angle due to Wind Drift (mRAD) :",round(float((math.atan(Shoot_Wind.WDz/D_Tir)*1000)),2))

ToF=float(tph)
print("Time of Flight (s) :",round(ToF,3))

print("Bullet Stability Factor Sg =",round(Sg,2))
if Sg<1:
    print("Sg <1 Bullet is NOT stable")
if Sg<1.5:
    print("0< Sg <1.5 Bullet is Marginaly Stable")
if Sg >=1.5:
    print("Sg >1.5 Bullet is Stable")

Y_Drop=impact_NoC.Y

Y_a_NoHa_NoCo=float((math.atan(impact_NoHa_NoCo.Y/impact_NoHa_NoCo.X)*1000))

Y_a_NoHa_Co=float((math.atan(impact_Co.Y/impact_Co.X)*1000))

Y_a_NoC=float((math.atan(impact_HA.Y/impact_HA.X)*1000))

Y_a_Ha_Co=float((math.atan(impact_HA_Co.Y/impact_HA_Co.X)*1000))

if H_Angle_d == 0:
    print("")
    print("Calculated bullet Impact parameters (With only Drag and Gravity influences without Coriolis and for horizontal shooting Ha=0) ")
    print("")
    print("Calculated Impact Speed Module |V| (m/s) :",round(impact_NoHa_NoCo.V,3))

    print("Calculated Y impact position with Horizontal Angle = 0 (m) :",round(impact_NoHa_NoCo.Y,3)," (cm) :",round(impact_NoHa_NoCo.Y*100,3))

else:
    print("")
    print("Calculated bullet Impact parameters (With only Drag and Gravity influences without Coriolis for NON horizontal shooting Ha>0 ) ")
    print("")
    print("Target Distance corrected for Horizontal Angle (deg) : ",round(H_Angle_d,1)," (m) :",round(D_Tir_HA,0))

    print("Calculated Impact Speed Module |V| (m/s) (With only Drag and Gravity influences) :",round(impact_HA.V,3))

    print("Calculated impact Y position corrected for Non Horizontal Angle No Coriolis (m) :",round(impact_HA.Y,3)," (cm) :",round(impact_HA.Y*100,3))
    Y_Drop_HA=impact_HA.Y
    print("Calculated impact Y position corrected for Non Horizontal Angle With Coriolis (m) :",round(impact_HA_Co.Y,3)," (cm) :",round(impact_HA_Co.Y*100,3))

print("Calculated Z impact position Coriolis Ha ? (m) :",round(impact_NoC.Z,3))

print("Calculated Z impact position Ha and Coriolis  (m) :",round(impact_HA_Co.Z,3))

print("Calculated Y impact Angle No Ha and No Coriolis (mRAD) :",round(Y_a_NoHa_NoCo,3))
print("Calculated Y impact Angle No Ha and Coriolis (mRAD) :",round(Y_a_NoHa_Co,3))
print("Calculated Y impact Angle Ha and No Coriolis (mRAD) :",round(Y_a_NoC,3))
print("Calculated Y impact Angle Ha and Coriolis (mRAD) :",round(Y_a_Ha_Co,3))

print("Elevation to be applied due to gravity drag, No Ha and No Coriolis (clicks) :",round(-Y_a_NoHa_NoCo/c_mRAD,1))
print("Elevation to be applied due to gravity drag No Ha and Coriolis (clicks) :",round(-Y_a_NoHa_Co/c_mRAD,1))
print("Elevation to be applied due to gravity drag Ha and No Coriolis (clicks) :",round(-Y_a_NoC/c_mRAD,1))
print("Elevation to be applied due to gravity drag Ha and Coriolis (clicks) :",round(-Y_a_Ha_Co/c_mRAD,1))

A_SpinDrift=0
if (zero.D_zero != 0 ):
    A_SpinDrift=float(((SD-(zero.wind_cm*(D_Tir/zero.D_zero)))/100))

print("Spin Drift including zero correction (m) : ",round(A_SpinDrift,3),"(cm): ",round((A_SpinDrift*100),2))

W_SpinDrift=-(math.atan(A_SpinDrift/D_Tir)*1000)/c_mRAD

print("Windage correction due to Spin Drift (clicks) :",round(W_SpinDrift,1))

A_AJ=AJE
print("Aerodynamic Jump (m) :",round(A_AJ,3)," (cm) ",round(A_AJ*100,2))

print("Elevation correction due to Aerodynamic Jump (clicks) :",round(float((-math.atan(A_AJ/D_Tir)*1000)/c_mRAD),2))

print("")

E_c_NoC=float(-Y_a_NoHa_NoCo/c_mRAD)
W_c_NoC=float(W_SpinDrift)

print("=========== CORRECTIONS TO BE APPLIED WITHOUT OPTION =================")

print("Elevation (gravity, drag No Ha and No Coriolis) to be applied (clicks) +=>Up -=>Down:",round(E_c_NoC,1))

print("Windage (Spin Drift only including zero correction) to be applied (clicks) +=>Rigt -=>Left:",round(W_c_NoC,1))

print("=======================================================")


print("============================================")
print("Calculated values depending on choosen options")
print("============================================")

# Elevation correction depending on Target Distance (gravity, drag), Range Wind, Horizontal Angle and Coriolis

if Co == "N":
    if H_Angle_d == 0:
        E_c=float(-Y_a_NoHa_NoCo/c_mRAD)
    else:
        E_c=float(-Y_a_NoC/c_mRAD)
    Z_Co=float(0)
else:
    if H_Angle_d == 0:
        E_c=float(-Y_a_NoHa_Co/c_mRAD)
    else:
        E_c=float(-Y_a_Ha_Co/c_mRAD)
    Z_Co=float(impact.Z)

print("Elevation to be applied due to Target Distance (gravity, drag), Range Wind, Horizontal Angle and Coriolis (clicks) :",round(E_c,1))

# Elevation correction depending on Aerodynamic Jump
if AJ == "N":
    AJE=float(0)

# print("Aerodynamic Jump (m) :",round(AJE,3)," (cm) ",round(AJE*100,2))
print("Calculated shift along Y axis due to Aerodynamic Jump (m) :",round(AJE,3)," (cm) :",round(AJE*100,2))
E_c_a_AJ=float((math.atan(AJE/D_Tir)*1000))
print("Calculated Angle along Y axis due to Aerodynamic Jump (mRAD) :",round(E_c_a_AJ,2))
E_c_AJ=float(-(E_c_a_AJ/c_mRAD))
print("Calculated Correction due to Aerodynamic Jump (click) :",round(E_c_AJ,1))

if AJ == "Y":
    E_c+=E_c_AJ
#    print("Angle due to gravity, drag and Aerodynamic Jump (mRAD) :",round(Y_a,3))

print("Elevation to be applied due to due to Target Distance (gravity, drag), Range Wind, Horizontal Angle, Coriolis and Aerodynamic Jump (clicks) :",round(E_c,1))

#Windage calculation
W_c=float(0)
W_c_Co=float(0)

#Windage Due to Coriolis
W_c_Co=-(math.atan(Z_Co/D_Tir)*1000/c_mRAD)

#Windage due to Spin Drift

ZC_SpinDrift=0
if (zero.D_zero != 0 ):
    ZC_SpinDrift=float(((SD-(zero.wind_cm*(D_Tir/zero.D_zero)))/100))

if SD_p == "N":
    SD=float(0)
    ZC_SpinDrift=float(0)
    ZC_SpinDrift_a=float(0)
    W_c_SD=float(0)
if SD_p == "Y":
    ZC_SpinDrift_a=float(math.atan(ZC_SpinDrift/D_Tir)*1000)
    W_c_SD=float(-(ZC_SpinDrift_a/c_mRAD))

print("Calculated Z shift due to Coriolis (m) :",round(Z_Co,5)," (cm) : ",round((Z_Co*100),2))

print("Windage to be applied due to due to Coriolis (clicks) :",round(W_c_Co,1))
W_c+=W_c_Co

print("Calculated Z shift due to Spin Drift including zero correction (m) :",round(ZC_SpinDrift,5)," (cm) : ",round((ZC_SpinDrift*100),2))

print("Windage to be applied due to due to Spin Drift (clicks) :",round(W_c_SD,1))
W_c+=W_c_SD

print("Calculated Z shift due to Cross Wind (m) :",round(Shoot_Wind.WDz,5)," (cm) : ",round((Shoot_Wind.WDz*100),2))

Wind_CW_a=float(math.atan(Shoot_Wind.WDz/D_Tir)*1000)

print("Calculated Z Angle due to Cross Wind (mRAD) :",round(Wind_CW_a,2))

Wind_CW_a_c=-(Wind_CW_a/c_mRAD)

print("Windage to be applied due to due to Cross Wind (clicks) :",round(Wind_CW_a_c,1))

W_c+=Wind_CW_a_c

print("Windage to be applied due to due to Spin Drift and Cross Wind (clicks) :",round(W_c,1))

J=float((ICAO_Atm.P/Shoot_Atm.P)*math.sqrt(Shoot_Atm.TV/ICAO_Atm.TV))

print("Impedance multiplicator  ",round(J,3))


# Mach_ICAO=float(bullet.V0/SoundSpeed_ICAO)
Mach_ICAO=float(amo.MuzzleSpeed/SoundSpeed_ICAO)

Cd_G1_ICAO=CdG1Search(Mach_ICAO)
Cd_G7_ICAO=CdG7Search(Mach_ICAO)

match BC_Gx:
    case "G1":
        BC_G1_ICAO=BC_Search(bullet.V0)
        BC_G7_ICAO=float(BC_G1_ICAO/(Cd_G1_ICAO/Cd_G7_ICAO))
    case "G7":
        BC_G7_ICAO=BC_Search(bullet.V0)
#        print("Mach_ICAO :",Mach_ICAO," BC_G7_ICAO : ",BC_G7_ICAO," Cd_G1_ICAO :",Cd_G1_ICAO," Cd_G7_ICAO :",Cd_G7_ICAO)
        BC_G1_ICAO=float(BC_G7_ICAO/(Cd_G7_ICAO/Cd_G1_ICAO))

BC_G1=float(BC_G1_ICAO*J)
print("At Muzzle Speed")
print("Ballistic Coefficient G1 ICAO",round(BC_G1_ICAO,3),"Ballistic Coefficient G1 Current Atm",round(BC_G1,3))
BC_G7=float(BC_G7_ICAO*J)
print("Ballistic Coefficient G7 ICAO",round(BC_G7_ICAO,3),"Ballistic Coefficient G7 Current Atm",round(BC_G7,3))

print("============================================================")
print("Calculation of the maximum value for Y along the trajectory ")
print("============================================================")


MaxY(E_c, zero.Prec)

#impact.show()

print("Max Z (m):",round(Y_max,4)," for distance (m) :",round(XY_max,3)," at time (s) :",round(TY_max,3))

print("")

print("=========== CORRECTIONS TO BE APPLIED =================")

print("Elevation to be applied (clicks) +=>Up -=>Down:",round(E_c,1))

print("Windage to be applied (clicks) +=>Rigt -=>Left:",round(W_c,1))

print("=======================================================")

#
if sys.argv[16] == "N":
    time.sleep(1)
    exit(0)

if sys.argv[16] == "C":
    print("===================================")
    print("Creating the Shot Card  ")
    print("===================================")

    FNoutput= "./output-"+str(int(bullet.d_i*1000))+"_"+amo.Amo_Brand+"_"+str(int(bullet.m_gr))+"_"+str(int(D_Tir))+"_"+str(todays_date)+".html"
    print("Output File Name : ",FNoutput)

    OutputFile = open(FNoutput, "w")

    Total_E=float(0)
    Total_W=float(0)

    OutputFile.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.x//EN\"\n")
    OutputFile.write("\"http://www.w3.org/TR/html4/strict.dtd\">\n")
    OutputFile.write("<html>\n")
    OutputFile.write("<head>\n")
    OutputFile.write("<style>\n")
    OutputFile.write("table, tr, th, td {\n")
    OutputFile.write("\tborder: 1px solid black;\n")
    OutputFile.write("\tborder-collapse: collapse;\n")
    OutputFile.write("}\n")
    OutputFile.write("<title>PBS Output</title>\n")
    OutputFile.write("</style>\n")
    OutputFile.write("</head>\n")
    OutputFile.write("<body>\n")
    OutputFile.write("<table>\n")
    OutputFile.write("\t<tr>\n")
    OutputFile.write("\t\t<th colspan=\"8\">")

    Ligne1="PBS v"+PBS_version+" "+PBS_year+" Shooting Card "

    OutputFile.write(Ligne1)
    OutputFile.write("</th>\n")
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")

    Ligne2="\t\t<td align=\"center\" colspan=\"8\">"+" RIFLE : "+rifle.Rifle_Brand+" "+rifle.Rifle_Model+" "+str(round(bullet.d_i,3))+" (inch)"+" - Rifle Bore "+rifle.Twist_C+"Twist 1:"+str(round(rifle.Twist,0))+" (inch)"+"- Sight Height : "+str(montage.SightHeight*1000)+" (mm)"+"</td>\n"
    OutputFile.write(Ligne2)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")
    Ligne200="\t\t<td align=\"center\" colspan=\"8\">"+" BULLET : "+str(round(bullet.d_i,3))+" (inch)"+str(int(bullet.m_gr))+"(gr) "+bullet.Bullet_Brand+" "+bullet.Bullet_Model+" - Muzzle Speed "+str(int(bullet.MuzzleSpeed))+" (m/s) in ICAO Atmosphere "+" - Ballistic Coefficient in current conditions : G1 "+str(round(BC_G1,3))+" - G7 "+str(round(BC_G7,3))+"</td>\n"
    OutputFile.write(Ligne200)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")
    if Sg<1:
        Stability="Sg <1 Bullet is NOT stable"
    if Sg<1.5:
        Stability="0< Sg <1.5 Bullet is Marginaly Stable"
    if Sg >=1.5:
        Stability="Sg >1.5 Bullet is Stable"
    Ligne300="\t\t<td align=\"center\" colspan=\"8\">"+" Time Of Flight (s) : "+str(round(ToF,2))+" Bullet Stability Factor Sg = "+str(round(Sg,2))+" - "+Stability+"</td>\n"
    OutputFile.write(Ligne300)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")
    Ligne400="\t\t<td align=\"center\" colspan=\"8\">"+" ALTITUDE (m) : "+str(Shoot_Atm.Alt)+" LATITUDE (°) : "+str(round(Lat_D,2))+" "+NS+" - Shooting Direction, relative to the North (°): "+str(Az_d)+" - Coriolis effects ="+Co+"</td>\n"
    OutputFile.write(Ligne400)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")
    Ligne3="\t\t<th align=\"center\" colspan=\"3\">"+" ELEVATION "+"</th>\n"
    OutputFile.write(Ligne3)

    Ligne4="\t\t<th align=\"center\" colspan=\"5\">"+" WINDAGE "+"</th>\n"
    OutputFile.write(Ligne4)

    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")
    Ligne5="\t\t<th align=\"center\" colspan=\"2\">"+" DATA "+"</th>\n"+"\t\t<th align=\"center\" colspan=\"1\">"+" CLICKS "+"</th>\n""\t\t<th align=\"center\" colspan=\"4\">"+" DATA "+"</th>\n"+"\t\t<th align=\"center\" colspan=\"1\">"+" CLICKS "+"</th>\n"
    OutputFile.write(Ligne5)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")

    if H_Angle_d==0:
        if Co == "Y":
            Total_E+=-round((Y_a_NoHa_Co/c_mRAD),1)
        if Co == "N":
            Total_E+=-round((Y_a_NoHa_NoCo/c_mRAD),1)
        Delta_Y_Drop=0
    if H_Angle_d > 0:
        if Co == "Y":
            Total_E+=-round((Y_a_Ha_Co/c_mRAD),1)
        if Co == "N":
            Total_E+=-round((Y_a_NoC/c_mRAD),1)
        Delta_Y_Drop=(Y_Drop-Y_Drop_HA)
        print("Shooting Angle correction : ",(math.atan(Delta_Y_Drop/D_Tir_c)*1000/c_mRAD))
    Partial_E=float(Total_E-(math.atan(Delta_Y_Drop/D_Tir_c)*1000/c_mRAD))

    Ligne6 = "\t\t<td align=\"center\">"+" TARGET DISTANCE (m) Include Gravity,Drag,Range Wind, Coriolis & Non-ICAO Influences"+"</td>"+"\t\t<td align=\"center\">"+str(int(D_Tir))+"</td>"+"\t\t<td align=\"center\">"+str(round((Partial_E),1))+"</td>"+"\t\t<td align=\"center\" colspan=\"1\">"+" BULLET MAX HEIGHT (m) "+"</td>"+"\t\t<td align=\"center\">"+ str(round(Y_max,2))+"</td>"+"\t\t<td align=\"center\" colspan=\"1\">"+" @ DISTANCE (m) "+"</td>"+"\t\t<td align=\"center\">"+ str(round(XY_max,1))+"</td>"+"\t\t<td align=\"center\" colspan=\"1\">"+" N A "+"</td\n>"

    OutputFile.write(Ligne6)
    OutputFile.write("\t</tr>\n")

    if H_Angle_d==0:
        Delta_Y_Drop=0
    else:
        Delta_Y_Drop=(Y_Drop-Y_Drop_HA)

    OutputFile.write("\t<tr>\n")
    Ligne7 = "\t\t<td align=\"center\">"+" SHOOTING ANGLE, relative to the horizontal (°) "+"</td>"+"\t\t<td align=\"center\">"+str(round(H_Angle_d,2))+"</td>"+"\t\t<td align=\"center\">"+str(round((math.atan(Delta_Y_Drop/D_Tir_c)*1000/c_mRAD),1))+"</td>"+"\t\t<td align=\"center\" colspan=\"4\">"+" SPIN DRIFT "+"</td>"+"\t\t<td align=\"center\">"+ str(round(W_SpinDrift,1)) +"</td\n>"
    OutputFile.write(Ligne7)
    OutputFile.write("\t</tr>\n")

    Total_W+=W_SpinDrift

    OutputFile.write("\t<tr>\n")
    Ligne8 = "\t\t<td align=\"center\">"+" P (hPa), 1013.25 is the reference "+"</td>"+"\t\t<td align=\"center\">"+str(round(Shoot_Atm.P/100,2))+"</td>"+"\t\t<td align=\"center\">"+" N A "+"</td>"+"\t\t<td align=\"center\">"+"WIND Direction (hour)  "+"</td>"+"\t\t<td align=\"center\">"+str(round(Shoot_Wind.Wa_o,1))+"</td>"+"\t\t<td align=\"center\">"+"WIND Speed (m/s)  "+"</td>"+"\t\t<td align=\"center\">"+str(round(Shoot_Wind.Ws,1))+"\t\t<td align=\"center\">"+" N A "+"</td\n>"
    OutputFile.write(Ligne8)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")
    Ligne9 = "\t\t<td align=\"center\">"+" T (°C), 15°C is the reference "+"</td>"+"\t\t<td align=\"center\">"+str(round(Shoot_Atm.T_C,1))+"</td>"+"\t\t<td align=\"center\">"+" N A "+"</td>"+"\t\t<td align=\"center\" colspan=\"5\">"+"   "+"</td\n>"

    OutputFile.write(Ligne9)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("\t<tr>\n")
    Ligne10 = "\t\t<td align=\"center\" colspan=\"2\">"+" RANGE WIND INFLUENCE "+"</td>"+"\t\t<td align=\"center\">"+" N A "+"</td>"+"\t\t<td align=\"center\" colspan=\"4\">"+" CROSS WIND INFLUENCE "+"</td>"+"\t\t<td align=\"center\">"+str(round((-math.atan(Shoot_Wind.WDz/D_Tir)*1000/c_mRAD),1))+"</td\n>"
    OutputFile.write(Ligne10)
    OutputFile.write("\t</tr>\n")

    Total_W+=-(math.atan(Shoot_Wind.WDz/D_Tir)*1000/c_mRAD)

    OutputFile.write("\t<tr>\n")
    Ligne11 = "\t\t<td align=\"center\" colspan=\"2\">"+" AERODYNAMIC JUMP (Due to Range Wind)"+"</td>"+"\t\t<td align=\"center\">"+str(round((-math.atan(A_AJ/D_Tir_c)*1000/c_mRAD),1))+"</td>"+"\t\t<td align=\"center\" colspan=\"4\">"+"Coriolis Lateral impact"+"</td>"+"\t\t<td align=\"center\">"+str(round(-(math.atan(Z_Co/D_Tir)*1000/c_mRAD),1))+"</td\n>"

    OutputFile.write(Ligne11)
    OutputFile.write("\t</tr>\n")

    Total_W+=-(math.atan(Z_Co/D_Tir)*1000/c_mRAD)

    Total_E+=round(-(math.atan(A_AJ/D_Tir_c)*1000/c_mRAD),1)

    OutputFile.write("\t<tr>\n")
    Ligne12 = "\t\t<th colspan=\"2\">"+" TOTAL ELEVATION "+"</th>"+"\t\t<td align=\"center\">"+str(round(Total_E,1))+"</td>"+"\t\t<th colspan=\"4\">"+"TOTAL WINDAGE : "+"\t\t<td align=\"center\">"+str(round(Total_W,1))+"</th\n>"

    OutputFile.write(Ligne12)
    OutputFile.write("\t</tr>\n")

    OutputFile.write("</table>\n")
    OutputFile.write("</body>\n")
    OutputFile.write("</html>\n")
    OutputFile.close()

else:
    print("===================================")
    print("Creating the Abacus ")
    print("===================================")

    FNabacus= "./Abacus-"+str(int(bullet.d_i*1000))+"_"+amo.Amo_Brand+"_"+str(int(bullet.m_gr))+"_"+str(int(D_Tir))+"_"+str(todays_date)+".html"
    print("Abacus File Name : ",FNabacus)

    AbacusFile = open(FNabacus, "w")

    AbacusFile.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.x//EN\"\n")
    AbacusFile.write("\"http://www.w3.org/TR/html4/strict.dtd\">\n")
    AbacusFile.write("<html>\n")
    AbacusFile.write("<head>\n")
    AbacusFile.write("<style>\n")
    AbacusFile.write("table, tr, th, td {\n")
    AbacusFile.write("\tborder: 1px solid black;\n")
    AbacusFile.write("\tborder-collapse: collapse;\n")
    AbacusFile.write("}\n")
    AbacusFile.write("<title>PBS Abacus</title>\n")
    AbacusFile.write("</style>\n")
    AbacusFile.write("</head>\n")
    AbacusFile.write("<body>\n")
    AbacusFile.write("<table>\n")
    AbacusFile.write("\t<tr>\n")
    AbacusFile.write("\t\t<th colspan=\"16\">")

    Ligne1="PBS v"+PBS_version+" "+PBS_year+" Generic Abacus - "+str(int(bullet.d_i*1000))+"(inch) "+str(int(bullet.m_gr))+"(gr) "+" Bullet - "+bullet.Bullet_Brand+" "+bullet.Bullet_Model+" - Rifle Bore "+rifle.Twist_C+"Twist 1:"+str(round(rifle.Twist,0))+" (inch) - Muzzle Speed "+str(int(bullet.MuzzleSpeed))+" (m/s) in ICAO Atmosphere "+"- Sight Height : "+str(montage.SightHeight*1000)+" (mm)"+" - Ballistic Coefficient in current conditions : G1 "+str(round(BC_G1,3))+" - G7 "+str(round(BC_G7,3))+" - Time of Flight (s) "+str(round(ToF,3))

    AbacusFile.write(Ligne1)
    AbacusFile.write("</th>\n")

    print("===================================")
    print("Summary ")
    print("===================================")
    print("PBS Generic Abacus :",bullet.d_i,"(inch)",bullet.m_gr,"(gr)","Bullet :",bullet.Bullet_Brand,"-",bullet.Bullet_Model,"-",rifle.Twist_C,"Twist 1:",round(rifle.Twist,0),"(inch) Muzzle Speed",bullet.MuzzleSpeed,"(m/s) in ICAO Atmosphere","- Sight Height :",montage.SightHeight*1000,"(mm)")

    print("===================================")
    print("Shooting parameters ")
    print("===================================")

    A_D_Tir_R=D_Tir
    print("Shooting distance (m) :",int(A_D_Tir_R))

    A_H_Angle_d_R=H_Angle_d
    print("Shooting vertical angle (deg) :",A_H_Angle_d_R)

    A_Shoot_AtmP_R=Shoot_Atm.P
    print("Shooting absolute pressure (hPa) :",round(A_Shoot_AtmP_R/100,2))

    A_Shoot_AtmT_C_R=Shoot_Atm.T_C
    print("Shooting temperature (°C) :",int(A_Shoot_AtmT_C_R))

    A_Ws_R=Shoot_Wind.Ws
    A_Wa_o_R=Shoot_Wind.Wa_o
    print("Wind intensity (m/s) :",round(A_Ws_R,1)," direction (hour) :",round(A_Wa_o_R,1))

#    print("Simulation ")

#    print("===================================")
#    print("Results ")
#    print("===================================")

    A_ToF_R=tph
    print("Time of Flight (s) :",round(A_ToF_R,2))
    A_Y_max_R=Y_max
    A_XY_max_R=XY_max
    A_TY_max_R=TY_max

#    print("===================================")
#    print("Elevation ")
#    print("===================================")

    A_E_R=E_c
#
    print("Elevation to be applied due to gravity and drag (clicks) Drum +=>Up -=>Down:",round(A_E_R,0))

    A_AJ_R=-((A_AJ*1000/D_Tir_c)/c_mRAD)

#    print("===================================")
#    print("Windage ")
#    print("===================================")

    A_Wc_R=W_c
    A_SpinDrift_R=-((A_SpinDrift*1000/A_D_Tir_R)/c_mRAD)

#    print("===================================")
#
    print("Simulation for Distances ")
#    print("===================================")

    Ligne2=""
    Ligne3=""
    Ligne4=""
    Ligne5=""

    A_D_Tir=A_D_Tir_R

#
    print("Simulation for Vertical Angles ")

    while A_D_Tir < (A_D_Tir_R+100):

#        print("===================================")
#        print("Shooting distance (m) ",A_D_Tir)
#        print("================================")

#        print("================================")

        A_H_Angle_d=float(0)

        while A_H_Angle_d <= 30:

            A_H_Angle=float(A_H_Angle_d*(math.pi/180))
            A_D_Tir_HA=A_D_Tir*math.cos(A_H_Angle)

            PRS_Solver(A_D_Tir_HA,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)

            Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))

            tph=float(tph-h)

            Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))

            SD=SpinDrift(Sg,tph)

            WDx=float(0)
            WDy=float(0)
            WDz=float(0)

            WDy=float(Shoot_Wind.Wsy*(tph-(round(A_D_Tir,0)/bullet.V0)))
            WDz=float(Shoot_Wind.Wsz*(tph-(round(A_D_Tir,0)/bullet.V0)))

            Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
            A_E=-Y_a/c_mRAD

            if A_H_Angle_d == 0:
                A_Delta_E=0
                A_E_RVA=A_E
                Ligne4+="\t\t<th>"+str(int(A_H_Angle_d))+"</th>\n"
                Ligne5+="\t\t<td align=\"center\">"+str(int(A_Delta_E))+"</td>\n"
            if A_H_Angle_d > 0:
                A_Delta_E = A_E - A_E_RVA
                Ligne4+="\t\t<th>"+str(int(A_H_Angle_d))+"</th>\n"
                Ligne5+="\t\t<td align=\"center\">"+str(int(A_Delta_E))+"</td>\n"
            A_H_Angle_d+=10
#            print("==============================")

        Ligne2+="\t\t<td align=\"center\" colspan=\"4\">"+str(int(A_D_Tir))+" = "+str(int(round(A_E_RVA,0)))+"</td>\n"
        Ligne3+="\t\t<td align=\"center\" colspan=\"4\">"+"Vertical Shooting Angle (deg)"+"</td>\n"

        A_D_Tir+=25

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne2)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne3)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne4)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne5)
    AbacusFile.write("\t</tr>\n")

#    print("===================================")
#
    print("Simulation for Pressure ")
#    print("===================================")
#    print("Shooting distance (m) ",D_Tir)
#    print("===================================")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write("\t\t<th colspan=\"16\">")
    Ligne10="Local Absolute Pressure (hPA)"
    AbacusFile.write(Ligne10)
    AbacusFile.write("</th>\n")

    Ligne20=""
    Ligne30=""
    Ligne40=""
    Ligne50=""

    A_P_B=Shoot_Atm.P
    A_P=float(108800)

    while A_P >= 86300:
        Shoot_Atm.P=A_P
        Shoot_Atm.update()

        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)

        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD
        A_Delta_E = A_E - A_E_R
        Ligne40+="\t\t<th>"+str(int(round((Shoot_Atm.P/100),2)))+"</th>\n"
        Ligne50+="\t\t<td align=\"center\">"+str(int(A_Delta_E))+"</td>\n"

        A_P-=1500

#        print("==============================")

    Shoot_Atm.P=A_P_B
    Shoot_Atm.update()

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne40)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne50)
    AbacusFile.write("\t</tr>\n")

#    print("===================================")
#
    print("Simulation for Temperature ")
#    print("===================================")
#    print("Shooting distance (m) ",D_Tir)
#    print("Local absolute Pressure (hPa) ",round((Shoot_Atm.P/100),2))
#    Shoot_Atm.show()
#    print("===================================")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write("\t\t<th colspan=\"16\">")
    Ligne100="Air Temperature (°C)"
    AbacusFile.write(Ligne100)
    AbacusFile.write("</th>\n")

    Ligne200=""
    Ligne300=""
    Ligne400=""
    Ligne500=""
    Ligne600=""
    Ligne700=""

    A_T_B=Shoot_Atm.T_C
    A_T=float(57.5)

    while A_T >= -60:

        Shoot_Atm.T_C=A_T
        Shoot_Atm.update()

        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)

        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD
        A_Delta_E = A_E - A_E_R

        if Shoot_Atm.T_C >=20:
            Ligne200+="\t\t<th>"+str(round((Shoot_Atm.T_C),1))+"</th>\n"
            Ligne300+="\t\t<td align=\"center\">"+str(int(A_Delta_E))+"</td>\n"
        else:
            if Shoot_Atm.T_C >=-20:
                Ligne400+="\t\t<th>"+str(round((Shoot_Atm.T_C),1))+"</th>\n"
                Ligne500+="\t\t<td align=\"center\">"+str(int(A_Delta_E))+"</td>\n"
            else:
                Ligne600+="\t\t<th>"+str(round((Shoot_Atm.T_C),1))+"</th>\n"
                Ligne700+="\t\t<td align=\"center\">"+str(int(A_Delta_E))+"</td>\n"

        A_T-=2.5

#        print("==============================")

    Shoot_Atm.T_C=A_T_B
    Shoot_Atm.update()

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne200)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne300)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne400)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne500)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne600)
    AbacusFile.write("\t</tr>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write(Ligne700)
    AbacusFile.write("\t</tr>\n")

#    print("===================================")
#
    print("Simulation for Wind ")
#    print("===================================")
#    print("Shooting distance (m) ",D_Tir)
#    print("Wind Speed (m/s) ",round(Shoot_Wind.Ws,2))
#    print("Wind Direction (h) ",round(Shoot_Wind.Wa_o,2))
#    Shoot_Wind.show()
#    print("===================================")

    ResetData()

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write("\t\t<th colspan=\"8\">")
    Ligne100=" Wind Speed (m/s) - Wind Direction (hour) ->"
    AbacusFile.write(Ligne100)
    AbacusFile.write("</th>\n")

    Ligne200=""
    Ligne200+="\t\t<th>"+"I / V"+"</th>\n"
    Ligne200+="\t\t<th>"+"II / IV"+"</th>\n"
    Ligne200+="\t\t<th>"+"III"+"</th>\n"
    Ligne200+="\t\t<th>"+"VI"+"</th>\n"
    Ligne200+="\t\t<th>"+"IX"+"</th>\n"
    Ligne200+="\t\t<th>"+"VIII / X"+"</th>\n"
    Ligne200+="\t\t<th>"+"VII / XI"+"</th>\n"
    Ligne200+="\t\t<th>"+"XII"+"</th>\n"
    AbacusFile.write(Ligne200)
    AbacusFile.write("\t</tr>\n")


    Shoot_Wind.Ws=2

    while Shoot_Wind.Ws <= 10:
        AbacusFile.write("\t<tr>\n")
        AbacusFile.write("\t\t<th colspan=\"8\">")
        Ligne100=str(round(Shoot_Wind.Ws,0))
        AbacusFile.write(Ligne100)
        AbacusFile.write("</th>\n")

        Ligne200=""

    # Simulation for I / V hour
        Shoot_Wind.Wa_o=1
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)
        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))
        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

    # Simulation for II / IV hour
        Shoot_Wind.Wa_o=2
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)
        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))
        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

    # Simulation for III hour
        Shoot_Wind.Wa_o=3
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)

        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))

        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

    # Simulation for VI hour
        Shoot_Wind.Wa_o=6
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)
        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))
        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

    # Simulation for IX hour
        Shoot_Wind.Wa_o=9
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)
        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))
        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

    # Simulation for VIII / X hour
        Shoot_Wind.Wa_o=8
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)
        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))
        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

    # Simulation for VII / XI hour
        Shoot_Wind.Wa_o=7
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)
        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))
        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

    # Simulation for XII hour
        Shoot_Wind.Wa_o=12
        Shoot_Wind.update()
        PRS_Solver(D_Tir,h,bullet.MuzzleSpeed,NS,Az_d,Az,Co,bullet.m_gr,montage.SightHeight,rifle.Twist,rifle.T_cal,rifle.Twist_C,rifle.Rifle_Brand,rifle.Rifle_Model,bullet.d_i,bullet.bl_cal,zero.alpha_d,Shoot_Atm.Alt,Gamma,ICAO_Atm.Alt,ICAO_Atm.P,ICAO_Atm.T_C,ICAO_Atm.T_K,ICAO_Atm.RH,ICAO_Atm.Rho,ICAO_Atm.PVS,ICAO_Atm.PV,ICAO_Atm.TV,Shoot_Atm.P,Shoot_Atm.T_C,Shoot_Atm.T_K,Shoot_Atm.RH,Shoot_Atm.Rho,Shoot_Atm.PVS,Shoot_Atm.PV,Shoot_Atm.TV,Omega,SoundSpeed_ICAO)
        Sg=float((30*bullet.m_gr)/(math.pow(rifle.T_cal,2)*math.pow(bullet.d_i,3)*bullet.bl_cal*(1+math.pow(bullet.bl_cal,2))*(math.pow((bullet.V0/853.4),(1/3))*((Shoot_Atm.T_K*101325)/(288.15*Shoot_Atm.P)))))
        tph=float(tph-h)
        Vmod=float(math.sqrt(Vxtph*Vxtph+Vytph*Vytph+Vztph*Vztph))
        SD=SpinDrift(Sg,tph)
        Shoot_Wind.WDx=float(0)
        Shoot_Wind.WDy=float(0)
        Shoot_Wind.WDz=float(0)
        #Using Dedion Formula
        Shoot_Wind.WDy=float(Shoot_Wind.Wsy*(tph-(round(D_Tir,0)/bullet.V0)))
        Shoot_Wind.WDz=float(Shoot_Wind.Wsz*(tph-(round(D_Tir,0)/bullet.V0)))
        Y_a=float(math.atan(Ytph/D_Tir_c)*1000)
        A_E=-Y_a/c_mRAD

        A_Delta_E = A_E - A_E_R
        Z_a=float((math.atan((Shoot_Wind.WDz+Ztph)/D_Tir))*1000)
        W_c=float(-Z_a/c_mRAD)

        AJ_MOA=0.01*Sg-0.0024*bullet.bl_cal+0.032
        AJ_mRAD=AJ_MOA*(((1/60)*(math.pi/180)*1000)/(1609.44/3600))
        AJD_z=-AJ_mRAD*Shoot_Wind.Wsz
        AJE=(D_Tir*math.tan(AJD_z/1000))

        Ligne200+="\t\t<td align=\"center\">"+"E: "+str(round(A_Delta_E,0))+" W: "+str(round(W_c,0))+" AJ: "+str(round(float((((-AJE)*1000)/D_Tir)/c_mRAD),0))+"</td>\n"

        AbacusFile.write(Ligne200)
        AbacusFile.write("\t</tr>\n")

        Shoot_Wind.Ws+=2

    AbacusFile.write("\t<tr>\n")

    AbacusFile.write("\t\t<th colspan=\"16\">")

    Ligne1000="Spin Drift (click) : "+str(int(round(A_SpinDrift_R,0)))+" Maximum Y (m) : "+str(round(A_Y_max_R,2))+" At (m) : "+str(round(A_XY_max_R,0))+" Time to get there (s) : "+str(round(A_TY_max_R,2))

    AbacusFile.write(Ligne1000)
    AbacusFile.write("</th>\n")

    AbacusFile.write("\t<tr>\n")
    AbacusFile.write("\t\t<th colspan=\"16\">")
    Ligne2000="How to use this Abacus ? Read HowToPBS_Abacus.pdf in  https://github.com/fabienfigueras/TLD"
    AbacusFile.write(Ligne2000)
    AbacusFile.write("</th>\n")

    AbacusFile.write("</table>\n")
    AbacusFile.write("</body>\n")
    AbacusFile.write("</html>\n")
    AbacusFile.close()