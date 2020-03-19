
# coding: utf-8

# In[ ]:


from arcpy import env
from arcpy.sa import*

# env.workspace =  r"C:\Users\gis\Result_folder"
env.workspace = r"C:\Users\gis\LC08_L1TP_167051_20170409_20170414_01_T1"

#This tool calculates the LAnd Surface temperature(LST) from landsat 8. It has an option
#to produce either a seperate LST for Band 10 and Band 11 or an average of these two Bands


# arcpy.CheckOutExtension("spatial")

# arcpy.env.overwriteOutput = True
B10 = arcpy.GetParameterAsText(0)
B11 = arcpy.GetParameterAsText(1)
B6 = arcpy.GetParameterAsText(2)
NIR = arcpy.GetParameterAsText(3)
Red = arcpy.GetParameterAsText(4)
average = arcpy.GetParameterAsText(5)
Landsat_7_Or_below = arcpy.GetParameterAsText(6)
workspace = arcpy.GetParameterAsText(7)


def TOA(B10):
    
    #Calculation of TOA (Top of Atmospheric) spectral radiance.
    #TOA (L) = ML * Qcal + AL
    #ML = Band-specific multiplicative rescaling factor from the metadata (RADIANCE_MULT_BAND_x,
    #where x is the band number).
    #Qcal = corresponds to band 10.
    #AL = Band-specific additive rescaling factor from the metadata (RADIANCE_ADD_BAND_x,
    #where x is the band number).
    # therefore the equation is TOA = 0.0003342 * “Band 10” + 0.1
    
    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    arcpy.SetProgressorLabel("Calculating Radiance for Band_10...")
    
    Outfolder = os.path.join(workspace,"Radiance_B10.TIF")
    LMAX = (0.0003342)
    LMIN = (0.1)

    ToA = arcpy.sa.Float(LMAX * Raster(B10)+ LMIN) 
    ToA.save(Outfolder)
    return Outfolder

def TOA11(B11):

    #Calculation of TOA (Top of Atmospheric) spectral radiance.
    #TOA (L) = ML * Qcal + AL
    #ML = Band-specific multiplicative rescaling factor from the metadata (RADIANCE_MULT_BAND_x,
    #where x is the band number).
    #Qcal = corresponds to band 10.
    #AL = Band-specific additive rescaling factor from the metadata (RADIANCE_ADD_BAND_x,
    #where x is the band number).
    # therefore the equation is TOA = 0.0003342 * “Band 10” + 0.1

    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    arcpy.SetProgressorLabel("Calculating Radiance for Band_11...")
    
    Outfolder = os.path.join(workspace,"Radiance_B11.TIF")
    
    LMAX = (0.0003342)
    LMIN = (0.1)

    ToA = arcpy.sa.Float(LMAX * Raster(B11)+ LMIN) 
    ToA.save(Outfolder)
    return Outfolder
    

def BT(TO):
    
    #TOA to Brightness Temperature conversion
    #BT = (K2 / (ln (K1 / L) + 1)) − 273.15    import os
    #K1 = Band-specific thermal conversion constant from the metadata
    #(K1_CONSTANT_BAND_x, where x is the thermal band number).
    #K2 = Band-specific thermal conversion constant from the metadata
    #(K2_CONSTANT_BAND_x, where x is the thermal band number).
    # L = TOA
    #Therefore, to obtain the results in Celsius, the radiant temperature is
    #adjusted by adding the absolute zero (approx. -273.15°C).
    #BT = (1321.0789 / Ln ((774.8853 / “%TOA%”) + 1)) – 273.15

    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    TOA = Raster(TO)

    arcpy.SetProgressorLabel("Calculating Brightness temperature for Band_10...")

    Outfolder = os.path.join(workspace,"BrighnessT_B10.TIF")
    K2 = float(1321.0789)
    K1 = float(774.8853)
    fh = float(273.15)
    
    Logs = arcpy.sa.Ln((K1 / TOA) + 1)
    div = arcpy.sa.Float(K2 / Logs - fh)

    div.save(Outfolder)
    return  Outfolder

def BT11(TO):

    #TOA to Brightness Temperature conversion
    #BT = (K2 / (ln (K1 / L) + 1)) − 273.15    import os
    #K1 = Band-specific thermal conversion constant from the metadata
    #(K1_CONSTANT_BAND_x, where x is the thermal band number).
    #K2 = Band-specific thermal conversion constant from the metadata
    #(K2_CONSTANT_BAND_x, where x is the thermal band number).
    # L = TOA
    #Therefore, to obtain the results in Celsius, the radiant temperature is
    #adjusted by adding the absolute zero (approx. -273.15°C).
    #BT = (1321.0789 / Ln ((774.8853 / “%TOA%”) + 1)) – 273.15
    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    TOA = Raster(TO)

    arcpy.SetProgressorLabel("Calculating Brightness temperature for Band_11...")

    Outfolder = os.path.join(workspace,"BrighnessT_B11.TIF")
    K2_11 = float(1201.1442)
    K1_11 = float(480.8883)
    fh = float(273.15)
    
    Logs = arcpy.sa.Ln((K1_11 / TOA) + 1)
    div = arcpy.sa.Float(K2_11 / Logs - fh)

    div.save(Outfolder)
    return  Outfolder


def NDVI (NIR, Red):

    #Calculate the NDVI
    #NDVI = (Band 5 – Band 4) / (Band 5 + Band 4)
    #Note that the calculation of the NDVI is important because, subsequently,
    #the proportion of vegetation (Pv), which is highly related to the NDVI,
    #and emissivity (ε), which is related to the Pv, must be calculated.
    #NDVI = Float(Band 5 – Band 4) / Float(Band 5 + Band 4)
    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    arcpy.SetProgressorLabel("Calculating NDVI...")

    Outfolder = os.path.join(workspace,"NDVI.TIF")

    Num = arcpy.sa.Float(Raster(NIR) - Raster(Red))
    Denom = arcpy.sa.Float(Raster(NIR) + Raster(Red))
    NIR_eq = arcpy.sa.Divide(Num, Denom)

    NIR_eq.save(Outfolder)
    return Outfolder

def PV (NDV):

    #Calculate the proportion of vegetation Pv
    #Pv = Square ((NDVI – NDVImin) / (NDVImax – NDVImin))
    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    arcpy.SetProgressorLabel("Calculating Proporsion of Vegetation...")

    Outfolder = os.path.join(workspace,"ProporV.TIF")
    
    NDVI = Raster(NDV)
    out_ras = arcpy.sa.Float((NDVI - NDVI.minimum)  / (NDVI.maximum - NDVI.minimum))**2

    out_ras.save(Outfolder)
    return Outfolder

def Emissiv(Pv):

    #Calculate Emissivity ε
    #ε = 0.004 * Pv + 0.986
    #Simply apply the formula in the raster calculator, the value of 0.986 corresponds
    #to a correction value of the equation.

    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    arcpy.SetProgressorLabel("Calculating Emissivity...")
    
    PV = Raster(Pv)
    Outfolder = os.path.join(workspace,"Emisivity.TIF")
    coc = 0.004
    cotn = 0.986
    Emis = sa.Float(coc * PV + cotn)
    
    Emis.save(Outfolder)
    return (Outfolder)

def LST(Bt, Emissi):
    
    #Calculate the Land Surface Temperature
    #LST = (BT / (1 + (0.00115 * BT / 1.4388) * Ln(ε)))

    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True

    arcpy.SetProgressorLabel("Calculating Land surface temperature for Band 10...")

    BT = Raster(Bt)
    Emissiv = (Emissi)
    Outfolder = os.path.join(workspace,"LST_B10.TIF")
    LST = arcpy.sa.Float(BT/(1 + (0.00115 * BT/ 1.4388) * Ln(Emissiv)))
    
    LST.save(Outfolder)
    return(Outfolder)

def LST11(Bt, Emissi):

    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True

    arcpy.SetProgressorLabel("Calculating Land surface temperature for Band 11...")

    BT = Raster(Bt)
    Emissiv = (Emissi)
    Outfolder = os.path.join(workspace,"LST_B11.TIF")
    LST = arcpy.sa.Float(BT/(1 + (0.00115 * BT/ 1.4388) * Ln(Emissiv)))
    
    LST.save(Outfolder)
    return (Outfolder)

def Average (B0, B1):

    # This section calculates the average of the two bands
    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True

    arcpy.SetProgressorLabel("Calculating Average Land surface temperature for Band 10 and 11...")

    B10 = Raster(B0)
    B11 = Raster(B1)
    
    Outfolder = os.path.join(workspace,"Average_BT.TIF")

    average = arcpy.sa.Float(B10 + B11)
    Average = arcpy.sa. Divide(average,2)

    Average.save(Outfolder)
    return(Outfolder)


def Radiance6(B6):

    #This section Calculate the radiance for Band_6. Unlike Band_10 and 11 in landsat_8
    # Band_6 does not require NDVI or Proportion of vegetation.
    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    arcpy.SetProgressorLabel("Radiance for Band_6...")
    
    Outfolder = os.path.join(workspace + "Radiance_B6.TIF")
    Radiance_Max = 15.303
    Radiance_Min = 1.235
    Quantile_Max = 255
    Quantile_Min = 1
    
    RMx = Radiance_Max 
    RMn = Radiance_Min
    QMx = Quantile_Max
    QMn = Quantile_Min
    

    TOA6 = arcpy.sa.Float((RMx - RMn)/(QMx - QMn) * (Raster(B6) - QMn) + RMn) 
    TOP6.save(Outfolder)
    return Outfolder


def LST6(TOA6):

    #This section calculates the Land surface temperature for band_6
    
    import os 
    from arcpy import *
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("spatial")

    arcpy.SetProgressorLabel("Calculating Land surface temperature for Band 6...")
    
    Outfolder = os.path.join(workspace + "LST_B6.TIF")
    
    K2 = 1260.56
    K1 = 607.76
    fh = 273.15
    
    Logs = arcpy.sa.ln((K1/TOA6) + 1)
    div = arcpy.sa.Float(K2/logs - fh)
    
    div.save(Outfolder)
    return Outfolder



#/////////////////////////////////////////////////////////Main//////////////////////////    



if average.lower() == "true":
    
    averages = Average((BT(TOA(B10))), (BT11(TOA11(B11))))
    LST(averages, (Emissiv (PV(NDVI(NIR,Red)))))

elif Landsat_7_Or_below.lower () == "true":

    LST(Radiance(B6))
    
else:

    LST((BT(TOA(B10))), (Emissiv (PV(NDVI(NIR,Red)))))
    LST11((BT11(TOA11(B11))), (Emissiv (PV(NDVI(NIR,Red)))))
    

#///////////////////////////////////////////////////Results//////////////////////
    
#This tool will produce:

    #Radiance

    #Brightness temperature

    #NDVI

    #Proportion of vegetation

    #Emissivity

    #Land Surface temperature

        









