import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
import SpectraByte as sb
import os
from tkinter import filedialog
from scipy.signal import find_peaks

# Create tkinter root window to query for data folder
Root = tk.Tk()
Root.withdraw()
Root.attributes('-topmost', True)
Parent = filedialog.askdirectory(title="Select a folder:")
if not Parent:
    print("No folder selected.")
else:
    print(f"Selected folder: {Parent}")
Root.destroy()
print("Parent", Parent)
os.chmod(Parent, 777)
print(os.listdir(Parent))

for Folder in os.listdir(Parent):
    print(os.listdir(Parent + '/' + Folder))
    print(os.getcwd())
    os.chdir(Parent)
    print("Folder", Folder)
    # Create file paths to data text files
    HeaderPath = Folder + '/Header.txt'
    XFilePath = Folder + '/X-Axis.txt'
    YFilePath = Folder + '/Y-Axis.txt'
    DataPath = Folder + '/Data.txt'
    DataFile = open(DataPath, "w").close()
    DataFile = open(DataPath, "w")

    # Extract Data from Header
    HeaderFile = open(HeaderPath)
    HeaderLines = [Line.rstrip('\n') for Line in HeaderFile]

    StepsX = (int)(HeaderLines[1][9:])
    StepsY = (int)(HeaderLines[2][9:])
    StepSize = (float)(HeaderLines[3][16:])

    ScanWidth = (int)(StepSize * StepsX)
    ScanHeight = (int)(StepSize * StepsY)
    SpectrumWidth = (int)(HeaderLines[4][16:])
    HeaderFile.close()

    # Set Laser Wavelength
    # FDU    = 487.79   nm
    # LEHIGH = 531.9585 nm
    LaserWavelength = 487.79

    # Array Initialization
    ScanData = np.zeros((StepsY, StepsX, SpectrumWidth))
    NormalizedScan = np.zeros((StepsY, StepsX, SpectrumWidth))
    CenterOfGravity = np.zeros((StepsY, StepsX))
    SpectrumMax = np.zeros((StepsY, StepsX))
    CrystalDetection = np.zeros((StepsY, StepsX))
    GlassDifference = np.zeros((StepsY, StepsX))
    CrystalDifference = np.zeros((StepsY, StepsX))
    Wavelengths = []
    Wavenumbers = []
    YData = []

    # Extract Wavelength Data 
    XFile = open(XFilePath)
    for i in range(SpectrumWidth):
        Wavelengths.append(float(XFile.readline()))
    XFile.close()

    # Find Wavenumber for Spectrum
    for i in range(len(Wavelengths)):
        Wavenumbers.append(10000000*(LaserWavelength-Wavelengths[i])/(LaserWavelength*Wavelengths[i]))

    # Find Intensity at Each Point on Each Spectrum
    YFile = open(YFilePath)
    for i in range(StepsX*StepsY*SpectrumWidth):
        YData.append((float)(YFile.readline()))
    YFile.close()

    # Raw Data Extraction, Normalization, Center of Gravity and Finding Spectrum Maxes
    for i in range(StepsY):
        for j in range(StepsX):
            Weight = 0
            Area = np.zeros(SpectrumWidth)
            for k in range(SpectrumWidth):
                ScanData[i][j][k] = YData[k+(i+j*StepsY)*SpectrumWidth]
                Area[k] = ScanData[i][j][k] * Wavelengths[k]
                if(ScanData[i][j][k] > SpectrumMax[i][j]):
                    SpectrumMax[i][j] = ScanData[i][j][k]
            for m in range(SpectrumWidth):
                NormalizedScan[i][j][m] = ScanData[i][j][m] / SpectrumMax[i][j]

    # Crystal Detection
    MaxDifference = 0
    CrystalLocation = np.zeros(2)
    for i in range(StepsY):
        for j in range(StepsX):
            for k in range(SpectrumWidth):
                GlassDifference[i][j] = GlassDifference[i][j] + abs(NormalizedScan[i][j][k] - NormalizedScan[0][0][k])
            if (GlassDifference[i][j] > MaxDifference):
                MaxDifference = GlassDifference[i][j]
                CrystalLocation[0] = i
                CrystalLocation[1] = j
    print(CrystalLocation[0], CrystalLocation[1])
    for i in range(StepsY):
        for j in range(StepsX):
            for k in range(SpectrumWidth):
                CrystalDifference[i][j] = CrystalDifference[i][j] + abs(NormalizedScan[i][j][k] - NormalizedScan[int(CrystalLocation[0])][int(CrystalLocation[1])][k])
            if(GlassDifference[i][j] > CrystalDifference[i][j]):
                CrystalDetection[i][j] = 1

    # Peak Finder on Specific Spectrum
    X = int(CrystalLocation[1])
    Y = int(CrystalLocation[0])
    Peaks, _ = find_peaks(ScanData[int(CrystalLocation[0]),int(CrystalLocation[1]),:], height = 500, distance = 10, prominence = 50)
    GraphPeaks = np.zeros((len(Peaks)))
    for i in range(len(Peaks)):
        GraphPeaks[i] = Wavelengths[Peaks[i]]

    # Determine which peak to look at
    PeakNumber = 0
    MaxPeak = 0
    for i in range(len(Peaks)):
        if(ScanData[Y,X,i] > MaxPeak):
            MaxPeak = ScanData[Y,X,i]
            PeakNumber = i

    # Center of Gravity for Each Peak
    CenterOfGravityMulti = []
    for i in range(len(GraphPeaks)):
        if(i == PeakNumber):
            CenterOfGravityMulti.append(sb.CenterOfGravity(ScanData, Peaks[i], 20, Wavelengths, StepsY, StepsX, SpectrumWidth))
        else:
            CenterOfGravityMulti.append(np.zeros((StepsY, StepsX)))

    # Intensity Colormap
    XAxis,YAxis = np.meshgrid(np.linspace(0, ScanWidth, StepsX), np.linspace(0, ScanHeight, StepsY))
    IntensityFigure = plt.figure(1)
    IntensityFigure.clf()
    Axes = plt.gca()
    # 586 for Pr Fluorescence at FDU, 359 for LaBGeO5 Raman at Lehigh
    PlotColor = Axes.pcolormesh(XAxis,YAxis,ScanData[:,:,586].squeeze(), cmap = plt.get_cmap('plasma'))
    ColorBar = plt.colorbar(PlotColor)
    ColorBar.ax.set_ylabel('Intensity (a.u.)')
    Axes.set_title('Flourescence Intensity')
    Axes.set_xlabel('Position (µm)')
    Axes.set_ylabel('Position (µm)')
    Axes.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    Axes.xaxis.set_label_position('top')
    plt.ylim(np.max(plt.ylim()), np.min(plt.ylim()))
    Axes.figure.set_size_inches(ScanWidth*.1,ScanHeight*.1) # Not actually correct
    plt.savefig(Folder + '/Intensity.png',bbox_inches='tight')

    # Specific Center of Gravity Colormap
    XAxis,YAxis = np.meshgrid(np.linspace(0, ScanWidth, StepsX), np.linspace(0, ScanHeight, StepsY))
    CenterOfGravityFigure = plt.figure(2)
    CenterOfGravityFigure.clf()
    Axes = plt.gca()
    PlotColor = Axes.pcolormesh(XAxis,YAxis,CenterOfGravityMulti[PeakNumber].squeeze(), cmap = plt.get_cmap('viridis'))
    ColorBar = plt.colorbar(PlotColor)
    ColorBar.ax.set_ylabel('Wavelength (nm)')
    Axes.set_title('Center of Gravity Shift')
    Axes.set_xlabel('Position (µm)')
    Axes.set_ylabel('Position (µm)')
    Axes.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    Axes.xaxis.set_label_position('top')
    plt.ylim(np.max(plt.ylim()), np.min(plt.ylim()))
    Axes.figure.set_size_inches(ScanWidth*.1,ScanHeight*.1) # Not actually correct
    plt.savefig(Folder + '/CenterOfGravity.png',bbox_inches='tight')

    # Crystal Detection Colormap
    XAxis,YAxis = np.meshgrid(np.linspace(0, ScanWidth, StepsX), np.linspace(0, ScanHeight, StepsY))
    CrystalDetectionFigure = plt.figure(3)
    CrystalDetectionFigure.clf()
    Axes = plt.gca()
    PlotColor = Axes.pcolormesh(XAxis,YAxis,CrystalDetection[:,:].squeeze(), cmap = plt.get_cmap('cividis'))
    ColorBar = plt.colorbar(PlotColor)
    ColorBar.ax.set_ylabel('Crystal (a.u.)')
    Axes.set_title('Crystal vs Glass')
    Axes.set_xlabel('Position (µm)')
    Axes.set_ylabel('Position (µm)')
    Axes.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    Axes.xaxis.set_label_position('top')
    plt.ylim(np.max(plt.ylim()), np.min(plt.ylim()))
    Axes.figure.set_size_inches(ScanWidth*.1,ScanHeight*.1) # Not actually correct
    plt.savefig(Folder + '/CrystalDetection.png',bbox_inches='tight')

    # Individual Spectrum Graph
    SpectrumPeaksFigure = plt.figure(4)
    SpectrumPeaksFigure.clf()
    plt.subplots(figsize=(6,4))
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Intensity (a.u.)')
    plt.plot(Wavelengths,ScanData[Y,X,:].squeeze(), color = 'tab:blue')
    plt.plot(GraphPeaks, ScanData[Y,X,Peaks], 'x')
    plt.savefig(Folder + '/'+str(X/2)+'x'+str(Y/2)+'Spectrum.png', bbox_inches='tight')

    # Compare Glass Shift to Crystal Shift
    GlassCount, CrystalCount = 0, 0
    GlassAverage, CrystalAverage = 0, 0
    for i in range(StepsY):
        for j in range(StepsX):
            if(CrystalDetection[i][j] == 0):
                GlassAverage = GlassAverage + CenterOfGravityMulti[PeakNumber][i,j]
                GlassCount = GlassCount + 1
            if(CrystalDetection[i][j] == 1):
                CrystalAverage = CrystalAverage + CenterOfGravityMulti[PeakNumber][i,j]
                CrystalCount = CrystalCount + 1
    if(CrystalCount > 1 and GlassCount > 1):
        CrystalAverage = CrystalAverage / CrystalCount
        GlassAverage = GlassAverage / GlassCount
    DataFile.write('Looking at peak ' + str(PeakNumber) + '\n')
    DataFile.write('Crystal Average: ' + str(CrystalAverage) + '\n')
    DataFile.write('Glass Average: ' + str(GlassAverage) + '\n')

    # Compare Inner Glass to Outer Glass
    InnerCrystalDetection = np.zeros((StepsY, StepsX))
    InnerCrystalCount, EdgeCrystalCount = 0,0
    InnerCrystalAverage, EdgeCrystalAverage = 0,0
    YMax, YMin, XMax = 0,0,0
    XMin = StepsX
    for i in range(StepsY):
        for j in range(StepsX):
            if(CrystalDetection[i][j] == 1):
                if(i > YMax):
                    YMax = i
                if(YMin == 0):
                    YMin = i
                if(j > XMax):
                    XMax = j
                if(j < XMin):
                    XMin = j
                
                if(i < StepsY - 1 and j < StepsX - 1):
                    if(CrystalDetection[i+1][j] == 0 or CrystalDetection[i-1][j] == 0 or CrystalDetection[i][j+1] == 0 or CrystalDetection[i][j-1] == 0 or CrystalDetection[i][j-1] == 0):
                        InnerCrystalDetection[i][j] = 1
                        EdgeCrystalCount = EdgeCrystalCount + 1
                        EdgeCrystalAverage = EdgeCrystalAverage + CenterOfGravityMulti[PeakNumber][i,j]
                    else:
                        InnerCrystalDetection[i][j] = 2
                        InnerCrystalCount = InnerCrystalCount + 1
                        InnerCrystalAverage = InnerCrystalAverage + CenterOfGravityMulti[PeakNumber][i,j]
    InnerCrystalAverage = InnerCrystalAverage / InnerCrystalCount
    EdgeCrystalAverage = EdgeCrystalAverage / EdgeCrystalCount
    DataFile.write('Inner Average ' + str(InnerCrystalAverage) + '\n')
    DataFile.write('Edge Average ' + str(EdgeCrystalAverage) + '\n')

    print('X Max:', XMax, 'X Min:', XMin, 'Y Max:', YMax, 'Y Min:', YMin)
    DataFile.write('Center to Edge Distance (X): ' + str(float((XMax - XMin)/4)) + '\n')
    DataFile.write('Center to Edge Distance (Y): ' + str(float((YMax - YMin)/4)) + '\n')
    CenterX, CenterY = int(XMin + ((XMax - XMin)/2)), int(YMin + ((YMax - YMin)/2))
    InnerCrystalDetection[CenterY, CenterX] = 3

    # Area of Crystal
    CrystalArea = 0
    for i in range(StepsY):
        for j in range(StepsX):
            if(CrystalDetection[i][j] == 1):
                CrystalArea = CrystalArea + 1
    CrystalArea = CrystalArea * 0.25
    DataFile.write('Crystal Area (µm²): ' + str(CrystalArea) + '\n')

    # Crystal Edge Detection Colormap
    XAxis,YAxis = np.meshgrid(np.linspace(0, ScanWidth, StepsX), np.linspace(0, ScanHeight, StepsY))
    CrystalDetectionFigure = plt.figure(3)
    CrystalDetectionFigure.clf()
    Axes = plt.gca()
    PlotColor = Axes.pcolormesh(XAxis,YAxis,InnerCrystalDetection[:,:].squeeze(), cmap = plt.get_cmap('cividis'))
    ColorBar = plt.colorbar(PlotColor)
    ColorBar.ax.set_ylabel('Crystal (a.u.)')
    Axes.set_title('Crystal vs Glass')
    Axes.set_xlabel('Position (µm)')
    Axes.set_ylabel('Position (µm)')
    Axes.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    Axes.xaxis.set_label_position('top')
    plt.ylim(np.max(plt.ylim()), np.min(plt.ylim()))
    Axes.figure.set_size_inches(ScanWidth*.1,ScanHeight*.1) # Not actually correct
    plt.savefig(Folder + '/InnerCrystalDetection.png',bbox_inches='tight')

    print("Done")