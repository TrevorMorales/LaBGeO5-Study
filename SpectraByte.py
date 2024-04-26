import numpy as np

def Normalize(Spectrum):
    Max = FindMax(Spectrum)
    for i in len(Spectrum):
          Spectrum[i] = Spectrum[i]/Max

def WavelengthToIndex(Wavelength, Wavelengths):
     for i in range(len(Wavelengths)):
          if((int)(Wavelengths[i]) == Wavelength):
               return i

def FindMax(Spectrum):
    Max = 0
    for i in range(len(Spectrum)):
        if(Spectrum[i] > Max):
            Max = Spectrum[i]
    return Max

def CenterOfGravity(ScanData, PeakIndex, IndexRange, Wavelengths, StepsY, StepsX, SpectrumWidth):
    StartIndex = (int)((PeakIndex)-(IndexRange / 2))
    CenterOfGravity = np.zeros((StepsY, StepsX))
    for i in range(StepsY):
        for j in range(StepsX):
            Weight = 0
            Area = np.zeros((SpectrumWidth))
            for k in range(SpectrumWidth):
                Area[k] = ScanData[i][j][k] * Wavelengths[k]
            for l in range(IndexRange):
                CenterOfGravity[i][j] = CenterOfGravity[i][j] + Area[l + StartIndex]
                Weight = Weight + ScanData[i][j][l + StartIndex]
            CenterOfGravity[i][j] = CenterOfGravity[i][j] / Weight
    return CenterOfGravity