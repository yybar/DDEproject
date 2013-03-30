#!usr/bin/env python
# 
# waveFrontWServer.py
# function that calls the ZEMAX file, extracts wavefront data
# and plots it out using pylab
#
# (c) Robert Upton
# January 10 2011
import os, sys, re
import numpy as num
import pylab as py
import pdb
import DDEServerClass as DDEServ                                            # import the toolbox as a class
import ImageScienceToolBox as ISTB
dataISTB = ISTB.ImageScienceClass()

def wavefrontCalc(name):
    '''This function calculates the required wavefront map from ZEMAX and then
    pads it according to padFac.  This is the kind of operation that is useful
    for PSF modeling, or OTF modeling of the optical system
    wfmap, rmsWFE = waveFrontCalc("fname", padFac)
    '''
    zmxName =  'z:\DDEProject\ZemaxModels\WOOLPERT/' + name + '.ZMX'            # name that is the filename of Zemax file
    txtName  = 'z:\DDEProject\data\WFmap.dat'
    cfgName  = 'z:\DDEProject\ZemaxModels\WOOLPERT/' + name + '.cfg'

    self = DDEServ.DDEToolBox(zmxName)                                          # open up the DDE server connection to ZEMAX using __init__(self,name)
    zCoeffUd = []
    h0Vec     = []
    dx  = 0.01
    dy  = 0.01
    dz  = 0.01
    dtx = 1e-02
    dty = 1e-02
    nZernike = 13

    for ii in num.arange(8):                                                    # generate 11x3 matrix
        self = self.zLoadFile()                                                 # loading the file containing
        self = self.zSetField(1,0,0)                                            # set the field prior to calculating wfm

# get the wfmap
# perform the modal decomposition
# perturb mirror though 3 degrees-of-freedom
        if ii == 0:
            dof = [0, 0, 0, 0, 0 ,0]
            surf = 12
        if ii == 1:
            dof = [dx, 0, 0, 0, 0, 0]
            surf = 16
        if ii == 2:
            dof = [0, dy, 0, 0, 0, 0]
            surf = 16
        if ii == 3:
            dof = [0, 0, dz, 0, 0, 0]
            surf = 16
        if ii == 4:
            dof = [0, 0, 0, dtx, 0, 0]
            surf = 16
        if ii == 5:
            dof = [0, 0, 0, 0, dty, 0]
            surf = 16
        if ii == 6:
            dof = [dx, 0, 0, 0, 0, 0]
            surf = 52
        if ii == 7:
            dof = [0, dy, 0, 0, 0, 0]
            surf = 52
        if ii > 3 and ii < 6:
            self.zOffAxisTilt(dof[3:6], [0,-14], -50.818, -1)
            vertexVec = self.vertexVec
            dof[0] = dof[0] + vertexVec[0]
            dof[1] = dof[1] + vertexVec[1]
            dof[2] = dof[2] + vertexVec[2]

        self = self.zRBChange(surf, dof)
# get the new wfmap
        self     = self.zWFmap(txtName,cfgName)      # calculate the wavefront map
        self.zWFsensorModes(nZernike)
        if ii == 0:
            h0Vec.append(self.coeffZ)                # as designed performance
            zCoeff0 = self.coeffZ
        if ii >= 1:
            ll = 0
            for kk in self.coeffZ:
            #    pdb.set_trace()
                delta = (kk-zCoeff0[ll])/num.max(dof)
                zCoeffUd.append(delta) # wavefront sensing response
                ll = ll + 1
                print ii

# shutdown the connection
    zCoeff    = num.array(zCoeffUd)
    zCoeffVec = num.reshape(zCoeff,(len(zCoeffUd)/ii,ii),order = 'c')

    self.zCoeffVec = zCoeffVec
    self.h0Vec     = num.transpose(h0Vec)
    self     = self.zDDEShutdown()                                              # shutdown the DDE connection
    py.axis('off')

    return self

def GLASanalysis(name,nZernike,diameter):
    '''This function calculates the required wavefront map from ZEMAX and then
    pads it according to padFac.  This is the kind of operation that is useful
    for PSF modeling, or OTF modeling of the optical system
    wfmap, rmsWFE = waveFrontCalc("fname", padFac)
    '''
    zmxName =  'z:\DDEProject\ZemaxModels\GLAS/' + name + '.ZMX'            # name that is the filename of Zemax file
    txtName  = 'z:\DDEProject\data\WFmap.dat'
    cfgName  = 'z:\DDEProject\ZemaxModels\GLAS/' + name + '.cfg'

    self = DDEServ.DDEToolBox(zmxName)                                      # open up the DDE server connection to ZEMAX using __init__(self,name)
    self = self.zLoadFile()                                                 # loading the file containing
    self = self.zSetField(1,0,0)
    self = self.zWFmap(txtName,cfgName)                                     # calculate the wavefront map
    self.zWFsensorModes(nZernike)
    
    wfmap = self.wfmap
    wfFFT = num.fft.fft2(wfmap)
    wfFFT = num.fft.fftshift(wfFFT)
    wfFFT = num.multiply((num.conj(wfFFT)),(wfFFT)).real
    self.wfFFT = wfFFT

    dataISTB.encEnergy(wfFFT,1)
    EE1 = dataISTB.encE
    dataISTB.encEnergy(wfFFT,0)
    EE2 = dataISTB.encE
    encE = EE1-EE2
    self.EE = encE

    #pdb.set_trace()
    dataISTB.structFcn(wfmap,diameter)
    Dfunc = dataISTB.Dfunc
    sep   = dataISTB.sep
    self.Dfunc = Dfunc
    self.sep     = sep

    return self

if __name__ == "__main__":
    """call the wavefrontCalc function
    """
    self =  wavefrontCalc('WOOLPERTreceiverReflectiveNegSecFreeSpaceTradTolAnalysisTemp')
    sVec = num.shape(self.zCoeffVec)
    fwrite  = open('z:\DDEProject\data\HmatrixM2LN.dat','w')
    for ii in num.arange(sVec[0]):
        for jj in num.arange(sVec[1]):
            fwrite.write('%8.3e  '%(self.zCoeffVec[ii][jj]))
        fwrite.write('\n')
    fwrite.close()

    fwrite = open('z:\DDEProject\data\h0matrixM2LN.dat','w')
    for ii in num.arange(len(self.h0Vec)):
        fwrite.write('%6.3e '%(self.h0Vec[ii]))
    fwrite.close()

    py.ion()
    py.figure(1)
    py.imshow(self.wfmap)
    py.colorbar()
    raw_input('%Press <return> to close')
    



