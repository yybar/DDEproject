#!usr/bin/env python
#
# Script to demonstrate the DDE server/client relationship between Python and Zemax
# Use of win32ui and dde modules
# The server is established and named OptCClient.  The server is encapsulated in ddeServer.
# A conversation is created with ZEMAX and is assigned the variable ddeConv.
# The link is terminated through the method ddeServer.Shutdown()
#
# (c) Robert Upton
# January 09, 2012

import win32ui
import dde
import pdb
import numpy as num
import ImageScienceToolBox as ISTB
dataISTB = ISTB.ImageScienceClass()
global ddeServer
global ddeConv

class DDEToolBox():
# This is a class of methods that define how python will talk to ZEMAX
# Included in this class is the instance of the class that defines how
# the class is used, as well as loading and manipulating lens files

#-----------------------------------------------------------------------#
    def __init__(self,name):
        '''when class initialized introduce the variable self that will
           have a number of attributes appended to it.  These include:
           the name of the file to be perturbed   = self.name
           the global variables ddeServer and ddeConv
        '''
# define the name
        self.name = name
# establish the link and create a conversation
        ddeServer = dde.CreateServer()
        ddeServer.Create("OptCClient")
        ddeConv   = dde.CreateConversation(ddeServer)

        ddeConv.ConnectTo("ZEMAX", "System")
        if ddeConv.Connected() == 1:
            print "Connection established"

        self.ddeConv   = ddeConv
        self.ddeServer = ddeServer

#-----------------------------------------------------------------------#
    def zLoadFile(self):
        ''' Load a ZEMAX file retVal = zLoadFile(fname)'''

# load the file in question and don't forget to push the lens from
# the ZEMAX lens buffer
        strVal = 'LoadFile, ' + self.name

        self.ddeConv.Request(strVal)
        strVal = 'PushLens'
        retVal = self.ddeConv.Request(strVal)
        self.retVal = retVal
        return self

 #-----------------------------------------------------------------------#
    def zRBChange(self, surf, dof):
        ''' introduce a 6 component dof vector that changes the rigid body
            state of a surface
         '''
# load the file in question and don't forget to push the lens from
# the ZEMAX lens buffer

# set the x, y, and z surface changes in mm
# set the x-tilt, y-tilt, and z-tilt in degrees
        kk = 0
        for ii in num.arange(6):
            if ii <> 2:                                   # xdec, ydec, xtilt, ytilt, ztilt
                kk = kk+1                                 # dof index
                strVal = 'SetSurfaceParameter, ' \
                       + str(surf) + ',' + str(kk) + ',' + str(dof[ii])
            elif ii == 2:                                 # z thickness perturbation
                 strVal = 'SetSurfaceData, '\
                        + str(surf) + ',' + str(3) + ',' + str(dof[2])
            self.ddeConv.Request(strVal)
        print dof
        strVal = 'PushLens'   # do I need to do this for each element in dof?
        retVal = self.ddeConv.Request(strVal)
        return self

 #-----------------------------------------------------------------------#
    def zOffAxisTilt(self, tiltVec, apLoc, Rb, k):
        '''calculate the required vertex motion of the off-axis conic
        '''
# alias cos, sin
        cos = num.cos
        sin = num.sin
        pi  = num.pi
# determine position of segment in 3D space
        xPos = apLoc[0]
        yPos = apLoc[1]
        zPos=1/Rb*(xPos**2+yPos**2)/\
        (1+num.sqrt(1-(1+k)*(xPos**2+yPos**2)/Rb**2))                           # Z sag position

# set up rotation matrices %
        tx=pi/180*tiltVec[0] 
        ty=pi/180*tiltVec[1]
        tz=pi/180*tiltVec[2]       # deg to rad
        Rlx=[[1, 0, 0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]]          # x rotation matrix
        Rly=[[cos(ty), 0, sin(ty)], [0, 1, 0], [-sin(ty), 0, cos(ty)]]          # y rotation matrix
        Rlz=[[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0, 0, 1]]          # z rotation matrix

        Rlyz = num.dot(Rly,Rlz)
        Rl   = num.dot(Rlx,Rlyz)                                                         # compound rotation matrix
        posVec = [xPos, yPos, zPos]
        posVec = num.transpose(posVec)

        vertexVec=posVec-num.dot(Rl,posVec)                                                   # vertex dispacements
        self.vertexVec = vertexVec

#-----------------------------------------------------------------------#
    def zSetField(self, fieldNum, xF, yF):
        ''' Change the field'''
   
# set up the change field variable    
        strVal = 'SetField,' + str(fieldNum) + ', ' + str(xF) + ', ' + str(yF) + ', ' + str(1)
    
        self.ddeConv.Request(strVal)
        strVal = 'PushLens'
        retVal = self.ddeConv.Request(strVal)
        self.retVal = retVal
        return self
 
#-----------------------------------------------------------------------#
    def zWFmap(self, txtFname, setFname):
        '''get the wavefront map data
           retVal=DDEServerToolBox.zWFmap(ddeConv,"dataTxtFname","configFileName")
        '''
# form the GetTextFile porgammable extension keyword
        strVal = 'GetTextFile, ' + txtFname + ', Wfm, ' + setFname + ', ' + str(0)
        retVal = self.ddeConv.Request(strVal)

        fid  = open(txtFname,'r')
        data = fid.readlines()
        fid.close()
# make data into an array
        wfmap = num.zeros((len(data[16:]),len(data[16:])))
        ii    = 0
        for lines in data[16:]:
            test = lines.split()
            for jj in num.arange(len(wfmap)-1):
                wfmap[ii,jj] = test[jj]
            ii = ii + 1

        self.wfmap = wfmap

        return self

#-----------------------------------------------------------------------#
    def zWFsensorModes(self, nZernike):
        '''Function that generates the wavefront sensor modes that 
           are nZernike modes excluding piston, tip, and tilt
        '''
        wfmap = self.wfmap
        nsamp = len(wfmap)
        dataISTB.pol2cartMat(nsamp,1)                # generate the mask for the Zdecomp
        mask = dataISTB.maskPol
        dataISTB.waveDecomp(wfmap,nZernike,1,nsamp)  # perform the Zernike decomposition
        coeffZ = dataISTB.coeffZ
        self.coeffZ = coeffZ

 #-----------------------------------------------------------------------#
    def zDeformAsphere(self, asphereCoeff, surf):
        '''Function that introduces deformations to a surface as a general asphere.
           ZEMAX has to already be set-up with an aspheric surface
        '''

# set the surface changes due to asphere in mm
        kk = 0
        for ii in num.arange(len(asphereCoeff)):
            kk = kk+1                                 # dof index
            strVal = 'SetSurfaceParameter, ' \
                   + str(surf) + ',' + str(kk) + ',' + str(asphereCoeff[ii])
            self.ddeConv.Request(strVal)

        strVal = 'PushLens'   # do I need to do this for each element in dof?
        retVal = self.ddeConv.Request(strVal)
        return self

#-----------------------------------------------------------------------#
    def zDDEShutdown(self):
        '''terminate the link'''
        retVal = self.ddeServer.Shutdown()
        self.retVal = retVal
        return self
