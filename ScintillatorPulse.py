"""
 Class for a detector pulse (Voltage as function of time)
 Init with "decaytime", "risetime", "amplitude constant", "baseline",
 and a time offset.
 """
import math
import numpy as np

class ScintillatorPulse:
  decayt   = 0.0 #Decay time
  riset    = 0.0 #Rise time
  aConst   = 0.0 #Area constant
  baseline = 0.0 #Baseline value
  tOffset  = 0.0 #time offset

  def __init__(self, decayt, riset, aConst, baseline, tOffset):
    self.decayt = decayt
    self.decayc = 1.0/decayt
    self.riset = riset
    self.aConst = aConst
    self.baseline = baseline
    self.tOffset = tOffset
    self.limitL  = tOffset - 5*self.riset
    self.limitR  = tOffset + 100*self.decayt

  def __pulseShape(self, xt):
    """
     A realistic scintillator pulse modelled on a convolution
     of a gaussian and an exponential decay starting at time zero.
    """
    return 0.5*self.aConst*self.riset*self.decayc\
         *math.exp((self.decayc)*(self.riset**2)*0.5)\
         *math.exp(-self.decayc*(xt-self.tOffset))\
         *(1-math.erf(((self.decayc*self.riset**2)-xt+self.tOffset)\
         /(self.riset*math.sqrt(2))))

  def normalise (self, norm=1.0):
    self.aConst = self.aConst/self.getArea()*norm

  def getMax (self):
    step = (self.limitR - self.limitL)/10000
    max = 0.0
    max_t = 0.0
    for t in np.arange(self.limitL,self.limitR, step):
      value = self.__pulseShape(t)
      if max > value:
        return [max_t, max]
      else:
        max = value
        max_t = t

  def getRealRisetime(self):
    peak = self.getMax()
    step = (peak[0] - self.limitL)/10000
    havet10 = 0
    t10 = 0.0
    for t in np.arange(self.limitL, peak[0],step):
      value = self.__pulseShape(t)
      if (havet10 != 0 and value > peak[1]*0.1) :
        havet10 = 1
        t10 = t
      elif (value > peak[1]*0.9) :
        return (t - t10)

   def getRealDecayTime(self):
     peak = self.getMax()
     step = (self.limitR - peak[0])/10000
     havet90 = 0
     t90 = 0.0
     for t in np.arange(peak[0], self.limitR, step):
       value = self.__pulseShape(t)
       if (havet90 != 0 and value < peak[1]*0.9) :
         havet90 = 1
         t90 = t
       elif (value < peak[1]*0.1) :
         return (t - t90)

  def getArea(self, tLow=0.0, tHigh=0.0, nSamples=0):
    area = 0.0
    # Set some reasonable default values if no ranges
    # and/or granularity is given.
    if (tLow==0.0 and tHigh == 0.0):
      tLow = self.limitL
      tHigh = self.limitR
    if (nSamples > 0):
      step = (tHigh - tLow)/nSamples
    else:
      step = (tHigh-tLow)/10000
    # Use average of "left side" and "right side"
    # integral
    for t in np.arange(tLow, tHigh-step, step):
      area += self.__pulseShape(t)*step
    for t in np.arange(tLow+step, tHigh, step):
      area += self.__pulseShape(t)*step
    return area/2
