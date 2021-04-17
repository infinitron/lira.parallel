from collections import namedtuple
from typing import NewType

PointPhys=namedtuple('PointPhys','x,y')
PointCel=namedtuple('PointCel','ra,dec')
FluxObj=namedtuple('FluxObj','photon_flux,  monoenergy')

LiraInputConfig=NewType('LiraInputConfig',dict)

