'''
Liam O'Toole
binaryrk.py
Simulate a binary star system
23 October, 2014
'''


#from __future__ import division
from visual import *

#initial conditions
G = 6.67e-11
AU = 1.5e11
year= 365.25*24*60*60

#create stars
starA = sphere(pos=(0,0,0), radius=7E9, mass=2.0e30, color=color.red, text = 'A')
starB = sphere(pos=(AU,0,0), radius=3.5E9, mass=2.0e30, color=color.yellow, text = 'B')
CM = sphere(pos = (starA.pos*starA.mass + starB.pos*starB.mass)/(starA.mass+starB.mass), radius=1.0e9, color=color.green)

#create coordinate axes
L = 2e11
xaxis = curve(pos=[(0,0,0),(L,0,0)], color=(0.5,0.5,0.5))
yaxis = curve(pos=[(0,0,0),(0,L,0)], color=(0.5,0.5,0.5))
zaxis = curve(pos=[(0,0,0),(0,0,L)], color=(0.5,0.5,0.5))

#initial conditions
G = 6.67e-11
AU = 1.5e11
year= 365.25*24*60*60

starA.vel = vector(0,1.4*pi*AU/year,pi*AU/year)
starB.vel = vector(0,-1.4*pi*AU/year,-pi*AU/year)
dt = 1e4
t = 0
Rmax = 0
Rmin = 5*AU

starA.trail = curve(color=starA.color)
starB.trail = curve(color=starB.color)
CM.vel = vector(0,0,0)

ecc = 0
ecc_s= str("%.3f"%ecc)
ecc_label = label(pos = starA.pos, text='ecc:'+ecc_s, xoffset=20, yoffset=20)


#define the gravitational accerleration as a function
def acc(starpos):
    return -G*starother.mass*(starpos - starother.pos)/r**3
def rk4(star):
    k1v = acc(star.pos)*dt
    k1x = star.vel*dt

    k2v = acc(star.pos + k1x/2.0)*dt
    k2x = (star.vel + k1v/2.0)*dt

    k3v = acc(star.pos + k2x/2.0)*dt
    k3x = (star.vel + k2v/2.0)*dt

    k4v = acc(star.pos + k3x)*dt
    k4x = (star.vel + k3v)*dt

    star.vel += (k1v + 2.0*k2v + 2.0*k3v +k4v)/6.0
    star.pos += (k1x + 2.0*k2x + 2.0*k3x +k4x)/6.0


#simulation
while True:
    if scene.kb.keys:
        s=scene.kb.getkey()
        if s== 'b':
            starB.vel *=1.1
        if s == 'r':
            starB.vel *= 0.9
        if s == 'a':
            a=(Rmax+Rmin)/2.0
            print "Rmax= ", Rmax, "\t Rmin = ", Rmin,"\t","a= ",(Rmax+Rmin)/2.0
            T = 2*pi*sqrt((1/(G*starA.mass+G*starB.mass))*a**3)
            print "T=", T
    rate(100)
    r = mag(starA.pos - starB.pos)
    for i in range(0,2):
        if i == 1:
            r = mag(starA.pos - starB.pos)
            starother=starB
            rk4(starA)
        else:
            r = mag(starA.pos - starB.pos)
            starother = starA
            rk4(starB)
    starA.trail.append(pos=starA.pos)
    starB.trail.append(pos=starB.pos)

    CM.pos = (starA.pos*starA.mass + starB.pos*starB.mass)/(starA.mass+starB.mass)
    M = starA.mass + starB.mass
    mu = starA.mass*starB.mass/(starA.mass+starB.mass)
    L = mu*cross(starB.pos - starA.pos, starB.vel-starA.vel)
    E = 0.5*mu*mag2(starB.vel-starA.vel) - G*M*mu/r
    ecc = sqrt(1.0+2*mag2(L)*E/((G*M)**2*mu**3))
    ecc_s= str("%.3f" %ecc)
    ecc_label.pos = CM.pos
    ecc_label.text = 'ecc: '+ ecc_s
    R = mag(starA.pos-starB.pos)
    if R>Rmax:
        Rmax = R
    if R<Rmin:
        Rmin = R
    t+=dt








    
