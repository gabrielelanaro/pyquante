# mobimol v0.6 for s60, a simple molecule viewer, (c) V.Ganesh, GPL
from math import sin,cos,acos,pi,sqrt
import appuifw,e32,sysinfo,string,dir_iter
from sysinfo import *
from appuifw import *
from e32 import *
from key_codes import *
from os.path import *

radii = {"H":0.23, "C":0.77, "N":0.75, "O":0.73, \
         "S":1.02, "P":1.06, "NA":1.54, "CL":0.99, "CA":1.74}
atmcol = {"H":(0,0,128), "C":(0,0,0), "N":(0,0,255), "O":(255,0,0,), \
          "S":(255,255,0), "P":(128,128,0), "NA":(0,234,58), "CL":(0,232,232), "CA":(230,230,230)}

class mat3d:
     def __init__(self):
         self.xx = self.xy = self.xz = self.xo = 0.0
         self.yx = self.yy = self.yz = self.yo = 0.0
         self.zx = self.zy = self.zz = self.zo = 0.0
         
     def unit(self):
         self.xo = 0; self.xx = 1; self.xy = 0; self.xz = 0
         self.yo = 0; self.yx = 0; self.yy = 1; self.yz = 0
         self.zo = 0; self.zx = 0; self.zy = 0; self.zz = 1

     def mult(self, rhs):
         lxx = self.xx * rhs.xx + self.yx * rhs.xy + self.zx * rhs.xz
         lxy = self.xy * rhs.xx + self.yy * rhs.xy + self.zy * rhs.xz
         lxz = self.xz * rhs.xx + self.yz * rhs.xy + self.zz * rhs.xz
         lxo = self.xo * rhs.xx + self.yo * rhs.xy + self.zo * rhs.xz + rhs.xo
         
         lyx = self.xx * rhs.yx + self.yx * rhs.yy + self.zx * rhs.yz
         lyy = self.xy * rhs.yx + self.yy * rhs.yy + self.zy * rhs.yz
         lyz = self.xz * rhs.yx + self.yz * rhs.yy + self.zz * rhs.yz
         lyo = self.xo * rhs.yx + self.yo * rhs.yy + self.zo * rhs.yz + rhs.yo
         
         lzx = self.xx * rhs.zx + self.yx * rhs.zy + self.zx * rhs.zz
         lzy = self.xy * rhs.zx + self.yy * rhs.zy + self.zy * rhs.zz
         lzz = self.xz * rhs.zx + self.yz * rhs.zy + self.zz * rhs.zz
         lzo = self.xo * rhs.zx + self.yo * rhs.zy + self.zo * rhs.zz + rhs.zo
         
         self.xx = lxx; self.xy = lxy; self.xz = lxz; self.xo = lxo
         self.yx = lyx; self.yy = lyy; self.yz = lyz; self.yo = lyo
         self.zx = lzx; self.zy = lzy; self.zz = lzz; self.zo = lzo

     def scale(self, xf, yf, zf):
         self.xx *= xf; self.xy *= xf; self.xz *= xf; self.xo *= xf
         self.yx *= yf; self.yy *= yf; self.yz *= yf; self.yo *= yf
         self.zx *= zf; self.zy *= zf; self.zz *= zf; self.zo *= zf

     def translate(self, x, y, z):
         self.xo += x; self.yo += y; self.zo += z

     def xrot(self, theta):
         theta *= (pi / 180)
         ct = cos(theta)
         st = sin(theta)

         Nyx = self.yx * ct + self.zx * st
         Nyy = self.yy * ct + self.zy * st
         Nyz = self.yz * ct + self.zz * st
         Nyo = self.yo * ct + self.zo * st

         Nzx = self.zx * ct - self.yx * st
         Nzy = self.zy * ct - self.yy * st
         Nzz = self.zz * ct - self.yz * st
         Nzo = self.zo * ct - self.yo * st
         
         self.yo = Nyo; self.yx = Nyx; self.yy = Nyy; self.yz = Nyz
         self.zo = Nzo; self.zx = Nzx; self.zy = Nzy; self.zz = Nzz

     def yrot(self, theta):
         theta *= (pi / 180)
         ct = cos(theta)
         st = sin(theta)

         Nxx = self.xx * ct + self.zx * st
         Nxy = self.xy * ct + self.zy * st
         Nxz = self.xz * ct + self.zz * st
         Nxo = self.xo * ct + self.zo * st
         
         Nzx = self.zx * ct - self.xx * st
         Nzy = self.zy * ct - self.xy * st
         Nzz = self.zz * ct - self.xz * st
         Nzo = self.zo * ct - self.xo * st
         
         self.xo = Nxo; self.xx = Nxx; self.xy = Nxy; self.xz = Nxz
         self.zo = Nzo; self.zx = Nzx; self.zy = Nzy; self.zz = Nzz
            
     def zrot(self, theta):
         theta *= (pi / 180)
         ct = cos(theta)
         st = sin(theta)

         Nyx = self.yx * ct + self.xx * st
         Nyy = self.yy * ct + self.xy * st
         Nyz = self.yz * ct + self.xz * st
         Nyo = self.yo * ct + self.xo * st
         
         Nxx = self.xx * ct - self.yx * st
         Nxy = self.xy * ct - self.yy * st
         Nxz = self.xz * ct - self.yz * st
         Nxo = self.xo * ct - self.yo * st
         
         self.yo = Nyo; self.yx = Nyx; self.yy = Nyy; self.yz = Nyz
         self.xo = Nxo; self.xx = Nxx; self.xy = Nxy; self.xz = Nxz

     def transform(self, atms, tatms):
         lxx = self.xx; lxy = self.xy; lxz = self.xz; lxo = self.xo
         lyx = self.yx; lyy = self.yy; lyz = self.yz; lyo = self.yo
         lzx = self.zx; lzy = self.zy; lzz = self.zz; lzo = self.zo
         
         for atm in atms:
             x = atm[1][0]; y = atm[1][1]; z = atm[1][2]
             tatms.append((atm[0], \
                                ((x * lxx + y * lxy + z * lxz + lxo), \
                                 (x * lyx + y * lyy + z * lyz + lyo), \
                                 (x * lzx + y * lzy + z * lzz + lzo)), \
                           atm[2]))

class filechooser:
     def __init__(self):
         self.current_path = "<root>"

     def show(self, filters=None):
         filelist = e32.drive_list()
         while (1):
            ll = []
            for file in filelist:
               if (filters != None):
                  if isfile(self.current_path + "\\" + file):
                     matched = 0; filext = splitext(file)[1]
                     for filter in filters:
                         if (filext == filter):
                            matched = 1; break
                     if (matched == 0): continue
               ll.append(unicode(file))
                 
            ind = selection_list(ll)
             
            if ind is not None:
              if (self.current_path == "<root>"):
                  self.current_path = ll[ind]
                  filelist = os.listdir(self.current_path)
                  continue
              if not isfile(self.current_path + "\\" + ll[ind]):
                  self.current_path = self.current_path + "\\" + ll[ind]
                  filelist = os.listdir(self.current_path)
                  continue
              return self.current_path + "\\" + ll[ind]
            else:
              self.current_path = "<root>"
              return None
                           
class molecule:
     def __init__(self):
         self.noOfAtoms = 0
         self.title = ""
         self.atoms = []
         
     def findBB(self):
         self.min = [0.0, 0.0, 0.0]
         self.max = [0.0, 0.0, 0.0]

         self.min[0] = self.max[0] = self.atoms[0][1][0]
         self.min[1] = self.max[0] = self.atoms[0][1][1]
         self.min[2] = self.max[0] = self.atoms[0][1][2]

         for i in range(0, self.noOfAtoms):
             atmxyz = self.atoms[i][1]

             if (atmxyz[0] < self.min[0]): self.min[0] = atmxyz[0]
             if (atmxyz[1] < self.min[1]): self.min[1] = atmxyz[1]
             if (atmxyz[2] < self.min[2]): self.min[2] = atmxyz[2]
             
             if (atmxyz[0] > self.max[0]): self.max[0] = atmxyz[0]
             if (atmxyz[1] > self.max[1]): self.max[1] = atmxyz[1]
             if (atmxyz[2] > self.max[2]): self.max[2] = atmxyz[2]

class mobimol:
     def __init__(self):
         self.script_lock = e32.Ao_lock()

     def run(self):
         self.displayModel = 0
         self.mat = mat3d(); self.amat = mat3d(); self.tmat = mat3d()
         self.mat.unit(); self.amat.unit(); self.tmat.unit()
         self.scalefudge = 1.0
         self.xrd = self.zrd = 0.0
         self.symbols = self.ids = 0;
         self.transformMode = 'rotate'
         self.mol = None
         
         self.canvas = appuifw.Canvas(redraw_callback=self.paintStuff)
         self.old_body = appuifw.app.body
         self.old_screen = appuifw.app.screen
         self.old_title = appuifw.app.title
         appuifw.app.exit_key_handler = self.do_exit
         self.initUi()
         self.script_lock.wait()
         appuifw.app.title = self.old_title

     def initUi(self):
         appuifw.app.title = u"mobimol v0.6, (c)V.Ganesh"
         appuifw.app.menu = [(u"Open", self.do_open),
                             (u"Transform mode",
                              ((u"Rotate", self.do_rotate_mode),
                               (u"Zoom", self.do_zoom_mode))),
                             (u"Label",
                              ((u"Symbols", self.do_label_symbols),
                               (u"IDs", self.do_label_ids))),
                             (u"Query",
                              ((u"Length", self.do_query_length),
                               (u"Angle", self.do_query_angle),
                               (u"Dihedral", self.do_query_dihedral))),
                             (u"Display model",
                              ((u"Line", lambda:self.do_model(0)),
                               (u"Line and Circle", lambda:self.do_model(1)))),
                             (u"Screen mode",
                              ((u"Full screen", self.do_full_screen),
                               (u"Normal screen", self.do_normal_screen))),
                             (u"About", self.do_about)
                             ]
         appuifw.app.exit_key_handler = self.exit_key_handler
         self.canvas.bind(EKeyRightArrow,lambda:self.add_zrd(-2))
         self.canvas.bind(EKeyLeftArrow,lambda:self.add_zrd(2))
         self.canvas.bind(EKeyUpArrow,lambda:self.add_xrd(-2))
         self.canvas.bind(EKeyDownArrow,lambda:self.add_xrd(2))
         appuifw.app.body = self.canvas

     def do_about(self):
         note(appuifw.app.title, 'info')

     def do_label_symbols(self):
         self.symbols = (not self.symbols)
         self.paintStuff()

     def do_label_ids(self):
         self.ids = (not self.ids)
         self.paintStuff()

     def do_rotate_mode(self):
         self.transformMode = 'rotate'         
         self.paintStuff()
         
     def do_zoom_mode(self):
         self.transformMode = 'zoom'
         self.paintStuff()
         
     def do_model(self, model):
         self.displayModel = model
         self.paintStuff()
         
     def acceptNumbers(self, label, items):
         itmTxt = u""
         for i in range(0, items-1): itmTxt = itmTxt + repr(i) + ","
         
         itmTxt = itmTxt + repr(i+1)
         nos = query(label, 'text', itmTxt)         

         if (nos == None):
            note(u"Cancelled", 'error')
            return None

         nos = string.split(nos, ',')
         if (len(nos) != items):
            note(u"Inadequate data", 'error')
            return None

         try:
            nos = map(string.atoi, nos)
            return nos
         except:
            note(u"Wrong data", 'error')
            return None
       
     def do_query_length(self):
         atms = self.acceptNumbers(u"Enter two atom indices:", 2)
         if (atms == None): return
         try:
             a1=self.mol.atoms[atms[0]][1]; a2=self.mol.atoms[atms[1]][1]
         except:
             note(u"Invalid indices", 'error')
             return 
         x = a1[0]-a2[0]; y = a1[1]-a2[1]; z = a1[2]-a2[2]
         dist = sqrt(x*x+y*y+z*z)
         note(u"Distance (" + repr(atms[0]) + ", " + repr(atms[1]) + "): " \
              + str(round(dist, 2)), 'info')

     def findAngle(self, v1, v2):
         v1mod = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])
         v2mod = sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])
         v1dv2 = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
         ang = (v1dv2/(v1mod*v2mod))
         if(round(ang,1)==1.0): return 0.0
         return acos(ang)

     def cross(self, v1, v2):
         return ((v1[1]*v2[2]-v1[2]*v2[1]), (v1[2]*v2[0]-v1[0]*v2[2]), (v1[0]*v2[1]-v1[1]*v2[0]))

     def normalize(self, v):
         vmod = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
         return (v[0]/vmod, v[1]/vmod, v[2]/vmod)
    
     def do_query_angle(self):
         atms = self.acceptNumbers(u"Enter three atom indices:", 3)
         if (atms == None): return
         try:
              a1=self.mol.atoms[atms[0]][1]; a2=self.mol.atoms[atms[1]][1]; a3=self.mol.atoms[atms[2]][1]
              v1 = (a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2])
              v2 = (a2[0]-a3[0], a2[1]-a3[1], a2[2]-a3[2])
              ang = self.findAngle(v1, v2)*180/pi
         except:
              note(u"Invalid indices", 'error')
              return 
         note(u"Angle (" + repr(atms[0]) + ", " + repr(atms[1]) + ", " + repr(atms[2]) + "): " \
              + str(round(ang, 2)), 'info')

     def do_query_dihedral(self):
         atms = self.acceptNumbers(u"Enter four atom indices:", 4)
         if (atms == None): return
         try:
              a1=self.mol.atoms[atms[0]][1]; a2=self.mol.atoms[atms[1]][1]
              a3=self.mol.atoms[atms[2]][1]; a4=self.mol.atoms[atms[3]][1]
              v1 = (a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2])
              v2 = (a2[0]-a3[0], a2[1]-a3[1], a2[2]-a3[2])
              n1 = self.normalize(self.cross(v1, v2))
              v1 = (a3[0]-a2[0], a3[1]-a2[1], a3[2]-a2[2])
              v2 = (a3[0]-a4[0], a3[1]-a4[1], a3[2]-a4[2])
              n2 = self.normalize(self.cross(v1, v2))
              dihed = self.findAngle(n1, n2)*180/pi
         except:
              note(u"Invalid indices", 'error')
              return 
         note(u"Dihedral (" + repr(atms[0]) + ", " + repr(atms[1]) + ", " + repr(atms[2]) \
              + ", " + repr(atms[3]) + "): " + str(round(dihed, 2)), 'info')
         
     def add_xrd(self, d):
         if (self.transformMode == 'rotate'):
            self.xrd = d
            self.tmat.unit()
            self.tmat.xrot(self.xrd)
            self.tmat.yrot(self.xrd)
            self.amat.mult(self.tmat)
         elif (self.transformMode == 'zoom'):
            self.scalefudge = self.scalefudge + (1.0/(d*10))

         self.paintStuff()
         
     def add_zrd(self, d):
         if (self.transformMode == 'rotate'):
            self.zrd = d
            self.tmat.unit()
            self.tmat.zrot(self.zrd)
            self.amat.mult(self.tmat)
         elif (self.transformMode == 'zoom'):
            self.scalefudge = self.scalefudge + (1.0/(d*10))

         self.paintStuff()
              
     def do_open(self):
         self.filename = filechooser().show([".xyz"])

         if (self.filename == None): return
         self.open()
         
     def open(self):
         global radii
         xyzf = open(self.filename, "r")
         lines = xyzf.readlines()
         xyzf.close()

         self.mol = molecule()
         self.mol.noOfAtoms = string.atoi(lines[0])
         self.mol.title = string.strip(lines[1])
         for i in range(2, len(lines)):
             words = string.split(lines[i])
             self.mol.atoms.append((string.upper(words[0]), map(string.atof, words[1:]), []))

         self.mol.findBB()
         
         for i in range(0, self.mol.noOfAtoms):
             atm1 = self.mol.atoms[i]
             for j in range(0, i):
                 atm2 = self.mol.atoms[j] 
                 x = atm1[1][0]-atm2[1][0]
                 y = atm1[1][1]-atm2[1][1]
                 z = atm1[1][2]-atm2[1][2]
                 dist = sqrt(x*x+y*y+z*z)
                 radsum = radii[atm1[0]] + radii[atm2[0]]
                 if (((radsum - 0.4) < dist) and (dist < (radsum + 0.4))):
                     atm1[2].append(j)
                     atm2[2].append(i)
                     
         lines = None
         
     def paintStuff(self, rect=None):
         global atmcol
         if (self.mol == None):
              self.canvas.clear((255,255,255))
              self.canvas.text((1,14),u"Options->Open to get started",0x000000)
              return
              
         xw = self.mol.max[0] - self.mol.min[0]
         yw = self.mol.max[1] - self.mol.min[1]
         zw = self.mol.max[2] - self.mol.min[2]
              
         if (yw > xw): xw = yw
         if (zw > xw): xw = zw
         
         f1 = self.canvas.size[0] / xw; f2 = self.canvas.size[1] / xw;
         
         if (f1 < f2): f2 = f1
         xfac = 0.7 * f2 * self.scalefudge;                
         
         self.mat.unit()
         self.mat.translate(-(self.mol.min[0] + self.mol.max[0]) / 2.0,  \
                            -(self.mol.min[1] + self.mol.max[1]) / 2.0,  \
                            -(self.mol.min[2] + self.mol.max[2]) / 2.0);
         self.mat.mult(self.amat);
              
         self.mat.scale(xfac, -xfac, 16.0 * xfac / self.canvas.size[0]);
         self.mat.translate(self.canvas.size[0] / 2.0, self.canvas.size[1] / 2.0, 8.0);
         
         tatoms = []
         self.mat.transform(self.mol.atoms, tatoms);
         
         self.canvas.clear((255,255,255))
         self.canvas.text((10,14),u""+self.mol.title+": [" + self.transformMode + " mode]",0x000000)
         for i in range(0, self.mol.noOfAtoms):
              atm1 = tatoms[i]
              col1 = atmcol[atm1[0]]
              if (self.displayModel == 1):
                   r = 4
                   self.canvas.arc(((atm1[1][0]-r,atm1[1][1]-r),(atm1[1][0]+r,atm1[1][1]+r)),-pi,pi,outline=col1)
              for j in range(0, len(atm1[2])):
                   atm2 = tatoms[atm1[2][j]]
                   col2 = atmcol[atm2[0]]
                   mid = ((atm1[1][0]+atm2[1][0])/2, (atm1[1][1]+atm2[1][1])/2)
                   self.canvas.line(((atm1[1][0],atm1[1][1]),mid), col1)
                   self.canvas.line(((atm2[1][0],atm2[1][1]),mid), col2)
                   
                   if (self.symbols and self.ids):
                        self.canvas.text((atm1[1][0],atm1[1][1]),u" "+atm1[0]+repr(i),0x000000)
                   elif (self.symbols):
                        self.canvas.text((atm1[1][0],atm1[1][1]),u" "+atm1[0],0x000000)
                   elif (self.ids):
                        self.canvas.text((atm1[1][0],atm1[1][1]),u" "+repr(i),0x000000)                   
              
     def do_full_screen(self):
         appuifw.app.screen='full'

     def do_normal_screen(self):
         appuifw.app.screen='normal'
         
     def do_exit(self):
         self.exit_key_handler()

     def exit_key_handler(self):
         appuifw.app.exit_key_handler = None
         appuifw.app.body = self.old_body
         appuifw.app.screen = self.old_screen
         self.canvas.redraw_callback = None
         self.script_lock.signal()
         self.canvas = None

if __name__ == '__main__':
     mobimol().run()
