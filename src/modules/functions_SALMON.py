import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_crystal_structure_SALMON(fname):
    f = open(fname,'r')
    lines = f.readlines()
    f.close()
    Nlen = len(lines)
    print(Nlen, ': No. of lines of the file, ',fname)
    al = np.zeros(3,dtype='float64')
    for i in lines:
        index = i.find('natom')
        if (index > 0):
            natom = int(i.split('=')[1])
            print(natom)
        index = i.find(' al ')
        if (index > 0):
            temp = i.split('=')[1]
            temp = temp.replace(',', '')
            temp = temp.replace('d','e')
            al[0] = float(temp.split()[0])
            al[1] = float(temp.split()[1])
            al[2] = float(temp.split()[2])
    atomdict = {}
    xred = np.zeros(3,dtype='float64')
    xcart = np.zeros(3,dtype='float64')
    for i in range(Nlen):
        index = lines[i].find('atomic_red_coor')
        if (index > 0):
            for j in range(natom):
                atom = lines[i+j+1].split()[0]
                atom = atom.replace("'", "")
                xred[0]  = lines[i+j+1].split()[1]
                xred[1]  = lines[i+j+1].split()[2]
                xred[2]  = lines[i+j+1].split()[3]
                xcart = xred*al
                atomkindid = lines[i+j+1].split()[4]
                atomdict.update({str(j): [atom, xred, xcart, atomkindid]})
    return al, natom, atomdict
def plot_crystal_structure(al, natom, atomdict):
    #fig = plt.figure(figsize=(4, 8))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for iatom in range(natom):
        value = atomdict[str(iatom)]
        xcart = value[2]
        if (value[0]=='Si'):
            ax.scatter(xcart[0],xcart[1],xcart[2],c='r',s=120,marker="o")
        elif (value[0]=='O'):
            ax.scatter(xcart[0],xcart[1],xcart[2],c='g',s=80,marker="o")
    ax.plot([0.0, al[0]],[0.0, 0.0],[0.0, 0.0],'k--')
    ax.plot([0.0, 0.0],[0.0, al[1]],[0.0, 0.0],'k--')
    ax.plot([0.0, 0.0],[0.0, 0.0],[0.0, al[2]],'k--')
    ax.plot([0.0, al[0]],[al[1],al[1]],[0.0, 0.0],'k--')
    ax.plot([0.0, 0.0],[0.0, al[1]],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[0.0, 0.0],[0.0, al[2]],'k--')
    ax.plot([0.0, al[0]],[0.0, 0.0],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[0.0, al[1]],[0.0, 0.0],'k--')
    ax.plot([al[0], al[0]],[0.0, 0.0],[0.0, al[2]],'k--')
    ax.plot([0.0, 0.0],[al[1],al[1]],[0.0, al[2]],'k--')
    ax.plot([0.0, al[0]],[al[1],al[1]],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[0.0, al[1]],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[al[1],al[1]],[0.0, al[2]],'k--')
    #ax.set_xlim([-1.0, 6.8])
    #ax.set_ylim([-1.0, 6.8])
    #ax.set_zlim([-1.0, 6.8])
    ax.axis('equal')
    plt.show()
def atomdict2bonddict(al, natom, atomdict, rc = 3.4):
    bonddict = {}
    k = 0
    vec = np.zeros(3,dtype='float64')
    for iatom in range(natom):
            value = atomdict[str(iatom)]
            atomi = value[0]
            xcart = value[2]
            for jatom in range(natom):
            #for jatom in range(iatom, natom):
                value = atomdict[str(jatom)]
                atomj = value[0]
                ycart = value[2]
                for ix in range(-1,2):
                    for iy in range(-1,2):
                        for iz in range(-1,2):
                            vec[0] = ix*al[0] 
                            vec[1] = iy*al[1]
                            vec[2] = iz*al[2]
                            ycart = value[2] + vec
                            if ((np.linalg.norm(xcart - ycart) < rc)and(atomi != atomj)):
                                bonddict.update({str(k): [iatom, xcart, atomi, jatom, ix, iy, iz, ycart, atomj]})
                                #print(atomi, atomj)
                                k += 1
                    #print(atomi, atomj,'hit!')
                    #print(iatom,jatom)
    nbond = 1*k
    print('rc', rc)
    print('nbond', nbond)
    print(bonddict.keys())
    return nbond, bonddict
def plot_bonds(al, natom, atomdict, nbond, bonddict):
    #fig = plt.figure(figsize=(8, 4))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    bondsum = np.zeros(3,dtype='float64')
    for iatom in range(natom):
        value = atomdict[str(iatom)]
        xcart = value[2]
        if (value[0]=='Si'):
            ax.scatter(xcart[0],xcart[1],xcart[2],c='r',s=120,marker="o")
        elif (value[0]=='O'):
            ax.scatter(xcart[0],xcart[1],xcart[2],c='g',s=80,marker="o")
    for ibond in range(nbond):
        value = bonddict[str(ibond)]
        iatom = value[0]
        xcart = value[1]
        atomi = value[2]
        jatom = value[3]
        ix = value[4]
        iy = value[5]
        iz = value[6]
        ycart = value[7]
        atomj = value[8]
        #print(xcart, ycart)
        if ((atomi == 'Si')and(atomj == 'O')): #Si (xcart) <= O (ycart)
            dcart = -ycart + xcart       # xcart + dcart = ycart
            ax.quiver(ycart[0],ycart[1],ycart[2],dcart[0],dcart[1],dcart[2])
        elif ((atomi == 'O')and(atomj == 'Si')): #O (xcart) => Si (ycart)
            dcart = ycart - xcart      # xcart - dcart = ycart
            ax.quiver(xcart[0],xcart[1],xcart[2],dcart[0],dcart[1],dcart[2])
        else:
            print('Error')
        bonddict.update({str(ibond): [iatom, xcart, atomi, jatom, ix, iy, iz, ycart, atomj, dcart]})
        bondsum += dcart
    ax.plot([0.0, al[0]],[0.0, 0.0],[0.0, 0.0],'k--')
    ax.plot([0.0, 0.0],[0.0, al[1]],[0.0, 0.0],'k--')
    ax.plot([0.0, 0.0],[0.0, 0.0],[0.0, al[2]],'k--')
    ax.plot([0.0, al[0]],[al[1],al[1]],[0.0, 0.0],'k--')
    ax.plot([0.0, 0.0],[0.0, al[1]],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[0.0, 0.0],[0.0, al[2]],'k--')
    ax.plot([0.0, al[0]],[0.0, 0.0],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[0.0, al[1]],[0.0, 0.0],'k--')
    ax.plot([al[0], al[0]],[0.0, 0.0],[0.0, al[2]],'k--')
    ax.plot([0.0, 0.0],[al[1],al[1]],[0.0, al[2]],'k--')
    ax.plot([0.0, al[0]],[al[1],al[1]],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[0.0, al[1]],[al[2],al[2]],'k--')
    ax.plot([al[0], al[0]],[al[1],al[1]],[0.0, al[2]],'k--')
    #ax.set_xlim([-1.0, 6.8])
    #ax.set_ylim([-1.0, 6.8])
    #ax.set_zlim([-1.0, 6.8])
    ax.axis('equal')
    plt.show()
    print('bondsum', bondsum)
    return bonddict