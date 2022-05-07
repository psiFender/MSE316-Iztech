#!/usr/bin/env python
# -*- coding: utf-8 -*- 
"""
VASP'in EIGENVAL ciktisindan k-point ve band enerjilerini okur ve bandlari cizer.
ISPIN=2 ise up- ve down- spinler ayri ayri cizilir.
Fermi enerjisi DOSCAR'dan okunur
"""
from pylab import *
import numpy

def main():
    # Efermi oku
    Efermi = f_Efermi_oku('DOSCAR')
    # reciprocal lattice vector'ler vasprun.xml'den
    rec_basis = f_read_rec_basis()
    print 'rec_basis = ',rec_basis
    
    # EIGENVAL oku
    [ispin,nk,nband,kpoint,band]=f_read_EIGENVAL('EIGENVAL',Efermi,rec_basis)
    if ispin==2:
        [band_up,band_down]=band
    else:
        band_up=band
    # high symmetry point'ler vasprun.xml'den okunsun
    high_symmetry_points = f_read_HighSymmetryPoints(rec_basis)
    print 'high_symmetry_points = ',high_symmetry_points
    # çizim için gerekli olan ölçeklendirilmiş kpoint array yarat (k4plot)
    nk_per_section        = nk/(len(high_symmetry_points)-1)
    print 'nk_per_section =',nk_per_section
    [k4plot,section4plot] = f_k4plot(nk,kpoint,nk_per_section,high_symmetry_points,rec_basis)

    # ciz
    for b in range (0,nband):
        plot(k4plot,band_up[:,b],'-b')  
        if ispin==2:
            plot(k4plot,band_down[:,b],'-r')
            hold(True)
    hold(True)
    ylabel('Energy (eV)')
    if ispin==2:
        legend(['up','down'])
    # sections
    #print section4plot
    min_ener = band_up.min()
    max_ener = band_up.max()
    for s in range (1,len(high_symmetry_points)):
        xx = [section4plot[s],section4plot[s]]
        plot(xx,[min_ener,max_ener],'-k')
        hold(True)
    plot([0,k4plot[nk-1]],[0,0],'--k')
    ylim([-13,+12])
    xlim([0,k4plot[nk-1]])
    xticks(section4plot,['$\Gamma$','M','K','$\Gamma$'])
    #utoscale(enable=True,axis='x', tight=True)
    show()
    
    ### bandları dosyaya kaydet
    # kpoints
    dosya  = open('band_kpoints.dat','w');
    for kk in range (0,nk):
        satir = " %+1.10E \n" % k4plot[kk]
        dosya.write(satir)
    dosya.close()
    # energies
    dosya = open('band_energies.dat','w');
    for kk in range (0,nk):
        for bb in range (0,nband):
            aux = "%+1.6E " % band_up[kk,bb]
            dosya.write(aux)
        aux = "\n"
        dosya.write(aux)
    dosya.close()
    
    return 0

def f_read_HighSymmetryPoints(rec_basis):
    HighSymmetryPoints = []
    dosya  = open('vasprun.xml','r')
    buldum = 0
    while buldum==0:
        satir = dosya.readline()
        if '<kpoints>' in satir:
            buldum = 1
    dosya.readline();dosya.readline(); # gereksiz satirlar
    # read high symmetry points
    buldum = 0
    while buldum==0:
        aux = dosya.readline().split()
        if aux[0]=='<v>':
            aux1 = [aux[1],aux[2],aux[3][1:6]]
            aux2 = [float(s) for s in aux1]
            print 'hsp_reci=',aux2
            aux3 = frac2cart(aux2,rec_basis)
            HighSymmetryPoints.append(aux3)
        else:
            buldum=1
    hsp = numpy.asarray(HighSymmetryPoints)
    return hsp


def f_read_rec_basis():
    rec_basis =[]
    dosya  = open('vasprun.xml','r')
    buldum = 0
    while buldum==0:
        satir = dosya.readline()
        if 'finalpos' in satir:
            buldum = 1
    for i in range (0,8):
        dosya.readline(); # gereksiz satirlar    
    for i in range (0,3):
        aux  = dosya.readline().split()
        aux1 = [aux[1],aux[2],aux[3]]
        aux2 = [float(s) for s in aux1]
        rec_basis.append(aux2)
        basis = numpy.array(rec_basis)
    return basis
    
def f_read_EIGENVAL(dosyaadi,Efermi,rec_basis):
	# band verilerini dosyadan okuyacagiz
	dosya = open(dosyaadi,'r')
	# ispin oku
	ispin = [int(float(s)) for s in dosya.readline().split()][3]
	print 'ispin = ',ispin
	# sonraki 4 satir isimize yaramiyor
	[dosya.readline() for i in range (0,4)]
	# altinci satirdan nk ve nband oku
	[aux,nk,nband]=[int(float(i)) for i in dosya.readline().split()]
	print 'nk    =',nk
	print 'nband = ',nband
	# band ve kpoint verileri icin bos matris hazirla
	band_up   = zeros((nk,nband))
	if ispin==2:
		band_down = zeros((nk,nband))
	kpoint    = zeros((nk,3))
	# oku
	for k in range (0,nk):
		dosya.readline()
		aux         = array([float(s) for s in dosya.readline().split()][0:3])
		kpoint[k,:] = frac2cart(aux,rec_basis)
		for b in range (0,nband):
			aux            = [float(s) for s in dosya.readline().split()]
			band_up[k,b]   = aux[1]-Efermi
			if ispin==2:
				band_down[k,b] = aux[2]-Efermi
	# çıktı hazırla
	cikti = [ispin,nk,nband,kpoint,band_up]
	if ispin==2:
		bands = [band_up,band_down]
		cikti = [ispin,nk,nband,kpoint,bands]
	dosya.close()
	return(cikti)

def f_k4plot(nk,kpoints,nk_per_section,high_symmetry_points,rec_basis):
    say          = 0
    k4plot       = []
    section4plot = [0.]
    n_section = len(high_symmetry_points)
    for sec in range (0,n_section-1):
        aux1   = high_symmetry_points[sec+1]-high_symmetry_points[sec]
        aux2   = frac2cart(aux1,rec_basis)
        dk     = norm(aux2)/nk_per_section
        print high_symmetry_points[sec],norm(high_symmetry_points[sec]),dk
        for kk in range (0,nk_per_section):
            if say>0:
                aux = k4plot[say-1]+dk
            else:
                aux = 0.
            k4plot.append(aux)
            say = say + 1
            if kk==nk_per_section-1:
                section4plot.append(aux)
    k4plot       = k4plot#/aux
    section4plot = section4plot#/aux
    #print k4plot,len(k4plot)
    print 'section4plot=',section4plot
    return(k4plot,section4plot)
    
    
def frac2cart(frac,basis):
    aux1    = numpy.array(frac)
    cart    = numpy.dot(aux1,basis)
    return cart

def f_Efermi_oku(dosyaadi):
	dosya  = open(dosyaadi,'r')
	for s in range (0,5):
		dosya.readline()
	Efermi = float(dosya.readline().split()[3])
	return(Efermi)

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False	
	

if __name__ == '__main__':
	main()

