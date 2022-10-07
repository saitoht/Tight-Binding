import numpy as np
import scipy as sc
import scipy.linalg as sl
import matplotlib.pyplot as plt
import cmath as cm
import time, re, sys, itertools
import get_ham

wfi0 = 6
wfe0 = 23
var0 = -2
nelecu = 23 # number of majority electrons
nelecd = 19 # number of minority electrons
nwf = 23 # number of wannier orbitals: e.g. d+6p => 23
sw_red = False # e.g. consider onsite 2p contribution reduction or not
sw_spn = True # consider spin up and dn, or not 
sw_save = True # save data into textfile or not
filename = "Hopping"

def get_tij(fname,nwf):
    """ get hopping amplitude <wi(ri)|tij|wj(rj)> from file, usually named Hopping.up & Hopping.dn """
    hop_all = np.loadtxt(fname,dtype='float',comments='#',unpack=True,ndmin=0)
    hoplen = len(hop_all[0])
    rlen = nwf * nwf
    re_trij = hop_all[8] # real part of tij
    im_trij = hop_all[9] # img part of tij
    trij = re_trij + 1.j * im_trij
    trmat = np.zeros((rlen,nwf,nwf), dtype=np.complex128)
    for l in range(rlen):
        tr = trij[l*nwf*nwf:(l+1)*nwf*nwf]
        for m in range(nwf): 
            tr2 = tr[m*nwf:(m+1)*nwf]
            for n in range(nwf):
                trmat[l,m,n] = tr2[n] # tij
    xarr = hop_all[3,0:hoplen:nwf*nwf]
    yarr = hop_all[4,0:hoplen:nwf*nwf]
    zarr = hop_all[5,0:hoplen:nwf*nwf]
    rv = [xarr.tolist(),yarr.tolist(),zarr.tolist()] # r vector
    return(trmat, rv)

def get_hamkij(tr,rv,k):
    """ get k-dependent tight-binding Hamiltonian matrix """
    def ekr(k,r):
        kdotr = np.dot(np.array(k), np.array(r))
        ekr = np.cos(kdotr) + 1.j * np.sin(kdotr) # exp(ikr) plane wave
        return ekr

    rlen = nwf * nwf
    rv = np.array(rv)
    Hk = np.zeros((nwf,nwf), dtype=np.complex128)
    for n in range(nwf):
        for m in range(nwf):
            Hk_comp = 0.+0.j
            for l in range(rlen):
                phase = ekr(k,rv[:,l])
                Hk_comp = Hk_comp + tr[l,m,n] * phase # Hamiltonian for k and r
            Hk[m,n] = Hk_comp # Hamiltonian for k
    return Hk

def get_syml(foobar):
    """ get symmetry line from syml.foobar """
    with open('syml.'+str(foobar),'r') as f:
        txt = f.readlines()
        path0 = np.zeros((len(txt)-1, 2)).tolist()
        frac0 = np.zeros((len(txt)-1, 7)).tolist()
        for i in range(len(txt)-1):
            path = re.findall('[a-z]+',txt[i],flags=re.IGNORECASE)
            frac = re.findall(r'[-+]?\d*\.\d+|\d+', txt[i])
            path0[i] = path
            frac0[i] = frac
    print("Path: ", path0)
    print(frac0)
    path0 = np.array(path0)
    frac0 = np.array(frac0)
    klen0 = np.zeros(len(txt)-1).tolist()
    ks0 = np.zeros((len(txt)-1,3))
    ke0 = np.zeros((len(txt)-1,3))
    ksname0 = np.zeros(len(txt)-1).tolist()
    kename0 = np.zeros(len(txt)-1).tolist()
    for i in range(len(frac0)):
        klen0[i] = int(frac0[i,0])
        ksname0[i] = path0[i,0]
        kename0[i] = path0[i,1]
        for j in range(3):
            ks0[i,j] = float(frac0[i,j+1])
            ke0[i,j] = float(frac0[i,j+4])
    return(klen0, ks0, ke0, ksname0, kename0)

def mk_k(ks,ke,ksn,klen):
    """ make k point mesh for calculation """
    kdif = np.array(ke)-np.array(ks)
    kd_norm = 2.*np.pi*np.sqrt(np.dot(kdif, kdif))
    kvarr = [[],[],[]]
    for i in range(3):
        kvarr[i] = (np.linspace(2.*np.pi*ks[i], 2.*np.pi*ke[i], int(klen))).tolist()
    karr = np.linspace(ksn,ksn+kd_norm,int(klen)).tolist()
    kvarrT = np.array(kvarr).T
    return(kvarrT,karr,ksn+kd_norm)

def calc_chemical_potential(eigu, eigd, neu, ned, sw_semi):
    """ calculate chemical potential """
    ### now just for semiconductor
    ### only consider valence top level
    eigu = eigu.T
    eigd = eigd.T
    if sw_semi:
        eigut = eigu[neu-1].tolist()
        eigdt = eigd[ned-1].tolist()
        vtu = max(eigut)
        vtd = max(eigdt)
        vt = max([vtu,vtd])
    else:
        vt = 0.
    return vt

def reduce_onsite(hopu, hopd, wfi, wfe, var):
    """ reduce the onsite components of matrix elements """
    rlen = nwf*nwf
    print("-------------------------------------------")
    print("Reduced onsite: ")
    for m in range(wfi-1, wfe):
        hopu[0,m,m] = hopu[0,m,m] + var
        hopd[0,m,m] = hopd[0,m,m] + var
        print(m+1, hopu[0,m,m])
    return(hopu, hopd)

def main(foobar, fn):
    """ main program """
    (klen, ks, ke, ksname, kename) = get_syml(foobar)
    if sw_spn:
        (tr_up,r_up) = get_tij(str(fn)+'.up',nwf)
        (tr_dn,r_dn) = get_tij(str(fn)+'.dn',nwf)
        if sw_red:
            (tr_up,tr_dn) = reduce_onsite(tr_up,tr_dn,wfi0,wfe0,var0)
        else:
            pass
    else:
        (tr,r) = get_tij(str(fn)+'.dat',nwf)
    kvarrT = np.zeros(len(klen)).tolist()
    karr = np.zeros(len(klen)).tolist()
    ksn = np.zeros(len(klen)).tolist()
    ksn0 = np.dot(ks[0], ks[0])
    (kvarrT[0],karr[0],ksn[0])=mk_k(ks[0],ke[0],ksn0,klen[0])
    for i in range(1,len(klen)):
        (kvarrT[i],karr[i],ksn[i])=mk_k(ks[i],ke[i],ksn[i-1],klen[i])
    klen_all = sum(klen)
    if sw_spn:
        eig_up = np.zeros((klen_all,nwf))
        eig_dn = np.zeros((klen_all,nwf))
    else:
        eig = np.zeros((klen_all,nwf))
    count = 0 # count for calculation, e.g. 32 /50 (show progress)
    countj = 0
    kvarrT = np.array(kvarrT)
    print("------- Start diagonalozation! -------")
    for j in range(len(klen)):
        counti = 0
        for i in range(klen[j]):
            if countj != 0 and counti == 0:
                if sw_spn:
                    eig_up[count] = eig_up[count-1]
                    eig_dn[count] = eig_dn[count-1]
                else:
                    eig[count] = eig[count-1]
                count += 1
                counti += 1
                print(count, "/", klen_all)
                print("*skipped")
                continue
            else:
                pass
            st = time.time()
            if sw_spn:
                hamk_up = get_ham.get_hamkij(nwf, tr_up, r_up, kvarrT[j,i])
                hamk_dn = get_ham.get_hamkij(nwf, tr_dn, r_dn, kvarrT[j,i])
                (eigu,uniu) = sl.eigh(hamk_up)
                (eigd,unid) = sl.eigh(hamk_dn)
                eig_up[count] = eigu
                eig_dn[count] = eigd
            else:
                hamk = get_ham.get_hamkij(tr,r,kvarrT[j,i])
                (eig0,uni) = sl.eigh(hamk)
                eig[count] = eig0
            count += 1
            dt = time.time() - st
            print(count, "/", klen_all)
            if count == klen_all:
                pass
            else:
                print("*Estimated residual calculation time: ", dt * (klen_all - count), "sec")
            counti += 1
        countj += 1
    print("------- Calculation finished! -------")

    ### Saving data
    vtop = calc_chemical_potential(eig_up, eig_dn, nelecu, nelecd, True)
    karr_all = []
    for karri in karr:
        karr_all += karri
    xlist=karr_all
    if sw_spn:
        ylist_up = eig_up - vtop
        ylist_dn = eig_dn - vtop
    else:
        ylist = eig - vtop
        
    if sw_save:
        # Saving energy eigenvalues
        with open("eigu.dat", "w") as f:
            f.write("###  k    energy_eigenvalue \n")
            for j in range(nwf):
                for i in range(klen_all):
                    data = "{:.8f}   {:.8f} \n".format(float(karr_all[i]),float(ylist_up[i,j]))
                    f.write(data)
                        
        with open("eigd.dat", "w") as f:
            f.write("###  k    energy_eigenvalue \n")
            for j in range(nwf):
                for i in range(klen_all):
                    data = "{:.8f}   {:.8f} \n".format(float(karr_all[i]),float(ylist_dn[i,j]))
                    f.write(data)
    else:
        pass

    ### Plotting
    plt.rcParams['font.family'] = 'Helvetica'
    plt.rcParams["xtick.labelsize"]=15.0
    plt.rcParams["ytick.labelsize"]=15.0
    plt.rcParams["xtick.major.pad"] = 5
    plt.rcParams["ytick.major.pad"] = 5
    plt.rcParams["axes.labelsize"] = 18.0
    plt.rcParams["axes.linewidth"] = 1.0
    plt.rcParams["axes.labelpad"] = 6
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.major.width"] = 1.0
    plt.rcParams["ytick.major.width"] = 1.0
    plt.rcParams["xtick.minor.width"] = 0.5
    plt.rcParams["ytick.minor.width"] = 0.5
    plt.rcParams["xtick.major.size"] = 4.5
    plt.rcParams["ytick.major.size"] = 4.5
    plt.rcParams["xtick.minor.size"] = 3.0
    plt.rcParams["ytick.minor.size"] = 3.0

    ksn.insert(0,ksn0)
    ksname.append(kename[len(kename)-1])
    kpath = [ '$\Gamma$' if kp=='G' else kp for kp in ksname]
    xlab = ' '
    ylab = 'Energy (eV)'
    tit = 'Tight-Binding Band'
    ylimmin = -10
    ylimmax = 15
    if sw_spn:
        # Majority
        figu = plt.figure()
        ax = figu.add_subplot(111,xlabel=xlab, ylabel=ylab, title=tit+' Majority',xlim=(ksn[0],ksn[len(ksn)-1]),ylim=(ylimmin,ylimmax))
        ax.plot(xlist,ylist_up,c='black',lw=1.)
        ax.plot([:,:], [0., 0.], color="gray", lw=0.5, ls='dashed')
        for i in range(len(ksn)):
            ax.plot([ksn[i], ksn[i]], [ylimmin, ylimmax], color="gray", lw=0.5, ls='dashed')
        plt.xticks(ksn, kpath)
        figu.savefig("TBband_up.pdf")
        plt.show()

        # Minority
        figd = plt.figure()
        ax = figd.add_subplot(111,xlabel=xlab, ylabel=ylab, title=tit+' Minority',xlim=(ksn[0],ksn[len(ksn)-1]),ylim=(ylimmin,ylimmax))
        ax.plot(xlist,ylist_dn,c='black',lw=1.)
        ax.plot([:,:], [0., 0.], color="gray", lw=0.5, ls='dashed')
        for i in range(len(ksn)):
            ax.plot([ksn[i], ksn[i]], [ylimmin, ylimmax], color="black", lw=0.5, ls='dashed')
        plt.xticks(ksn, kpath)
        figd.savefig("TBband_dn.pdf")
        plt.show()
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111,xlabel=xlab, ylabel=ylab, title=tit ,xlim=(ksn[0],ksn[len(ksn)-1]),ylim=(ylimmin,ylimmax))
        ax.plot(xlist,ylist,c='black',lw=1.)
        ax.plot([:,:], [0., 0.], color="gray", lw=0.5, ls='dashed')
        for i in range(len(ksn)):
            ax.plot([ksn[i], ksn[i]], [ylimmin, ylimmax], color="black", lw=0.5, ls='dashed')
        plt.xticks(ksn, kpath)
        figu.savefig("TBband.pdf")
        plt.show()

if __name__=="__main__":
    args = sys.argv
    if len(args) == 1 or args[1] == "--help" or args[1] == "-h":
        help_message = """ Usage of TBband.py:
        python TBband.py {foobar} """
        print(help_message)
        exit()
    else:
        for i in args[1:]:
           starttime = time.time()
           foobar0 = args[1]
           main(foobar0, filename)
           endtime = time.time()
           print("Calculation time: ", endtime-starttime)
