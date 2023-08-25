# Plot of the main quantities. 

import numpy as np
import math 
import sys 
import os 
import matplotlib.pyplot as plt


def vel(dbrotAvg,nu,rho):

    plt.figure(figsize=(13,10))
    plt.xscale('log')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Uz']/dbrotAvg['u_tau'])  # Radial coordinate scaled for R=1
    plt.xlabel('$y+$',fontsize=16)
    plt.ylabel('$<u_z>+$',fontsize=16)
    plt.show()
    return


def kBudgets(dbrotAvg,nu,rho):


    plt.figure(figsize=(13,10))
    plt.xscale('log')
    fac = nu/np.power(dbrotAvg['u_tau'],4)
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['P_k']*fac,  label='Production')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['T_k']*fac,  label='Turbulent Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PS_k']*fac, label='Pressure Strain')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PT_k']*fac, label='Pressure Transport')
#    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Pi_k']*fac)
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['VD_k']*fac, label='Viscous Diffusion')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['D_k']*fac,  label='Dissipation')
#    plt.plot(dbrotAvg['r'], dbrotAvg['C_k']*fac)
    plt.xlabel('$y+$',fontsize=16)
    plt.legend(fontsize=16)
    plt.title('TKE Budget',fontsize=16)
    plt.show()
    return


def Re_stresstensor(dbrotAvg,nu,rho):


    plt.figure(figsize=(13,10))
    plt.xscale('log')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], (np.sqrt(dbrotAvg['urur']-dbrotAvg['Ur']*dbrotAvg['Ur']))/dbrotAvg['u_tau'],       label='$<u_{r,rms}^{\'+}>$')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], (dbrotAvg['urut']-dbrotAvg['Ur']*dbrotAvg['Ut'])/np.power(dbrotAvg['u_tau'],2),    label='$<u_r^{\'} u_{\\theta}^{\'+}>$')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], (dbrotAvg['uruz']-dbrotAvg['Ur']*dbrotAvg['Uz'])/np.power(dbrotAvg['u_tau'],2),    label='$<u_r^{\'} u_z^{\'+}>$')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], (np.sqrt(dbrotAvg['utut']-dbrotAvg['Ut']*dbrotAvg['Ut']))/dbrotAvg['u_tau'],       label='$<u_{\\theta,rms}^{\'+}>$')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], (np.sqrt(abs(dbrotAvg['uzuz']-dbrotAvg['Uz']*dbrotAvg['Uz'])))/dbrotAvg['u_tau'],  label='$<u_{z,rms}^{\'+}>$')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], (dbrotAvg['utuz']-dbrotAvg['Ut']*dbrotAvg['Uz'])/np.power(dbrotAvg['u_tau'],2),    label='$<u_t^{\'} u_z^{\'+}>$')
    plt.xlabel('$y+$',fontsize=16)
    plt.legend(fontsize=16)
    plt.title('Velocity 2nd-order moments',fontsize=16)
    plt.show()
    return


def urur(dbrotAvg,nu,rho):

    fac = nu/np.power(dbrotAvg['u_tau'],4)
    plt.figure(figsize=(13,10))
    plt.xscale('log')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Prr']*fac,  label='Production')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Trr']*fac,  label='Turbulent Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PSrr']*fac, label='Pressure Strain')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PTrr']*fac, label='Pressure Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['VDrr']*fac, label='Viscous Diffusion')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Drr']*fac,  label='Dissipation')
#    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Crr']*fac)
    plt.xlabel('$y+$',fontsize=16)
    plt.legend(fontsize=16)
    plt.title('Budget of $<u_r^{\'} u_r^{\'}>$',fontsize=16)
    plt.show()
    return


def utut(dbrotAvg,nu,rho):

    fac = nu/np.power(dbrotAvg['u_tau'],4)
    plt.figure(figsize=(13,10))
    plt.xscale('log')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Ptt']*fac,  label='Production')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Ttt']*fac,  label='Turbulent Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PStt']*fac, label='Pressure Strain')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PTtt']*fac, label='Pressure Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['VDtt']*fac, label='Viscous Diffusion')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Dtt']*fac,  label='Dissipation')
#    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Ctt']*fac)
    plt.xlabel('$y+$',fontsize=16)
    plt.legend(fontsize=16)
    plt.title('Budget of $<u_{\\theta}^{\'} u_{\\theta}^{\'}>$',fontsize=16)
    plt.show()
    return


def uzuz(dbrotAvg,nu,rho):

    fac = nu/np.power(dbrotAvg['u_tau'],4)
    plt.figure(figsize=(13,10))
    plt.xscale('log')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Pzz']*fac,  label='Production')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Tzz']*fac,  label='Turbulent Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PSzz']*fac, label='Pressure Strain')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PTzz']*fac, label='Pressure Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['VDzz']*fac, label='Viscous Diffusion')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Dzz']*fac,  label='Dissipation')
#    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Czz']*fac)
    plt.xlabel('$y+$',fontsize=16)
    plt.legend(fontsize=16)
    plt.title('Budget of $<u_z^{\'} u_z^{\'}>$',fontsize=16)
    plt.show()
    return


def uruz(dbrotAvg,nu,rho):

    fac = nu/np.power(dbrotAvg['u_tau'],4)
    plt.figure(figsize=(13,10))
    plt.xscale('log')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Prz']*fac,  label='Production')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Trz']*fac,  label='Turbulent Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PSrz']*fac, label='Pressure Strain')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['PTrz']*fac, label='Pressure Transport')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['VDrz']*fac, label='Viscous Diffusion')
    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Drz']*fac,  label='Dissipation')
#    plt.plot(dbrotAvg['r']*dbrotAvg['Re_tau'], dbrotAvg['Crz']*fac)
    plt.xlabel('$y+$',fontsize=16)
    plt.legend(fontsize=16)
    plt.title('Budget of $<u_r^{\'} u_z^{\'}>$',fontsize=16)
    plt.show()
    return

