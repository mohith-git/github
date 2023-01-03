# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:10:15 2022

@author: MOHITHTHAMANA
"""

from flask import Flask, redirect, url_for, render_template, request, Markup
import numpy as np
from plotly.offline import plot
from plotly.graph_objs import Scatter
import plotly.express as px
import plotly.graph_objects as go
import os
import pandas as pd
import argparse
import json
import requests
import time

import plotly.figure_factory as ff

STATIC_DIR = os.path.abspath(r'C:\Users\MOHITHTHAMANA\Desktop\Transmembrane\static')
app = Flask(__name__,static_folder=STATIC_DIR)
app = Flask(__name__)


#Basic functions:
def read_file(f):
    name = os.path.join(app.root_path,'pdbf',f+'.pdb')
    
    file = open(name)
    a=file.read()
    a=a.splitlines()
    return a

def read_xml(f):
        import xml.etree.ElementTree as ET
        seq_xml = []
        seq = []
        name = os.path.join(app.root_path,'xmlf',f+'.xml')
        tree = ET.parse(name)
        root = tree.getroot()
        for i in range(len(root)):
            if 'CHAIN' in root[i].tag and root[i].attrib['CHAINID'] in ['A','P','M']:
                for j in range(len(root[i])):
                    if j!=0:
                        seq_xml.append(root[i][j].attrib)
        for i in root:
            if 'CHAIN' in i.tag and i.attrib['CHAINID'] in ['A','P','M']:
               for j in i:
                 if 'SEQ' in j.tag:
                     seq = ''.join(str(j.text).split(' '))
        return([seq_xml,seq])



def distance_formula(i,j):
        i[6]=float(i[6])
        i[7]=float(i[7])
        i[8]=float(i[8])
        j[6]=float(j[6])
        j[7]=float(j[7])
        j[8]=float(j[8])
        d=(i[6]-j[6])**2+(i[7]-j[7])**2+(i[8]-j[8])**2
        return d

def arm_arm(f):

    arm,aromatic,aromatic_trp=[],[],[]


    fi = read_file(f) 
    for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2] in ['CG','CD1','CD2','CE1','CE2','CZ'] and aa[3] in ['PHE','TYR']:
                aromatic.append(aa)
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2] in ['CD2','CZ2','CZ3','CE2','CE3','CH2'] and aa[3] in ['TRP']:
                aromatic_trp.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
    
    for i in range(0,len(aromatic)-6,6):
            px=(float(aromatic[i][6])+float(aromatic[i+1][6])+float(aromatic[i+2][6])+float(aromatic[i+3][6])+float(aromatic[i+4][6])+float(aromatic[i+5][6]))/6
            py=(float(aromatic[i][7])+float(aromatic[i+1][7])+float(aromatic[i+2][7])+float(aromatic[i+3][7])+float(aromatic[i+4][7])+float(aromatic[i+5][7]))/6
            pz=(float(aromatic[i][8])+float(aromatic[i+1][8])+float(aromatic[i+2][8])+float(aromatic[i+3][8])+float(aromatic[i+4][8])+float(aromatic[i+5][8]))/6
            arm.append([aromatic[i][3],aromatic[i][5],px,py,pz])
    for i in range(0,len(aromatic_trp)-6,6):
            px=(float(aromatic_trp[i][6])+float(aromatic_trp[i+1][6])+float(aromatic_trp[i+2][6])+float(aromatic_trp[i+3][6])+float(aromatic_trp[i+4][6])+float(aromatic_trp[i+5][6]))/6
            py=(float(aromatic_trp[i][7])+float(aromatic_trp[i+1][7])+float(aromatic_trp[i+2][7])+float(aromatic_trp[i+3][7])+float(aromatic_trp[i+4][7])+float(aromatic_trp[i+5][7]))/6
            pz=(float(aromatic_trp[i][8])+float(aromatic_trp[i+1][8])+float(aromatic_trp[i+2][8])+float(aromatic_trp[i+3][8])+float(aromatic_trp[i+4][8])+float(aromatic_trp[i+5][8]))/6
            arm.append([aromatic_trp[i][3],aromatic_trp[i][5],px,py,pz])

    return [arm,aromatic,aromatic_trp]

def Short_range(f,cu):
        i1=0
        # read the pdb file
        fi = read_file(f)
        xa=[]
        sr=[]
        #
        #Storing the x,y,z coordinates of C-alpha atoms of each residue 
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        #
        #
        for i in xa:
            j1=0
            for j in xa:
                #calculating distance between any two residue
                d=distance_formula(i,j)
                #checking with corresponding conditions for short range contacts
                if 0<d**0.5<float(cu) and [j[3],j[5],i[3],i[5],round(d**0.5,2)] not in sr and abs(i1-j1)<=2 :
                    sr.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return sr
    
def Medium_range(f,cu):
        i1=0
        # read the pdb file
        fi = read_file(f)
        xa=[]
        mr=[]
        #
        #Storing the x,y,z coordinates of C-alpha atoms of each residue 
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        #
        # 
        for i in xa:
            j1=0
            for j in xa:
                #calculating distance between any two residue
                d=distance_formula(i,j)
                #checking with corresponding conditions for medium range contacts 
                if 0<d**0.5<float(cu) and [j[3],j[5],i[3],i[5],round(d**0.5,2)] not in mr and 3<=abs(i1-j1)<=4 :
                    mr.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return mr
    

def Long_range(f,cu):
        i1=0
        # read the pdb file
        fi = read_file(f)
        xa=[]
        lr=[]
        #
        #Storing the x,y,z coordinates of C-alpha atoms of each residue 
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        #
        # 
        for i in xa:
            j1=0
            for j in xa:
                #calculating distance between any two residue
                d=distance_formula(i,j)
                #checking with corresponding conditions for long range contacts
                if 0<d**0.5<float(cu) and [j[3],j[5],i[3],i[5],round(d**0.5,2)] not in lr and abs(i1-j1)>5 :
                    lr.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return lr

def LRO(f):
        i1=0
        # read the pdb file
        fi = read_file(f)
        xa=[]
        lro = []
        res=[]
        #
        #Storing the x,y,z coordinates of C-alpha atoms of each residue 
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        

        for i in range(len(xa)):
            
            s=0
            for j in range(i,len(xa)):
                #calculating distance between any two residue
                d=distance_formula(xa[i],xa[j])
                # Applying the distance cutoff to 8 angstorms and residue number separation more than 12 in number
                if d**0.5<8 and abs(i-j)>12 :
                    #summation over total length of the protein
                    s=s+(1/len(xa))
                
            lro.append(s)
            res.append([xa[i][3],xa[i][5],'A',round(s,3)])
            
        #output LRO for each residue , Total LRO 
        return res,round(sum(lro),3)
    
    
    
def Contact_Order(f):
        i1=0
        # read the pdb file
        fi = read_file(f)
        xa=[]
        
        l = 0
        summ=0
        #
        #Storing the x,y,z coordinates of C-alpha atoms of each residue
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        
        #
        #calculating the Contact Order
        c_order = 0
        for i in range(len(xa)):
          
          
          for j in range(i,len(xa)):
            #calculating distance between any two residue
            d=distance_formula(xa[i],xa[j])
            # Applying the distance cutoff to 6 angstorms and summation over the residue number separation
            if d**0.5<6 and i!=j:
                l=l+1
                c_order=c_order+(abs(i-j))
          
        # Contact Order = Summation / ((Total no of contacts) * (Length of the Protein))
        return(round((c_order)/(len(xa)*l),2))




def A8(f,inp,a):
        fi = read_file(f)
        
        xa=[]
        
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]==a:
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        
        #degree = {}
        contacts = []
        for i in xa:
            j1=0
            deg = 0
            res = ''
            for j in xa:
                d=distance_formula(i,j)
                if 0<d**0.5<=float(inp):
                    deg+=1
                    res+='{0}-{1}'.format(j[3],j[5])+' '

                j1=j1+1

            #degree['{0}-{1}'.format(i[3],i[5])] = deg
            contacts.append([i[3],i[5],'A',deg,res])
        #print('C-alpha:{0}'.format(inp))
        return(contacts)
    # def A14(self):
    #     inp = float(input())
    #     degree = {}
    #     for i in self.xa:
    #         j1=0
    #         deg = 0
    #         for j in self.xa:
    #             d=c.distance_formula(i,j)
    #             if 0<d**0.5<=inp:
    #                 deg+=1
    #             j1=j1+1
    #         degree['{0}-{1}'.format(i[3],i[5])] = deg
    #     print(degree)
def B8(f,inp):
        fi = read_file(f)
        
        xb=[]
        
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CB':
                xb.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        
        #degree = {}
        contacts = []
        for i in xb:
            j1=0
            deg = 0
            res = ''
            for j in xb:
                d=distance_formula(i,j)
                if d**0.5<=float(inp):
                    deg+=1
                    res+='{0}-{1}'.format(j[3],j[5])+' '

                j1=j1+1
            #degree['{0}-{1}'.format(i[3],i[5])] = deg
            contacts.append([i[3],i[5],'A',deg,res])
        #print('C-alpha:{0}'.format(inp))
        return(contacts)




def Disulphide_interaction(f):
        
        fi = read_file(f)
        cys1,cys = [],[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='SG':
                cys1.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        
        for i in cys1:
            if i[3]=='CYS':
                cys.append(i)
        
        p=[]
        i1=0
        for i in cys:
            j1=0
            for j in cys:
                d=distance_formula(i,j)
                if 0<d**0.5<=2.2 and [j[3],j[5],i[3],i[5],round(d**0.5,2)] not in p:
                    p.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return(p)
def Hydrophobic_interaction(f):
        fi = read_file(f)
        p,xb=[],[]
        i1=0
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CB':
                xb.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        for i in xb:
          if i[3] in ['ALA','VAL','LEU','TYR','PRO','TRP','MET','PHE','ILE']:
            j1=0
            for j in xb:
              if j[3] in ['ALA','VAL','LEU','TYR','PRO','TRP','MET','PHE','ILE']:
                d=distance_formula(i,j)
                if 0<d**0.5<=5 and [j[3],j[5],i[3],i[5],d**0.5] not in p:
                    p.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return(p)



def Cation_pi(f):
        fi = read_file(f)
        arm = arm_arm(f)
        cation = []
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2] in ['NH1','NH2','NZ','ND1','NE2'] and aa[3] in ['ARG','LYS','HIS']:
                cation.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break  
        cpi=[]
        i1=0
        for i in cation:
         if i[3]!='HIS':
            j1=0
            for j in arm[0]:
                d=(float(i[6])-j[2])**2+(float(i[7])-j[3])**2+(float(i[8])-j[4])**2
                if 0<d**0.5<10 and [j[0],j[1],i[3],i[5],i[2],d**0.5] not in cpi:
                    cpi.append([i[3],i[5],i[2],j[0],j[1],'-',round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return((cpi))




def Aromatic_Aromatic(f):
        p=[]
        
        arm = arm_arm(f)
        
        for i in arm[0]:
            for j in arm[0]:
                d=((i[2]-j[2])**2+(i[3]-j[3])**2+(i[4]-j[4])**2)
                if 4.5<d**0.5<7 and [j[0],j[1],i[0],i[1],round(d**0.5,2)] not in p:
                    p.append([i[0],i[1],j[0],j[1],round(d**0.5,2)])
        return p



    
def Main_Main_H_Bond(f):
        fi = read_file(f)
        m_n,m_s,m_o=[],[],[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['N','NE','NH1','NH2','NZ','ND1','NE2','NE1','ND2']:
                m_n.append(aa)
            elif aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['O','OE2','OE1','OD1','OD2','OG1','OG2','OH','OG']:
                m_o.append(aa)
            elif aa[0]=='ATOM' and aa[4]=='A' and aa[2]==['SG','SD']:
                m_s.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        h_bond = []
        i1=0
        for i in m_n:
         if i[2] == 'N':
            j1=0
            
            for j in m_o:
              if j[2] == 'O':
               if abs(int(i[5]) - int(j[5])) > 1:
                d=distance_formula(i,j)
                if 0<d**0.5<3.5 and [j[3],j[5],j[2],i[3],i[5],i[2],round(d**0.5,2)] not in h_bond:
                    h_bond.append([i[3],i[5],i[2],j[3],j[5],j[2],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        h_bond = sorted(h_bond,key=lambda x : x[1])
        return((h_bond))
    
    
def Main_Side_Chain_H_Bond(f):
        fi = read_file(f)
        m_n,m_s,m_o=[],[],[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['N','NE','NH1','NH2','NZ','ND1','NE2','NE1','ND2']:
                m_n.append(aa)
            elif aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['O','OE2','OE1','OD1','OD2','OG1','OG2','OH','OG']:
                m_o.append(aa)
            elif aa[0]=='ATOM' and aa[4]=='A' and aa[2]==['SG','SD']:
                m_s.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        
        atoms= [*m_n,*m_o,*m_s]
        h_bond = []
        real = []
        i1=0
        for i in atoms:
           
            j1=0
            
            
            for j in atoms:
             
               if abs(int(i[5]) - int(j[5])) > 1:
                
                 d=distance_formula(i,j)
                
                 if 0<d**0.5<3.5:
                    if i[2] == 'O':
                      h_bond.append([j[3],j[5],j[2],i[3],i[5],i[2],round(d**0.5,2)])
                    else:
                        if j[2] == 'N':
                          h_bond.append([j[3],j[5],j[2],i[3],i[5],i[2],round(d**0.5,2)])   
                     
               j1=j1+1
            i1=i1+1
        for i in h_bond:
           if i[2] in ['N','O'] and i[5]in ['N','O']:
               pass
           elif i[2] in ['N','O'] or i[5]in ['N','O']:
               real.append(i)
        real = sorted(real,key=lambda x : x[1])
        return(((real)))   
    
    
    
def Side_Side_Chain_H_Bond(f):
        fi = read_file(f)
        m_n,m_s,m_o=[],[],[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['N','NE','NH1','NH2','NZ','ND1','NE2','NE1','ND2']:
                m_n.append(aa)
            elif aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['O','OE2','OE1','OD1','OD2','OG1','OG2','OH','OG']:
                m_o.append(aa)
            elif aa[0]=='ATOM' and aa[4]=='A' and aa[2]==['SG','SD']:
                m_s.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        atoms= [*m_n,*m_o,*m_s]
        h_bond = []
        
        i1=0
        for i in atoms:
           if i[2] not in ['N','O']:
            j1=0
            
            
            for j in atoms:
              if j[2] not in ['N','O']:
               if abs(int(i[5]) - int(j[5])) > 1:
                d=distance_formula(i,j)
                if 0<d**0.5<3.5 and [j[3],j[5],j[2],i[3],i[5],i[2],round(d**0.5,2)] not in h_bond:
                    h_bond.append([i[3],i[5],i[2],j[3],j[5],j[2],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        
        return(((h_bond)))   
    
def Ionic_Interaction(f):
        fi = read_file(f)
        cation,anion=[],[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2] in ['NH1','NH2','NZ','ND1','NE2'] and aa[3] in ['ARG','LYS','HIS']:
                cation.append(aa)
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2] in ['OE2','OE1','OD1','OD2'] and aa[3] in ['ASP','GLU']:
                anion.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break  
        
        ion = []
        i1=0
        for i in cation:
            j1=0
           
            for j in anion:
                d=distance_formula(i,j)
                if 0<d**0.5<6:
                    ion.append([i[3],i[5],i[2],j[3],j[5],j[2],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
         
        
        
        return(ion)
    
def Surr_Hydrophob(f):
         # read the pdb file
         fi = read_file(f)
         xa=[]
         
         #
         #Storing the x,y,z coordinates of C-alpha atoms of each residue
         for i in fi:
             aa=i.split()
             if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                 xa.append(aa)
             if aa[0] == 'MODEL' and aa[1] == '2':
                 break
         #
         #
         # Hydrophobicity index (kcal/mol) given by Tanford and Jones (1971) 
         hydro_index = {'ALA': 0.87 , 'ASP': 0.66, 'CYS': 1.52, 'GLU': 0.67, 'PHE': 2.87 , 'GLY':0.1, 'HIS': 0.87, 'ILE': 3.15,
                   'LYS': 1.64, 'LEU': 2.17,'MET':1.67 ,'ASN': 0.09,'PRO': 2.77,'GLN': 0 ,'ARG':0.85, 'SER': 0.07,
                   'THR':0.07, 'VAL':1.87, 'TRP': 3.77, 'TYR': 2.67}
         #
         #
         #
         surr = []
        
         i1=0
         for i in xa:
             j1=0
             interact = []
             score = 0
             for j in xa:
                 #calculating distance between any two residue
                 d=distance_formula(i,j)
                 if 0<d**0.5<8:
                     #Storing the interacting residues using a standard cutoff distance of 8 Angstorms
                     interact.append(j[3])
             #
             #
             # Surrounding Hydrophobicity = summation (total number of surrounding residues of type j around ith residue of the protein * hydrophobicity index of ith residue)  
             for k in set(interact):
                    if len(k)==4:
                     score+=(interact.count(k)*hydro_index[k[1:]])
                    else:
                     score+=(interact.count(k)*hydro_index[k]) 

             j1=j1+1
             surr.append([i[3],i[5],round(score,2)])
                
             i1=i1+1
         
         return((surr))
# @app.route("/",methods=['POST','GET'])
# def first():
#     if request.method=='POST':
#       f = request.form.get('Pdb')
      
#       return render_template('layout.html')
#     return render_template('first.html')



    
def Regions(f):
        extr,intr,helix= [],[],[]
        xm = read_xml(f)

        for i in xm[0]:
            if i['type'] == '1':
                for k in range(int(i['seq_beg']),int(i['seq_end'])+1):
                    extr.append(str(k))
            elif i['type'] in ['B','H','I'] : 
                for k in range(int(i['seq_beg']),int(i['seq_end'])+1):
                    helix.append(str(k))
            else:
                for k in range(int(i['seq_beg']),int(i['seq_end'])+1): 
                    intr.append(str(k))
        #print(helix)
        
        return [extr,helix,intr,xm[1]]




def c_m(f):
        cm=[]
        i1=0
        fi = read_file(f)
        xa=[]
        
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        for i in xa:
            j1=0
            for j in xa:
                d=distance_formula(i,j)
                if 0<=d**0.5<8 and [j[3],j[5],i[3],i[5],d**0.5] not in cm:#basically we are appending so changing its direction do no append again
                    cm.append([i[3],i[5],j[3],j[5],d**0.5])
                j1=j1+1
            i1=i1+1
        x=[]
        y=[]
        z=[]
        for i in cm:
            x.append(int(i[1]))
            y.append(int(i[3]))
            z.append(int(i[4]))
        for i in cm:
            if i[3]!=i[1]:
                x.append(int(i[3]))
                y.append(int(i[1]))
                z.append(int(i[4]))

        return [x,y,z]

   
def c_n(f):
        x=[]
        y=[]
        z=[]
        fi = read_file(f)
        xa=[]
        
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4] in ['A','P','M'] and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        for i in xa:
            x.append(float(i[6]))
            y.append(float(i[7]))
            z.append(float(i[8]))

        return [x,y,z]

def pdb_to_hssp(pdb_file_path, rest_url):
    # Read the pdb file data into a variable
    files = {'file_': open(pdb_file_path, 'rb')}

    # Send a request to the server to create hssp data from the pdb file data.
    # If an error occurs, an exception is raised and the program exits. If the
    # request is successful, the id of the job running on the server is
    # returned.
    url_create = '{}api/create/pdb_file/dssp/'.format(rest_url)
    r = requests.post(url_create, files=files)
    r.raise_for_status()

    job_id = json.loads(r.text)['id']
    print("Job submitted successfully. Id is: '{}'".format(job_id))

    # Loop until the job running on the server has finished, either successfully
    # or due to an error.
    ready = False
    while not ready:
        # Check the status of the running job. If an error occurs an exception
        # is raised and the program exits. If the request is successful, the
        # status is returned.
        url_status = '{}api/status/pdb_file/dssp/{}/'.format(rest_url,
                                                                  job_id)
        r = requests.get(url_status)
        r.raise_for_status()

        status = json.loads(r.text)['status']
        print("Job status is: '{}'".format(status))

        # If the status equals SUCCESS, exit out of the loop by changing the
        # condition ready. This causes the code to drop into the `else` block
        # below.
        #
        # If the status equals either FAILURE or REVOKED, an exception is raised
        # containing the error message. The program exits.
        #
        # Otherwise, wait for five seconds and start at the beginning of the
        # loop again.
        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED']:
            raise Exception(json.loads(r.text)['message'])
        else:
            time.sleep(5)
    else:
        # Requests the result of the job. If an error occurs an exception is
        # raised and the program exits. If the request is successful, the result
        # is returned.
        url_result = '{}api/result/pdb_file/dssp/{}/'.format(rest_url,
                                                                  job_id)
        r = requests.get(url_result)
        r.raise_for_status()
        result = json.loads(r.text)['result']

        # Return the result to the caller, which prints it to the screen.
        return result

def b_site(f,a,b,c):
    fi = read_file(f)
    xa=[]
    xb=[]
    c = float(c)
    for i in fi:
             aa=i.split()
             if aa[0]=='ATOM' and aa[4]==a:
                 xa.append(aa)
             if aa[0] == 'MODEL' and aa[1] == '2':
                 break

    for i in fi:
             bb=i.split()
             if bb[0]=='ATOM' and bb[4]==b:
                 xb.append(bb)
             if bb[0] == 'MODEL' and aa[1] == '2':
                 break
    
    aaa=[]
    for i in xa:
        for j in xb:
            d = distance_formula(i,j)
            if d**(0.5) <= c and [j,i,round(d**0.5,2)] not in aaa:
                aaa.append([i,j,round(d**0.5,2)])
    
    return aaa





@app.route("/",methods=['POST','GET'])
def first():
     if request.method=='POST':
       f = request.form.get('Pdb')
       cut= request.form.get('cutoff')
       title = request.form.get('title')
       if title not in ['residue contact','interactions','binding_site']:
        return redirect(url_for(title,file = f,c=cut))
       else:
        return render_template('layout.html',x = title,y=f,z=cut)
       #return redirect(url_for('home',titles=title))
     return render_template('first.html')



@app.route("/home/")
@app.route("/home/",methods=['POST','GET'])
def home():
    
    if request.method=='POST':
      
      f = request.form.get('Pdb')
      #chain = request.form.get('chain')
      cutoff = request.form.get('cutoff')
      #ato = request.form.get('atom')
      reg = request.form.get('reg')
      ato = 'CA'
      
#       cpi=c.Cation_pi()
      base = (request.form.getlist('redirect'))
      
      
      
      #for i in base:
    try:     
      if base[0] in ['short range','medium range','long range','contact order','long range order','surrounding hydrophobicity','contacts']:
                if reg:
                #return redirect(url_for(i)+str(f)+'/'+str(cutoff[0])+'/'+str(cutoff[1]))
                  if cutoff:
                    return redirect(url_for('short_range',file=f,titles=base,cut =float(cutoff),reg = reg,ato = ato))
                  else:
                    return redirect(url_for('short_range',file=f,titles=base,cut =8,reg = reg,ato = ato))
                else:
                  if cutoff:
                    return redirect(url_for('short_range',file=f,titles=base,cut =float(cutoff),reg = 'all',ato = ato))
                  else:
                    return redirect(url_for('short_range',file=f,titles=base,cut=8,reg = 'all',ato = ato))
                    
      elif base[0] in ['ionic interaction', 'disulphide interaction', 'hydrophobic interaction', 'cation pi interaction', 'aromatic aromatic interaction', 'mmch', 'msch', 'ssch']:
                return redirect(url_for('medium_range',file=f,titles=base,reg = reg))
               
      #elif base[0] == 'contact map':
                #return redirect(url_for('contact_map',file=f))
      
      #elif base[0] == 'contact network':
                #return redirect(url_for('contact_network',file=f))
    except:
                b_s = (request.form.getlist('b_s'))
                
                return redirect(url_for('binding_site',file=f,cutoff=b_s[2],chain1=b_s[0],chain2=b_s[1]))


      #elif base[0] == 'ramachandran plot':
                #return redirect(url_for('ramachandran_plot',file=f))

    return render_template('layout.html')






#The properties
@app.route("/sr/")
@app.route("/sr/")
def short_range():
          t = request.args.getlist('titles')

          f = request.args['file']
          cu = request.args['cut']
          region = request.args['reg']
          ato = request.args['ato']
          sequ = ['']
          for i in Regions(f)[3]:
            if i !='\n':
               sequ.append(i)
          print((sequ))
          return render_template('sr.html',x=[(Short_range(f,cu),'short range'),(Medium_range(f,cu),'medium range'),(Long_range(f,cu),'long range'),(Contact_Order(f),'contact order'),(LRO(f),'long range order'),(Surr_Hydrophob(f),'surrounding hydrophobicity'),(A8(f,cu,ato),'contacts',cu)],a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2],d=sequ,title = t,r= region)
          
            
@app.route("/mr/")
@app.route("/mr/")
def medium_range():
          

        t = request.args.getlist('titles')

        f = request.args['file']
        #region = request.args['reg']
        #region = 'all'
        sequ = ['']
        for i in Regions(f)[3]:
            if i !='\n':
              sequ.append(i)
        
        return render_template('mr.html',x=[(Ionic_Interaction(f),'ionic interaction'),(Hydrophobic_interaction(f),'hydrophobic interaction'),(Cation_pi(f),'cation pi interaction'),(Disulphide_interaction(f),'disulphide interaction'),(Aromatic_Aromatic(f),'aromatic aromatic interaction'),('A')],a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2],d=sequ,title = t)


@app.route("/cmap/")
@app.route("/cmap/")
def contact_map():
     f = request.args['file']
     sequ = [' ']
     for i in Regions(f)[3]:
        
        if i !='\n':
            sequ.append(i)
     
     cm = c_m(f)
     cm = pd.DataFrame({'Amino acid':cm[0],' Amino acid':cm[1],'dist':cm[2]})
     fig = px.scatter(cm,x='Amino acid',y=' Amino acid',color ='dist',opacity=0.3,color_continuous_scale=px.colors.sequential.Bluyl,title='Contact Map     {}'.format(f),template='plotly_dark',width=750,height=550)
     fig.update_layout(
    font_family="Courier New",
    font_color="#66ff00",
    font_size=25,
    title_font_size=25,
    legend_title_font_color="green",
    
)
     my_plot_div = fig.to_html(full_html=True)
     #my_plot_div = plot([Scatter(x=cm[0], y=cm[1])],x='Amino_acid',y='Amino_acid', output_type='div')
     return render_template('c_map.html',div_placeholder=Markup(my_plot_div),f=f,a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2],d=sequ)


@app.route("/cnet/")
@app.route("/cnet/")
def contact_network():
     f = request.args['file']
     sequ = [' ']
     for i in Regions(f)[3]:
        
        if i !='\n':
            sequ.append(i)
     cm = c_n(f)
     #cm = pd.DataFrame({'x':cm[0],'y':cm[1],'z':cm[2]})
     fig = go.Figure(data=go.Scatter3d(
    x=cm[0], y=cm[1], z=cm[2],
    marker=dict(
        size=10,
        color=cm[2],
        colorscale='Viridis',
    ),
    line=dict(
        color='darkblue',
        width=1.5
    )
))
     fig.update_layout(
    width=800,
    height=700,
    autosize=False,
    scene=dict(
        camera=dict(
            up=dict(
                x=0,
                y=0,
                z=0
            ),
            eye=dict(
                x=0,
                y=1.0707,
                z=1,
            )
        ),
        aspectratio = dict( x=1, y=1, z=0.7 ),
        aspectmode = 'manual'
    ),
)
     my_plot_div = fig.to_html(full_html=True)
     #my_plot_div = plot([Scatter(x=cm[0], y=cm[1])],x='Amino_acid',y='Amino_acid', output_type='div')
     return render_template('c_net.html',div_placeholder=Markup(my_plot_div),f=f,a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2],d=sequ)



@app.route("/b_site/")
@app.route("/b_site/")
def binding_site():
          f = request.args['file']
          cut = request.args['cutoff']
          chain1 = request.args['chain1']
          chain2 = request.args['chain2']
          
          result = b_site(f,chain1,chain2,cut)
          sequ = ['']
          for i in Regions(f)[3]:
            if i !='\n':
              sequ.append(i)

          return render_template('b_site.html',x=result,a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2],d=sequ)



@app.route("/co/")
@app.route("/co/<f>")
def contact_order(f):
          
          return render_template('co.html',y=Contact_Order(f))




@app.route("/sh/")
@app.route("/sh/<f>")
def surr_hydro(f):
          
          return render_template('surr.html',y=Surr_Hydrophob(f),z='Surrounding_hydrophobicity')




@app.route("/lro/")
@app.route("/lro/<f>")
def long_range_order(f):
          stack = LRO(f)
          return render_template('ion.html',y=stack[0],x=stack[2],z='Long_Range_Order')



@app.route("/cd/")
@app.route("/cd/<f>")
def contact_degree(f):
          stack = Contact_Degree(f)
          return render_template('cd.html',y=stack,z='Contact_Degree')




@app.route("/cs/")
@app.route("/cs/<f>/<coff1>/<coff2>")
def contacts(f,coff1,coff2):
    
    return render_template('cs.html',z='Total No of Contacts',a=(coff1),b=A8(f,coff1),c=(coff2),d=B8(f,coff2),x='CA atoms',y='CB atoms')





@app.route("/ion/")
@app.route("/ion/<f>")
def ionic_interaction(f):
        return render_template('ion.html',y=Ionic_Interaction(f),z='Ionic_Interactions',a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2])

@app.route("/di/")
@app.route("/di/<f>")
def disulphide_interaction(f):
        return render_template('ion.html',y=Disulphide_interaction(f),z='Disulphide_Interactions',a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2])
    
@app.route("/hi/")
@app.route("/hi/<f>")
def hydrophobic_interaction(f):
        return render_template('ion.html',y=Hydrophobic_interaction(f),z='Hydrophobic_Interactions',a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2])
    

@app.route("/arm_arm/")
@app.route("/arm_arm/<f>")
def aromatic_aromatic(f):
        return render_template('ion.html',y=Aromatic_Aromatic(f),z='Aromatic_Aromatic_Interactions',a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2])


@app.route("/cpi/")
@app.route("/cpi/<f>")
def cation_pi(f):

        return render_template('ion.html',y=Cation_pi(f),z='Cation_pi_Interactions',a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2])
    

@app.route("/arm_s/")
@app.route("/arm_s/<f>")
def aromatic_sulphur(f):

        return render_template('ion.html',y=Aromatic_sulphur(f),z='Aromatic_Sulphur_Interactions',a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2])
    

@app.route("/mmch/")
@app.route("/mmch/<f>")
def mmch(f):
        return render_template('ion.html',y=Main_Main_H_Bond(f),z='Main_Main_Chain_Hydrogen_Interactions')


@app.route("/msch/")
@app.route("/msch/<f>")
def msch(f):
        return render_template('ion.html',y=Main_Side_Chain_H_Bond(f),z='Main_Side_Chain_Hydrogen_Interactions')


@app.route("/ssch/")
@app.route("/ssch/<f>")
def ssch(f):
        return render_template('ion.html',y=Side_Side_Chain_H_Bond(f),z='Side_Side_Chain_Hydrogen_Interactions')



@app.route("/re/")
@app.route("/re/<f>")
def regions(f):
    return render_template('region.html',a=Regions(f)[0],b=Regions(f)[1],c=Regions(f)[2],d=Regions(f)[3])


# @app.route("/di/")
# @app.route("/di/<f>")
# def disulphide_interaction(f):
#         return render_template('ion.html',y=Ionic_Interaction(f))
    
    
# @app.route("/di/")
# @app.route("/di/<f>")
# def disulphide_interaction(f):
#         return render_template('ion.html',y=Ionic_Interaction(f))
    


     


        

        
      
      
    # return render_template('layout.html')



if __name__ == "__main__":
    app.run(debug=True)
    
  
    
    
