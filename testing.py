# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:10:15 2022

@author: MOHITHTHAMANA
"""

from flask import Flask, redirect, url_for, render_template, request

import os
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

def distance_formula(i,j):
        i[6]=float(i[6])
        i[7]=float(i[7])
        i[8]=float(i[8])
        j[6]=float(j[6])
        j[7]=float(j[7])
        j[8]=float(j[8])
        d=(i[6]-j[6])**2+(i[7]-j[7])**2+(i[8]-j[8])**2
        return d

def Short_range(f):
        i1=0
        fi = read_file(f)
        xa=[]
        sr=[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        for i in xa:
            j1=0
            for j in xa:
                d=distance_formula(i,j)
                if 3.5<d**0.5<7 and [j[3],j[5],i[3],i[5],round(d**0.5,2)] not in sr and abs(i1-j1)<=2 :
                    sr.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return sr
    
def Medium_range(f):
        i1=0
        fi = read_file(f)
        xa=[]
        mr=[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        for i in xa:
            j1=0
            for j in xa:
                d=distance_formula(i,j)
                if 3.5<d**0.5<7 and [j[3],j[5],i[3],i[5],round(d**0.5,2)] not in mr and abs(i1-j1)<=2 :
                    mr.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return mr
    

def Long_range(f):
        i1=0
        fi = read_file(f)
        xa=[]
        lr=[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        for i in xa:
            j1=0
            for j in xa:
                d=distance_formula(i,j)
                if 3.5<d**0.5<7 and [j[3],j[5],i[3],i[5],round(d**0.5,2)] not in lr and abs(i1-j1)<=2 :
                    lr.append([i[3],i[5],j[3],j[5],round(d**0.5,2)])
                j1=j1+1
            i1=i1+1
        return lr

def LRO(f):
        i1=0
        fi = read_file(f)
        xa=[]
        lro = []
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        
        for i in xa:
            j1=0
            s=0
            for j in xa:
                d=distance_formula(i,j)
                if d**0.5<8 and abs(i1-j1)>12 :
                    s=s+(1/len(xa))
                j1=j1+1
            lro.append(s)
            i1=i1+1
        return lro
    
    
    
def Contact_Order(f):
        i1=0
        fi = read_file(f)
        xa=[]
        c_order=0
        l=0
        summ=0
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CA':
                xa.append(aa)
            if aa[0] == 'MODEL' and aa[1] == '2':
                break
        for i in xa:
            j1=0
            c_order=0
            l=0
            for j in xa:
                d=distance_formula(i,j)
                if d**0.5<8:
                    l=l+1
                c_order=c_order+(abs(j1-i1))
                j1=j1+1
            summ=summ+(c_order/l)
            i1=i1+1
        return(summ/len(xa))

# def A8(self):
#         inp = float(input())
#         degree = {}
#         for i in self.xa:
#             j1=0
#             deg = 0
#             for j in self.xa:
#                 d=c.distance_formula(i,j)
#                 if 0<d**0.5<=inp:
#                     deg+=1
#                 j1=j1+1
#             degree['{0}-{1}'.format(i[3],i[5])] = deg
#         #print('C-alpha:{0}'.format(inp))
#         return(degree)
#     # def A14(self):
#     #     inp = float(input())
#     #     degree = {}
#     #     for i in self.xa:
#     #         j1=0
#     #         deg = 0
#     #         for j in self.xa:
#     #             d=c.distance_formula(i,j)
#     #             if 0<d**0.5<=inp:
#     #                 deg+=1
#     #             j1=j1+1
#     #         degree['{0}-{1}'.format(i[3],i[5])] = deg
#     #     print(degree)
# def B8(self):
#         inp = float(input())
#         degree = {}
#         for i in self.xb:
#             j1=0
#             deg = 0
#             for j in self.xb:
#                 d=c.distance_formula(i,j)
#                 if 0<d**0.5<=inp:
#                     deg+=1
#                 j1=j1+1
#             degree['{0}-{1}'.format(i[3],i[5])] = deg
#         #print('C-beta:{0}'.format(inp))
#         return (degree)
#     # def B14(self):
#     #     inp = float(input())
#     #     degree = {}
#     #     for i in self.xb:
#     #         j1=0
#     #         deg = 0
#     #         for j in self.xb:
#     #             d=c.distance_formula(i,j)
#     #             if 0<d**0.5<=inp:
#     #                 deg+=1
#     #             j1=j1+1
#     #         degree['{0}-{1}'.format(i[3],i[5])] = deg
#     #     print(degree)
def Disulphide_interaction(f):
        
        fi = read_file(f)
        cys1,cys = [],[]
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='SG':
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
                if 0<d**0.5<=2.2 and [j[3],j[5],i[3],i[5],d**0.5] not in p:
                    p.append([i[3],i[5],j[3],j[5],d**0.5])
                j1=j1+1
            i1=i1+1
        return(p)
def Hydrophobic_interaction(f):
        fi = read_file(f)
        p,xb=[],[]
        i1=0
        for i in fi:
            aa=i.split()
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CB':
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
                    p.append([i[3],i[5],j[3],j[5],d**0.5])
                j1=j1+1
            i1=i1+1
        return(p)
        

    
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
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['NH1','NH2','NZ','ND1','NE2'] and aa[3] in ['ARG','LYS','HIS']:
                cation.append(aa)
            if aa[0]=='ATOM' and aa[4]=='A' and aa[2] in ['OE2','OE1','OD1','OD2'] and aa[3] in ['ASP','GLU']:
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
    
# def Surr_Hydrophob(f):
        
#         fi = read_file(f)
#         xa=[]
        
#         for i in fi:
#             aa=i.split()
#             if aa[0]=='ATOM' and aa[4]=='A' and aa[2]=='CA':
#                 xa.append(aa)
#             if aa[0] == 'MODEL' and aa[1] == '2':
#                 break
#         hydro_index = {'ALA': 0.87 , 'ASP': 0.66, 'CYS': 1.52, 'GLU': 0.67, 'PHE': 2.87 , 'GLY':0.1, 'HIS': 0.87, 'ILE': 3.15,
#        			'LYS': 1.64, 'LEU': 2.17,'MET':1.67 ,'ASN': 0.09,'PRO': 2.77,'GLN': 0 ,'ARG':0.85, 'SER': 0.07,
#        			'THR':0.07, 'VAL':1.87, 'TRP': 3.77, 'TYR': 2.67}
        
#         surr = []
        
#         i1=0
#         for i in xa:
#             j1=0
#             interact = []
#             score = 0
#             for j in xa:
#                 d=distance_formula(i,j)
#                 if 0<d**0.5<8:
#                     interact.append(j[3])
                
#             for k in set(interact):
#                    if k=='AHIS' or k=='BHIS':
#                     score+=(interact.count(k)*0.87)
#                    else:
#                     score+=(interact.count(k)*hydro_index[k]) 

#             j1=j1+1
#             surr.append([i[3],i[5],round(score,2)])
                
#             i1=i1+1
         
#         return((surr))
# @app.route("/",methods=['POST','GET'])
# def first():
#     if request.method=='POST':
#       f = request.form.get('Pdb')
      
#       return render_template('layout.html')
#     return render_template('first.html')


@app.route("/",methods=['POST','GET'])
def home():
    if request.method=='POST':
      f = request.form.get('Pdb')
      
      

#       cpi=c.Cation_pi()
      data = (request.form.getlist('redirect'))
      for i in data:
            
            
            return redirect(url_for(i)+str(f))
    return render_template('layout.html')

@app.route("/sr/")
@app.route("/sr/<f>")
def short_range(f):
          
          return render_template('sr.html',x=Short_range(f))

@app.route("/mr/")
@app.route("/mr/<f>")
def medium_range(f):
          
          return render_template('mr.html',x=Medium_range(f))


@app.route("/lr/")
@app.route("/lr/<f>")
def long_range(f):
          
          return render_template('lr.html',x=Long_range(f))



@app.route("/co/")
@app.route("/co/<f>")
def contact_order(f):
          
          return render_template('co.html',y=Contact_Order(f))

@app.route("/ion/")
@app.route("/ion/<f>")
def ionic_interaction(f):
        return render_template('ion.html',y=Ionic_Interaction(f),z='Ionic_Interactions')

@app.route("/di/")
@app.route("/di/<f>")
def disulphide_interaction(f):
        return render_template('ion.html',y=Disulphide_interaction(f),z='Disulphide_Interactions')
    
@app.route("/hi/")
@app.route("/hi/<f>")
def hydrophobic_interaction(f):
        return render_template('ion.html',y=Hydrophobic_interaction(f),z='Hydrophobic_Interactions')
    

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
  
    
    
