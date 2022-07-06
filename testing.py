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
    file = open("C:/Users/MOHITHTHAMANA/Desktop/Transmembrane/{}.pdb".format(f))
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

# @app.route("/ion/")
# def ionic_interaction():
#         return render_template('ion.html',y=)




     


        

        
      
      
    # return render_template('layout.html')



#if __name__ == "__main__":
app.run(debug=True)
    
  
    
    
