import xdrlib
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import time


class Equilibrium:
    def __init__(self,alpha):
        self.alpha=alpha
          
    def yeq(self,x):
        return ((self.alpha*x)/(1+x*(self.alpha-1)))
    
    def xeq(self,y):
        return (y/(-self.alpha*y+self.alpha+y))
    
    def Eqplot(self):
        x=np.linspace(0,1,100)
        y=self.yeq(x)
        plt.plot(x,y)
        plt.plot(x,x)
        

class Feed:
    Ftotal=0
    VolatilTotal=0
    Xftotal=0
    AmountOfFeeds=0
    TaggedMinimumRs=[]
    MinimumRs=[]
    def __init__(self,F,Xf,q):
        if q==0 or q==1:
            q+=0.001
        self.F=F
        self.Xf=Xf
        self.q=q
        Feed.Ftotal+=F
        Feed.VolatilTotal+=Xf*F
        x=sp.symbols('x')
        if q>0:
            self.xIntersect=max(sp.solve((x**2)*(q*(Eq.alpha-1))+x*(q-Xf*(Eq.alpha-1)-Eq.alpha*(q-1))-Xf,x))
        else:
            self.xIntersect=min(sp.solve((x**2)*(q*(Eq.alpha-1))+x*(q-Xf*(Eq.alpha-1)-Eq.alpha*(q-1))-Xf,x))
        self.yIntersect=((self.q/(self.q-1))*self.xIntersect)-(self.Xf/(self.q-1))
        print(f'Feed de Xf {self.Xf} equilibrio en x= {self.xIntersect:.2f}, y= {self.yIntersect:.2f}')
        self.rmin=(Props.xD-self.yIntersect)/(self.yIntersect-self.xIntersect)
        Feed.AmountOfFeeds+=1
        Feed.TaggedMinimumRs.append([f'f{Feed.AmountOfFeeds}',self.rmin])
        Feed.MinimumRs.append(self.rmin)
    def yFeed(self,x):
            print(f'Recta Alimento: {(self.q/(self.q-1)):.3f} *Xn-1 + {-(self.Xf/(self.q-1)):.3f} ')
            return ((self.q/(self.q-1))*x)-(self.Xf/(self.q-1))
    
    def MinR(self):
        return min(Feed.MinimumRs)
    
    def Feedplot(self):
        if self.q<1:
                xmin=self.xIntersect
                xmax=self.Xf
        else:
                xmin=self.Xf
                xmax=self.xIntersect
        x=np.linspace(float(xmin),float(xmax),5)    
        y=self.yFeed(x)
        plt.plot(x,y)
        plt.annotate(f'{self.xIntersect:.2f} , {self.yIntersect:.2f}',
                        [self.xIntersect,self.yIntersect],xytext=(self.xIntersect,self.yIntersect+0.1),
                        textcoords='axes fraction',arrowprops=dict(facecolor='black',width=3,headlength=8),
                        horizontalalignment='right', verticalalignment='top')
    
class Properties:
    def __init__(self,xD,xR,Rfactor):
        self.xD=xD
        self.xR=xR
        self.Rfactor=Rfactor
        self.workingR=0
        self.D=0
        self.R=0
        self.L=0
        self.V=0
        self.DxD=0
    def workR(self,Rmin):
        self.workingR=self.Rfactor*Rmin
        return self.Rfactor*Rmin
    def MaterialBalance(self,FeedTotal,TotalVolatil):
        self.D=(FeedTotal*self.xR-TotalVolatil)/(-self.xD+self.xR)
        self.R=FeedTotal-self.D
        self.L=self.D*self.workingR
        self.V=self.L/(self.workingR/(self.workingR+1))
        self.DxD=self.D*self.xD

Feeds=[]


##### PROPIEDADES A CAMBIAR #####
Eq=Equilibrium(2.46)                #VOLATILIDAD RELATIVA
Props=Properties(0.98,0.134,1.5)    # xD, xR, Veces*Rminima (ratio Rmin -> R de operacion)

###Poner los Feeds 
f1=Feed(1.196,0.6,1.25)
f2=Feed(1.355,0.3,0.3)              #PONER MAS FEEDS SI HACE FALTA TAL QUE, POR EJEMPLO: f2=Feed(CAUDAL MOLAR TOTAL, xF , q)
Feeds.append(f1)
Feeds.append(f2)                    #Llamar a los feeds: f + un numero empezando por 1 : f1,f2,f3,f4...
#SUS FEEDS EXTRA AQUI               #POR CADA FEED ADICIONAL PONER Feeds.append(NOMBRE DEL FEED) 

############################




Eq.Eqplot()

for feed in Feeds:
    feed.Feedplot()
    
print(f'R minima: {f1.MinR():.3f} , R de trabajo: {Props.workR(f1.MinR()):.3f}')

plt.xlim(0,1)
plt.ylim(0,1)


def Sort(sortees):
    return sorted(sortees,key=lambda x:x[1])

SortedR=Sort(Feed.TaggedMinimumRs)  #FEEDS ORDENADOS

Props.MaterialBalance(Feed.Ftotal,Feed.VolatilTotal)

def Crossing(xF,q,xD,R):
    xiF = (xF/(q-1)+xD/(R+1))/(q/(q-1)-R/(R+1))
    yiF = R/(R+1)*xiF + xD/(R+1)
    plt.annotate(f'{xiF:.2f} , {yiF:.2f}',
                        [xiF,yiF],xytext=(xiF,xiF+0.1),
                        textcoords='axes fraction',arrowprops=dict(facecolor='black',width=3,headlength=8),
                        horizontalalignment='right', verticalalignment='top')
    return xiF


def RectaOperacionSuperior(R,xD,xiF):
    x_rect=np.linspace(xiF,xD,5)
    y_rect= R /(R+1)*x_rect + xD/(R+1)
    print(f'Recta de Operacion: {R /(R+1):.3f} *Xn-1 + {xD/(R+1):.3f} ')
    return y_rect,x_rect

class RectInterm:
    L=Props.L
    V=Props.V
    DxD=Props.D*Props.xD
    def __init__(self,q,F,Xf):    
        self.q=q
        self.F=F
        self.FxF=F*Xf
        RectInterm.L+=F*q
        RectInterm.V-=F*(1-q)
        RectInterm.DxD-=(F*Xf)
        self.L=RectInterm.L
        self.V=RectInterm.V
        self.DxD=RectInterm.DxD
    def RIntermediaPlot(inicio,final):
        x_int=np.linspace(final,inicio,5)
        y_int=((RectInterm.L/RectInterm.V)*x_int)+(RectInterm.DxD/RectInterm.V)
        print(f'Recta de Operacion: {(RectInterm.L/RectInterm.V):.3f} *Xn-1 + {(RectInterm.DxD/RectInterm.V):.3f}')
        return y_int,x_int


def IntCros(rA,feed):
    
    result= (-rA.DxD*feed.q+rA.DxD-rA.V*feed.Xf)/(rA.L*feed.q-rA.L-feed.q*rA.V)
    yresult=feed.yFeed(result)
    plt.annotate(f'{result:.2f} , {yresult:.2f}',
                        [result,yresult],xytext=(result,yresult+0.1),
                        textcoords='axes fraction',arrowprops=dict(facecolor='black',width=3,headlength=8),
                        horizontalalignment='right', verticalalignment='top')
    return result


    

def Dibujo():
    n=0
    Lineas=[]
    Intersecciones=[]
    Graficos=[]
    Int1=Crossing(eval(SortedR[n][0]).Xf,eval(SortedR[n][0]).q,Props.xD,Props.workingR)
    Intersecciones.append(Int1)
    y_rect,x_rect= RectaOperacionSuperior(Props.workR(eval(SortedR[n][0]).MinR()),Props.xD,float(Intersecciones[0]))
    plt.plot(x_rect,y_rect)
    while n< len(SortedR):        
        a=RectInterm(eval(SortedR[n][0]).q,eval(SortedR[n][0]).F,eval(SortedR[n][0]).Xf)
        Lineas.append(a)
        try:
            Intersecciones.append(IntCros(Lineas[n],eval(SortedR[n+1][0])))
        except:
            Intersecciones.append(Props.xR)
        y_int,x_int=RectInterm.RIntermediaPlot(float(Intersecciones[n]),float(Intersecciones[n+1]))
        Graficos.append([x_int,y_int])
        plt.plot(Graficos[n][0],Graficos[n][1])
        n+=1

plt.xticks(np.arange(0, 1.1, 0.1))
plt.yticks(np.arange(0, 1.1, 0.1))
Dibujo()
plt.show()


