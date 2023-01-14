"""
Balls collision simulation

@authors: 
    Juan Jose Potes Gomez
    Julie Alejandra Ibarra
    Cristian Camilo Jimenez
"""

import pygame
from pygame.locals import *

from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np
import random
import copy

# Valores para los ejes en el plano 2D
ejex = 150
ejey = 100

# Coeficientes de restitucion
e_1 = 0.9
e_2 = 0.95
e_3 = 0.97

# Gravedad en cm/s
g = 981
                    
# Tiempo total que debe durar la simulacion
t_max = 120

# Declaracion de atributos de las particulas
rad_p = 3
pos_p = np.array([[-15.0 + rad_p + 0.1],[50.0]])
masa_p = 0.002 #Kg
ang_p = 60.0 * np.pi / 180 #rad
vel_p = 250 #cm/s
a_p = np.array([[0],[-g]])

n_particulas = 6
t_sal = 0
particulas = []

# Declaracion de atributos del hexagono particulas
rad_h = 3
pos_h = np.array([[0.0],[30.0]])
masa_h = 0.5 #Kg
ang_h = 0
vel_h = np.array([[0.0],[0.0]])
a_h = np.array([[0],[0]])

# Angulos de cada linea del hexagono
angulos = np.array([[120.0],[180.0],[240.0],[300.0],[0.0],[60.0]])
for i in range (0, len(angulos)):
    v = float(angulos[i])
    angulos[i] = v *np.pi / 180

# Variacion del tiempo
ht = 0.002
t = 0

# Funcion para pasar angulo de radianes a grados y viceversa
def conv_ang(angulo, tipo):
    if tipo == "rad":
        conv = (angulo * np.pi)/180
        return conv
    if tipo == "grad":
        conv = (angulo * 180)/np.pi
        return conv
        

# Funcion para pintar circulo en OpenGL
def circle(nsides, xc, yc, radio, r, g, b):
    n = 0
    glBegin(GL_LINE_STRIP)
    glColor3f(r, g, b)
    while(n <= nsides):
        angle = 2 * np.pi * n / nsides
        x = xc + radio * np.cos(angle)
        y = yc + radio * np.sin(angle)
        glVertex2f(x, y)
        n += 1
    glEnd()

# Funcion para pintar linea en OpenGL
def linea(x1, y1, x2, y2, r, g, b):
    glBegin(GL_LINES)
    glColor3f(r, g, b)
    glVertex2f(x1, y1)
    glVertex2f(x2, y2)
    glEnd()

# Funcion para pintar los ejes del plano de referencia
def plano():
    linea(0,-ejey,0,ejey,0.8,0.8,0.8)
    linea(-ejex,0,ejex,0,0.8,0.8,0.8)

# Funciones para metodo de Runge Kutta 4to orden
def f_ax(v,p,tiempo):
    resul = 0
    return resul

def f_ay(v,p,tiempo):
    resul = -g
    return resul

def f_vx(v,p,tiempo):
    resul = v
    return resul

def f_vy(v,p,tiempo):
    resul = v
    return resul

# Definicion de clase particula
class Particula:
    
    # Constructor de la clase
    def __init__(self, pos, rad, masa, ang, vel, acel, red, green, blue, tsalida):
        self.pos = pos
        self.rad = rad
        self.masa = masa
        self.ang = ang
        self.vel = vel
        self.v = np.array([[self.vel*np.cos(self.ang)],[self.vel*np.sin(self.ang)]])
        self.acel = acel
        self.red = red
        self.green = green
        self.blue = blue
        self.tsalida = tsalida
    
    # Metodo para mostrar la particula en la posicion respectiva
    def graficar(self):
        circle(20, self.pos[0],self.pos[1],self.rad,self.red,self.green,self.blue)
        # p_dx = self.pos[0] + self.v[0]/10
        # p_dy = self.pos[1] + self.v[1]/10
        # linea(self.pos[0],self.pos[1],p_dx,p_dy,self.red,self.green,self.blue)
        
        
    # Metodo que hace los calculos para determinar la nueva posicion con respecto a un intervalo de tiempo
    def avanzar(self,ht):
        global t
        
        #Runge Kutta 4to orden
        
        # Datos en Y
        v_y = self.v[1]
        y = self.pos[1]
        
        k1 = f_vy(v_y,y,t)
        l1 = f_ay(v_y,y,t)
        k2 = f_vy(v_y+(ht*l1/2),y+(ht*k1/2),t+(ht/2))
        l2 = f_ay(v_y+(ht*l1/2),y+(ht*k1/2),t+(ht/2))
        k3 = f_vy(v_y+(ht*l2/2),y+(ht*k2/2),t+(ht/2))
        l3 = f_ay(v_y+(ht*l2/2),y+(ht*k2/2),t+(ht/2))
        k4 = f_vy(v_y+(ht*l3),y+(ht*k3),t+ht)
        l4 = f_ay(v_y+(ht*l3),y+(ht*k3),t+ht)
        
        self.pos[1] = y + ((ht/6)*(k1 + 2*k2 + 2*k3 + k4))
        self.v[1] = v_y + ((ht/6)*(l1 + 2*l2 + 2*l3 + l4))
        
        # Datos en X
        v_x = self.v[0]
        x = self.pos[0]
        
        k1 = f_vx(v_x,x,t)
        l1 = f_ax(v_x,x,t)
        k2 = f_vx(v_x+(ht*l1/2),x+(ht*k1/2),t+(ht/2))
        l2 = f_ax(v_x+(ht*l1/2),x+(ht*k1/2),t+(ht/2))
        k3 = f_vx(v_x+(ht*l2/2),x+(ht*k2/2),t+(ht/2))
        l3 = f_ax(v_x+(ht*l2/2),x+(ht*k2/2),t+(ht/2))
        k4 = f_vx(v_x+(ht*l3),x+(ht*k3),t+ht)
        l4 = f_ax(v_x+(ht*l3),x+(ht*k3),t+ht)
        
        self.pos[0] = x + ((ht/6)*(k1 + 2*k2 + 2*k3 + k4))
        self.v[0] = v_x + ((ht/6)*(l1 + 2*l2 + 2*l3 + l4))

        self.vel = np.sqrt((self.v[0]**2)+(self.v[1]**2))
        self.ang = np.arctan(self.v[1]/self.v[0])
    
    # Metodo que realiza la colision entre 2 particulas
    def colision_part(self, part2):
        global e_1
        
        # Se calcula la distancia entre los centros de ambas particulas
        dist = np.sqrt((self.pos[0] - part2.pos[0]) ** 2 + (self.pos[1] - part2.pos[1]) ** 2)
        
        # Se evalua si hubo colision
        if dist <= (self.rad + part2.rad):
            
            # Angulo de choque
            angulo = np.arctan((part2.pos[1] - self.pos[1])/(part2.pos[0] - self.pos[0]))
            
            # Componentes de velocidad antes del choque en el plano rotado (Matriz de rotacion) para la particula 1
            vp1_ac = (self.v[0] * np.cos(angulo)) + (self.v[1] * np.sin(angulo))
            vn1 = - (self.v[0] * np.sin(angulo)) + (self.v[1] * np.cos(angulo))
            
            # Componentes de velocidad antes del choque en el plano rotado (Matriz de rotacion) para la particula 2
            vp2_ac = (part2.v[0] * np.cos(angulo)) + (part2.v[1] * np.sin(angulo))
            vn2 = - (part2.v[0] * np.sin(angulo)) + (part2.v[1] * np.cos(angulo))
            
            m1 = self.masa
            m2 = part2.masa
            
            # Se calculan las velocidades despues del choque (Solo se ve afectada la componente perpendicular al angulo de choque)
            vp1_dc = (((m1 - (e_1 * m2))/(m1 + m2)) * vp1_ac) + ((((1 + e_1) * m2)/(m1 + m2)) * vp2_ac)
            vp2_dc = ((((1 + e_1) * m1)/(m1 + m2)) * vp1_ac) + (((m2 - (e_1 * m1))/(m1 + m2)) * vp2_ac)
            
            # Componentes de velocidad despues del choque en el plano original (Matriz de rotacion inversa) para la particula 1
            self.v[0] = (vp1_dc * np.cos(angulo)) - (vn1 * np.sin(angulo))
            self.v[1] = (vp1_dc * np.sin(angulo)) + (vn1 * np.cos(angulo))
            
            # Magnitud y direccion de la velocidad
            self.vel = np.sqrt((self.v[0] ** 2) + (self.v[1] ** 2))
            self.ang = np.arctan(self.v[1] / self.v[0])
            
            # Componentes de velocidad despues del choque en el plano original (Matriz de rotacion inversa) para la particula 2
            part2.v[0] = (vp2_dc * np.cos(angulo)) - (vn2 * np.sin(angulo))
            part2.v[1] = (vp2_dc * np.sin(angulo)) + (vn2 * np.cos(angulo))
            
            # Magnitud y direccion de la velocidad
            part2.vel = np.sqrt((part2.v[0] ** 2) + (part2.v[1] ** 2))
            part2.ang = np.arctan(part2.v[1] / part2.v[0])
            
            # Correccion de posicion para que no se siga chocando
            if(part2.pos[0] >= self.pos[0]):
                part2.pos[0] += 0.1
                self.pos[0] -= 0.1
            else:
                part2.pos[0] -= 0.1
                self.pos[0] += 0.1
            
            if(part2.pos[1] >= self.pos[1]):
                part2.pos[1] += 0.1
                self.pos[1] -= 0.1
            else:
                part2.pos[1] -= 0.1
                self.pos[1] += 0.1
            
            
# Definicion de la clase Hexagono
class Hexagono():
    # Constructor de la clase
    def __init__(self, pos, rad, masa, ang, vel, acel, red, green, blue):
        self.pos = pos
        self.rad = rad
        self.masa = masa
        self.ang = ang
        self.v = vel
        self.acel = acel
        self.red = red
        self.green = green
        self.blue = blue
    
    # Metodo para mostrar la particula en la posicion respectiva
    def graficar(self):
        circle(6, self.pos[0], self.pos[1], self.rad, self.red, self.green, self.blue)
        # p_dx = self.pos[0] + self.v[0]/10
        # p_dy = self.pos[1] + self.v[1]/10
        # linea(self.pos[0],self.pos[1],p_dx,p_dy,self.red,self.green,self.blue)
        
    # Metodo que hace los calculos para determinar la nueva posicion con respecto a un intervalo de tiempo
    def avanzar(self,ht):
        global t
        
        # Datos en X
        self.pos[0] += self.v[0] * ht
    
    
    # Metodo que realiza la colision entre el hexagono y una particula
    def colision_h_part(self, part2):
        global e_2, angulos, angulos_p
        n_lado = 0
        tol = 0.1
        choque_esquina = False
        # Se calcula la distancia entre los centros de ambas particulas
        dist = np.sqrt((self.pos[0] - part2.pos[0]) ** 2 + (self.pos[1] - part2.pos[1]) ** 2)
        
        # Se verifica si pudo haber colision (Si se choca con el circulo que rodea el hexagono)
        if dist <= (self.rad + part2.rad):
            
            # Angulo de choque
            angulo = np.arctan((part2.pos[1] - self.pos[1])/(part2.pos[0] - self.pos[0]))
            angulo_g = conv_ang(angulo,"grad")
            
            # Si se choca con los vertices
            if((abs(angulo_g - 60) <= tol) or (abs(angulo_g + 60) <= tol) or (abs(part2.pos[1]-self.pos[1]) <= tol)):
                choque_esquina = True
            
            # Se evalua con que lado choca
            if(part2.pos[0] > (self.pos[0] + self.rad * np.cos(conv_ang(60,"rad")))):
                if(part2.pos[1] > self.pos[1]):
                    n_lado = 0
                elif(part2.pos[1] < self.pos[1]):
                    n_lado = 5
            elif(part2.pos[0] < (self.pos[0] + self.rad * np.cos(conv_ang(120,"rad")))):
                if(part2.pos[1] > self.pos[1]):
                    n_lado = 2
                elif(part2.pos[1] < self.pos[1]):
                    n_lado = 3
            elif((self.pos[0] + self.rad * np.cos(conv_ang(120,"rad"))) < part2.pos[0] < (self.pos[0] + self.rad * np.cos(conv_ang(60,"rad")))):
                if(part2.pos[1] > self.pos[1]):
                    n_lado = 1
                elif(part2.pos[1] < self.pos[1]):
                    n_lado = 4
            else:
                choque_esquina = True
            
            # Si no se choca con los vertices
            if(choque_esquina == False):
                # Se sacan los 2 vertices del hexagono que forman cada lado
                menor = 2 * np.pi * n_lado / 6
                mayor = 2 * np.pi * (n_lado + 1)/ 6
                esq1 = np.array([[self.pos[0] + self.rad * np.cos(menor)],[self.pos[1] + self.rad * np.sin(menor)]])
                esq2 = np.array([[self.pos[0] + self.rad * np.cos(mayor)],[self.pos[1] + self.rad * np.sin(mayor)]])
                
                # Se halla la pendiente y el punto de corte de la recta formada por los 2 vertices
                deno = esq2[0] - esq1[0]
                m = (esq2[1] - esq1[1])/deno
                b = - m * esq1[0] + esq1[1]
                
                # Formula de distancia de un punto a una recta
                dist2 = abs((m * part2.pos[0]) - part2.pos[1] + b)/(np.sqrt((m ** 2) + 1))
                
                # Se verifica si colisiona con la recta de ese lado del hexagono
                if(dist2 <= part2.rad):
                    # Angulo perpendicular al angulo de la linea
                    angulo3 = angulos[n_lado] + (np.pi/2)
                    
                    # Componentes de velocidad del hexagono en el plano rotado 
                    vp1_ac = (self.v[0] * np.cos(angulo3)) + (self.v[1] * np.sin(angulo3))
                    vn1 = - (self.v[0] * np.sin(angulo3)) + (self.v[1] * np.cos(angulo3))
                    
                    # Componentes de velocidad de la particula en el plano rotado 
                    vp2_ac = (part2.v[0] * np.cos(angulo3)) + (part2.v[1] * np.sin(angulo3))
                    vn2 = - (part2.v[0] * np.sin(angulo3)) + (part2.v[1] * np.cos(angulo3))
                    
                    m1 = self.masa
                    m2 = part2.masa
                    
                    # Se calculan las velocidades despues del choque (Solo cambia la velocidad perpendicular al plano)
                    vp1_dc = (((m1 - (e_2 * m2))/(m1 + m2)) * vp1_ac) + ((((1 + e_2) * m2)/(m1 + m2)) * vp2_ac)
                    vp2_dc = ((((1 + e_2) * m1)/(m1 + m2)) * vp1_ac) + (((m2 - (e_2 * m1))/(m1 + m2)) * vp2_ac)
                    
                    # Se pasan los componentes de velocidad del hexagono al plano normal (En el caso del hexagono, no se tiene en cuenta la velocidad en Y)
                    self.v[0] = (vp1_dc * np.cos(angulo3)) - (vn1 * np.sin(angulo3))
                    self.v[1] = 0
                    
                    # Magnitud y direccion de la velocidad
                    self.vel = np.sqrt((self.v[0] ** 2) + (self.v[1] ** 2))
                    self.ang = np.arctan(self.v[1] / self.v[0])
                    
                    # Se pasan los componentes de velocidad de la particula al plano normal (Rotacion inversa)
                    part2.v[0] = (vp2_dc * np.cos(angulo3)) - (vn2 * np.sin(angulo3))
                    part2.v[1] = (vp2_dc * np.sin(angulo3)) + (vn2 * np.cos(angulo3))
                    
                    # Magnitud y direccion de la velocidad
                    part2.vel = np.sqrt((part2.v[0] ** 2) + (part2.v[1] ** 2))
                    part2.ang = np.arctan(part2.v[1] / part2.v[0])
                    
                    # Correccion de posicion para que no se siga chocando
                    if(part2.pos[0] >= self.pos[0]):
                        part2.pos[0] += 0.1
                        self.pos[0] -= 0.1
                    else:
                        part2.pos[0] -= 0.1
                        self.pos[0] += 0.1
                    
                    if(part2.pos[1] >= self.pos[1]):
                        part2.pos[1] += 0.1
                    else:
                        part2.pos[1] -= 0.1
            else:
                # Si se choca con los vertices, el angulo de choque es el mismo angulo entre los centros
                
                # Componentes de velocidad antes del choque en el plano rotado (Matriz de rotacion) para la particula 1
                vp1_ac = (self.v[0] * np.cos(angulo)) + (self.v[1] * np.sin(angulo))
                vn1 = - (self.v[0] * np.sin(angulo)) + (self.v[1] * np.cos(angulo))
                
                # Componentes de velocidad antes del choque en el plano rotado (Matriz de rotacion) para la particula 2
                vp2_ac = (part2.v[0] * np.cos(angulo)) + (part2.v[1] * np.sin(angulo))
                vn2 = - (part2.v[0] * np.sin(angulo)) + (part2.v[1] * np.cos(angulo))
                
                m1 = self.masa
                m2 = part2.masa
                
                # Se calculan las velocidades despues del choque (Solo se ve afectada la componente perpendicular al angulo de choque)
                vp1_dc = (((m1 - (e_1 * m2))/(m1 + m2)) * vp1_ac) + ((((1 + e_1) * m2)/(m1 + m2)) * vp2_ac)
                vp2_dc = ((((1 + e_1) * m1)/(m1 + m2)) * vp1_ac) + (((m2 - (e_1 * m1))/(m1 + m2)) * vp2_ac)
                
                # Componentes de velocidad despues del choque en el plano original (Matriz de rotacion inversa) para la particula 1
                self.v[0] = (vp1_dc * np.cos(angulo)) - (vn1 * np.sin(angulo))
                self.v[1] = 0
                
                # Magnitud y direccion de la velocidad
                self.vel = np.sqrt((self.v[0] ** 2) + (self.v[1] ** 2))
                self.ang = np.arctan(self.v[1] / self.v[0])
                
                # Componentes de velocidad despues del choque en el plano original (Matriz de rotacion inversa) para la particula 2
                part2.v[0] = (vp2_dc * np.cos(angulo)) - (vn2 * np.sin(angulo))
                part2.v[1] = (vp2_dc * np.sin(angulo)) + (vn2 * np.cos(angulo))
                
                # Magnitud y direccion de la velocidad
                part2.vel = np.sqrt((part2.v[0] ** 2) + (part2.v[1] ** 2))
                part2.ang = np.arctan(part2.v[1] / part2.v[0])
                
                # Correccion de posicion para que no se siga chocando
                if(part2.pos[0] >= self.pos[0]):
                    part2.pos[0] += 0.1
                    self.pos[0] -= 0.1
                else:
                    part2.pos[0] -= 0.1
                    self.pos[0] += 0.1
                
                if(part2.pos[1] >= self.pos[1]):
                    part2.pos[1] += 0.1
                else:
                    part2.pos[1] -= 0.1
                        
hxgn = Hexagono(pos_h, rad_h, masa_h, ang_h, vel_h, a_h, 255, 255, 255)

# Funcion que grafica el recipiente
def mostrar_recipiente():
    linea(-5,0,15,0,255,255,255)
    linea(15,0,15,ejey,255,255,255)
    linea(-15,10,-15,ejey,255,255,255)
    linea(-15,10,-5,0,255,0,255)

# Funcion que verifica la colision entre el recipiente y una particula o hexagono
def colision_recipiente(part, tipo):
    global e_3
    
    # Si es particula o hexagono
    if(tipo == 1 or tipo == 2):
        # Verifica colision con pared izquierda
        if((abs(part.pos[0]-(-15)) <= part.rad)):
            # Nueva velocidad en X
            part.v[0] = - e_3 * part.v[0]
            while(abs(part.pos[0]-(-15)) <= part.rad):
                part.pos[0] += 0.01
        
        # Verifica colision con pared derecha
        if(abs(part.pos[0]-(15)) <= part.rad):
            # Nueva velocidad en X
            part.v[0] = - e_3 * part.v[0]
            while(abs(part.pos[0]-(15)) <= part.rad):
                part.pos[0] -= 0.01
    
    # Solo si es particula
    if(tipo == 1):
        # Verifica colision con el piso
        if(part.pos[1] <= part.rad):
            part.v[1] = - e_3 * part.v[1]
            while(part.pos[1] <= part.rad):
                part.pos[1] += 0.01
        
        # Verifica colision con recta inclinada
        # Pendiente y punto de corte de la recta del piso
        m = (10 - 0)/(-15 - (-5))
        b = - m * (-5) + 0
        
        # Angulo del plano
        d_ang = np.arctan(m)
        
        # Formula de distancia de un punto a una recta
        dist = abs((m*part.pos[0])-part.pos[1]+b)/(np.sqrt((m**2)+1))
        
        # Si se choca con la recta
        if(dist <= part.rad):
            # Angulo perpendicular al de la recta
            angulo = d_ang + (np.pi/2)
            
            # Componentes velocidad de la particula rotados
            vp1_ac = (part.v[0] * np.cos(angulo)) + (part.v[1] * np.sin(angulo))
            vn1 = - (part.v[0] * np.sin(angulo)) + (part.v[1] * np.cos(angulo))
            
            # Nuevo componente de velocidad perpendicular al plano
            vp1_dc = - e_3 * vp1_ac + 100
            
            # Nuevos componentes de velocidad en el plano normal (Rotacion inversa)
            part.v[0] = (vp1_dc * np.cos(angulo)) - (vn1 * np.sin(angulo)) 
            part.v[1] = (vp1_dc * np.sin(angulo)) + (vn1 * np.cos(angulo)) 
            
            # Magnitud y direccion de la velocidad
            part.vel = np.sqrt((part.v[0] ** 2) + (part.v[1] ** 2))
            part.ang = np.arctan(part.v[1] / part.v[0])
            
            dist2 = abs((m*part.pos[0])-part.pos[1]+b)/(np.sqrt((m**2)+1))
            if (dist2 <= part.rad):
                while(dist2 <= part.rad):
                    part.pos[0] += 0.01
                    part.pos[1] += 0.01
                    dist2 = abs((m*part.pos[0])-part.pos[1]+b)/(np.sqrt((m**2)+1))



# Se declaran las particulas en una lista
while len(particulas) < n_particulas:   
    red = random.random() + 0.3
    green = random.random() + 0.3
    blue = random.random() + 0.3
    part_base = Particula(pos_p, rad_p, masa_p, ang_p, vel_p, a_p, red, green, blue, t_sal)
    particulas.append(copy.deepcopy(part_base))    
    t_sal += 0.5


# Funcion principal
def main():
    global t
    running = True
    pygame.init()
    display=(900,600)
    pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
    gluOrtho2D(-ejex/2,ejex/2,-10,ejey-10)
    
    # Ciclo para que se vaya visualizando la simulacion
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                running=False
                break
            
        # Se ejecuta mientras que el tiempo total sea menor al tiempo maximo definido
        if( running == True and t <= t_max):
            # Se limpia la pantalla de OpenGL
            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
            
            # Se grafican el hexagono, las particulas y el recipiente
            hxgn.graficar()
            for i in range (0, len(particulas)):
                if(t >= particulas[i].tsalida):
                    particulas[i].graficar()
            mostrar_recipiente()
            
            # Se verifica si hay colision con el recipiente para todas las particulas y el hexagono
            for i in range (0, len(particulas)):
                colision_recipiente(particulas[i], 1)
            colision_recipiente(hxgn,2)
            
            # Se verifica si hay colision entre particulas
            for i in range (0, len(particulas)):
                for j in range(i+1 , len(particulas)):
                    if(t >= particulas[i].tsalida and t >= particulas[j].tsalida):
                        particulas[i].colision_part(particulas[j])
            
            # Se verifica si hay colision entre particulas y hexagono
            for i in range (0, len(particulas)):
                hxgn.colision_h_part(particulas[i])
            
            # Avanzan las particulas y el hexagono
            for i in range (0, len(particulas)):
                if(t >= particulas[i].tsalida):
                    particulas[i].avanzar(ht)
            hxgn.avanzar(ht)
            
            pygame.display.flip() #Mostrar pantalla
            pygame.time.wait(1)
            print("Tiempo procesado: ",round(t,3)," s")
            t+=ht
        else:
            running = False
            break
    print("Tiempo total simulado = ",round(t,1))
    pygame.quit()

main()


