"""
TP2 méthode LU sur des matrice tridiagonales
"""
import numpy as np


def TRIreduite(A):
    
    A = np.copy(A)
    n,m = np.shape(A)
    
    A_tilde = np.zeros((n,3))
    d0 = np.diag(A)
    d1 = np.diag(A, -1)
    d2 = np.diag(A, 1)
    
    c1 = np.insert(d1, 0, 0)
    c2 = np.append(d2,0)
    A_tilde = np.column_stack((c1, d0, c2))
    
    return A_tilde


def TRIcomplete(A_tilde):
    
    A_tilde = np.copy(A_tilde)
    n,m = np.shape(A_tilde)
    
    c0 = A_tilde[1:,0]
    c1 = A_tilde[:,1]
    c2 = A_tilde[0:n-1,2]  
    
    A = np.diag(c0, -1) + np.diag(c1) + np.diag(c2, 1)
    
    return(A)

def produitTRIvect(A,x) :
   
    A = np.copy(A)
    n,m = np.shape(A)
    p = np.size(x)
    
    if p != n : # vérification de l'existence du produit
        print('Le produit est impossible.')
        return
    
    Ax = np.zeros((p,1))
    
    for i in range (n) :
        for j in range(m):
            if i+j-1 >= 0 and i+j-1 <= n-1 :
                    Ax[i] += A[i,j]*x[i+j-1] # calcul de chaque terme du vecteur
    return Ax


def DecompositionLU(A):  
    
    A = np.copy(A)
    n,m = A.shape
    U = np.array(A,float)
    L = np.eye(n)
    
    for i in range(n):
        if U[i,i] == 0: 
            print('Un pivot est nul')
            return U
        for j in range(i+1,n):
            g = U[j,i]/U[i,i] 
            U[j,:] = U[j,:]-g*U[i,:] 
            L[j,i] = g 
    a = np.diag(L,-1) 
    b = np.diag(U)  
    c = np.diag(U,1) 
    
    B = np.zeros((n,3))
    B[1:,0] = a 
    B[:,1] = b 
    B[:n-1,2] = c 
    return B
 
 
def triLU (A) :

    A = np.copy(A)
    n,m = np.shape(A)
    M = np.array(A,float) 

    for i in range (1,n) :
            g = A[i,0]/M[i-1,1] 
            M[i,0] = g 
            M[i,1] = M[i,1]-g*M[i-1,2]
            
    return M

    
def TriLUResol(M,b):
    
    n = len(M)
    y = np.zeros((n,1)) 
    x = np.zeros((n,1))
    y[0] = b[0] 
    for i in range(1,n):
        y[i] = b[i] - M[i,0]*y[i-1] 
    x[n-1] = y[n-1]/M[n-1,1] 
    for i in range(2,n+1):
        x[n-i] = (y[n-i]-x[n-i+1]*M[n-i,2]) / M[n-i,1]
    return x


def TriResol(A,b):
    
    M = triLU(A)
    x = TriLUResol(M,b) 
    return x







    
#---------------------test---------------------

# création d'un système tri diagonale
nb = 10000 #nombre de variable
Ar = np.random.random(size=(nb,3)) # matrice tridiagonale reduite
b = np.random.random(size=(nb,))

x=TriResol(Ar, b)

print('la solution trouver avec la méthode de forme reduite :' )
print(x)

print (b)
print(produitTRIvect(Ar, x)) # on utilise la fonction produit trivect pour multiplier 
# notre matrice reduite au solution trouver pour retrouver b