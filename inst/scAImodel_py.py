
# python implementation of scAI optimization model
# import packages
import numpy as np
from numpy import linalg as LA
from numpy import matrix


def scAImodel_py(X1,X2,K,S=0.25,Alpha=1,Lambda=10000,Gamma=1,Maxiter = 500,Stop_rule=1,Seed=1,W1=None,W2=None,H=None,Z=None,R=None):
# parameters
    d1 = X1.shape
    d2 = X2.shape
    p = d1[0]
    n = d1[1]
    q = d2[0]
    np.random.seed(Seed)
    if W1 is None:
        W1 = np.random.rand(p,K)
        
    if W2 is None:
        W2 = np.random.rand(q,K)
    
    if H is None:
        H = np.random.rand(K,n)
    
    if Z is None:
        Z = np.random.rand(n,n)
    
    if R is None:
        R = np.random.binomial(1,S,size=(n,n))
    # main function
    XtX2 = np.dot(np.transpose(X2),X2)
    obj_old = 1
#    W1 = np.asarray(W1)
#    W2 = np.asarray(W2)
#    Z = np.asarray(Z)
#    R = np.asarray(R)
    for iter in range(1,Maxiter+1):
        #  print(iter)
        # normalize H
#        H = matrix(H)
#        lib = np.sum(H,axis = 1)
#        H = H/np.tile(lib,(1,n))
#        H = np.asarray(H)
        H = H/H.sum(axis=1,keepdims=1)
        # update W1
        HHt = np.dot(H,np.transpose(H))
        X1Ht = np.dot(X1,np.transpose(H))
        W1HHt = np.dot(W1,HHt)
        W1 = W1*X1Ht/(W1HHt+np.spacing(1))
        # update W2
        ZR = Z*R
        ZRHt = np.dot(ZR,np.transpose(H))
        X2ZRHt = np.dot(X2,ZRHt)
        W2HHt = np.dot(W2,HHt)
        W2 = W2*X2ZRHt/(W2HHt+np.spacing(1))
        # update H
        W1tX1 = np.dot(np.transpose(W1),X1)
        W2tX2 = np.dot(np.transpose(W2),X2)
        W2tX2ZR = np.dot(W2tX2,ZR)
        HZZt = np.dot(H,(Z+np.transpose(Z)))
        W1tW1 = np.dot(np.transpose(W1),W1)
        W2tW2 = np.dot(np.transpose(W2),W2)
        H = H*(Alpha*W1tX1+W2tX2ZR+Lambda*HZZt)/(np.dot(Alpha*W1tW1+W2tW2+2*Lambda*HHt+Gamma*np.ones([K,K]),H)+np.spacing(1))
        # update Z
        HtH = np.dot(np.transpose(H),H)
        X2tW2H = np.dot(np.transpose(W2tX2),H)
        RX2tW2H = R*X2tW2H
        XtX2ZR = np.dot(XtX2,ZR)
        XtX2ZRR = XtX2ZR*R
        Z = Z*(RX2tW2H+Lambda*HtH)/(XtX2ZRR+Lambda*Z+np.spacing(1))
        if Stop_rule == 2:
            obj = Alpha*pow(LA.norm(X1-np.dot(W1,H),ord = 'fro'),2)+pow(LA.norm(np.dot(X2,ZR)-np.dot(W2,H),ord = 'fro'),2)+Lambda*pow(LA.norm(Z-np.dot(np.transpose(H),H),ord = 'fro'),2)+Gamma*pow(LA.norm(np.dot(np.ones([1,K]),H),ord = 'fro'),2)
            if (obj_old-obj)/obj_old < 1e-6 and iter > 1:
                break
        iter = iter+1
            
# print ("\n")
#   print ("## Running scAI with seed %d ##" % Seed)
#   print ("\n")
    
    return W1, W2, H, Z, R
        
