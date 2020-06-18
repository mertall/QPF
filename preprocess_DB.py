import histine_Data
import numpy as np
def preprocess():
    uf_DB,names=histine_Data.grab_db()
    X=np.c_[uf_DB]
    data = []
    #for i in range(len(names)):
        #if i != len(names):
            #data.append(names[i] + ' ' + str(X[i]) + '; ' + names[i+1] + ' ' + str(X[i+1]))
        #else:
            #break
    print(names)
    print(X)
print(preprocess())