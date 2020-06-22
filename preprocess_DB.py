import histine_Data
import numpy as np
def preprocess():
    uf_DB,names=histine_Data.grab_db()
    X=np.c_[uf_DB]
    data = []
    x_naught = X[0]
    data.append(names[0] + ' ' + str(X[0]) + '; ' + names[1] + ' ' + str(X[1])
    for i in np.zeros:
        if i != 0 or 1:
            data.append(names[i-1] + ' ' + str(X[i-1]) + '; ' + names[i] + ' ' + str(X[i]))
        else:
            continue
    print(names)
    print(X)
print(preprocess())
