import os
import pickle
print(os.path.dirname(os.path.abspath(__file__)))
print(os.path.abspath(os.getcwd()))
with open('start1951.000.pkl','rb') as f:
    argin = pickle.load(f)
print(argin)
argin[0]=1950.0
print(argin)
with open('start1951.000.pkl','wb') as f:
    pickle.dump(argin,f)

print("Hello World 2!")

