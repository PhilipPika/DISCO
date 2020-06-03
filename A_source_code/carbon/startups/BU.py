import os
import pickle

with open('/Users/pippo/Desktop/DISCO_Test/A_source_code/carbon/startups/start1951.000.pkl','rb') as f:
    argin = pickle.load(f)
print(argin)
argin[0]=1950.0
print(argin)
with open('/Users/pippo/Desktop/DISCO_Test/A_source_code/carbon/startups/start1951.000.pkl','wb') as f:
    pickle.dump(argin,f)

print("Hello World 2!")

