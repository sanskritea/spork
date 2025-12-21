import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_json('C://Users/awschlab/Desktop/data/250803/Pulsed ODMR/82G_300nm_Pulsed ODMR231406final.json') # insert your file location here

freqs = np.array(data['datasets']['freqs']) # you'll need to make sure this matches the frequency columni n Jonathan's data, it might 'frequency' or 'freq' etc, same for the MW ON/OFF columns below

MW_ON = []
MW_OFF = []

for f in freqs:
	# average over all the data points taken for this particular frequency and add to the correct MW ON/OFF array
    MW_ON.append(np.mean(data['datasets']['mwCountsDict'][str(f)])) 
    MW_OFF.append(np.mean(data['datasets']['noMwCountsDict'][str(f)]))
    
diff = np.array(MW_ON) - np.array(MW_OFF)
    
plt.plot(freqs, MW_ON)
plt.plot(freqs, MW_OFF)
# plt.plot(freqs, diff)
plt.show()
