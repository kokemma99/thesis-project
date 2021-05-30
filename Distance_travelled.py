import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# polyester dataframe, this is also done for rayon and cotton. Only the names and the input file was changed
df_concentration = pd.read_csv('PE_Conc_Time.csv', index_col = 0 )
df_concentration.rename(columns={df_concentration.columns[0]: "Time" }, inplace = True)
df_concentration.drop(columns=['Time'])
columns = list(df_concentration)

x = 0  # keeping track of how many groups we have had of the 80
sumofconcentrations = 0 # variable for the sum of those 80 groups 
# iterating over the columns
for i, point in enumerate(columns):
    x = x + 1

    if x < 81:
        # 17280 is the last row of the entire csv file
        sumofconcentrations =  sumofconcentrations + df_concentration[point][17280]
        if x == 80:
            if i <= 79:
                # creating a name to be the key of the dictionary 
                name = str(i)
                data = {name: sumofconcentrations}
            else:
                name = str(i)
                data[name] = sumofconcentrations
            # back to 0, otherways all the microfibers will be summed up
            sumofconcentrations = 0
            x = 0

# plotting the data
x = list(range(0, 21))
y = list(data.values())
plt.plot(x,y, '-g')

# naming the labels and the title
plt.xlabel('Distance to the source (in km)')
plt.ylabel('Concentration of microfibers (in kg/m^3)')
plt.title('Distance travelled by polyester microfibers after a time period of two years')
plt.savefig('figure_polyester.pdf')
plt.show()
