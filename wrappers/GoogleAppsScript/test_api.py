import requests

# local url
# base_url = "http://127.0.0.1:5000"

base_url = "https://qnikil7.pythonanywhere.com"

# Test PropsSI
function = 'PropsSI'
output = 'D'
property_1 = 'T'
value_1 = 298.15
property_2 = 'P'
value_2 = 100e5
fluid = 'CO2'

string = base_url + '/' + function + '/' + output + '/' + property_1 + '/' + str(value_1) + '/' + property_2 + '/' + str(value_2) + '/' + fluid

response = requests.get(string)

print(response.json())

"""
Density of carbon dioxide at 100 bar and 25C
PropsSI('D', 'T', 298.15, 'P', 100e5, 'CO2')
Out[3]: 817.6273812375758
"""

# Test HAPropsSI
function = 'HAPropsSI'
output = 'H'
property_1 = 'T'
value_1 = 298.15
property_2 = 'P'
value_2 = 101325
property_3 = 'R'
value_3 = 0.5

string = base_url + '/' + function + '/' + output + '/' + property_1 + '/' + str(value_1) + '/' + property_2 + '/' + str(value_2) + '/' + property_3 + '/' + str(value_3)

response = requests.get(string)

print(response.json())

"""
Enthalpy (J per kg dry air) as a function of temperature, pressure,
   and relative humidity at STP
h = HAPropsSI('H','T',298.15,'P',101325,'R',0.5); print(h)
Out[4]: 50423.45039107799
"""
