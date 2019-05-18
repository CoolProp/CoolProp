
import CoolProp
import CoolProp.CoolProp as CP

fluids = [
  ["R50", "Methane"],
  ["R170", "Ethane"],
  ["R600", "n-Butane"],
  ["R704", "Helium"],
  ["R718", "Water"],
  ["R720", "Neon"],
  ["R728", "Nitrogen"],
  ["R729", "Air"],
  ["R732", "Oxygen"],
  ["R740", "Argon"],
  ["R1150", "Ethylene"],
  ["R1270", "Propylene"]]

for f in fluids:
    print("{0}:{1} - {2}:{3}".format(f[0], CP.PropsSI('Tcrit', 'T', 0, 'D', 0, f[0]), CP.PropsSI('Tcrit', 'T', 0, 'D', 0, f[1]), f[1]))
